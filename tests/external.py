"""C callback adapters for the external-model fixtures in the reference tests.

The callbacks deliberately use the public state/model API and raw supplied
buffers, exercising registration, symbolic binding, and solver dispatch.
"""
from ffi import *

Bind = C.CFUNCTYPE(None,P,P)
Update = C.CFUNCTYPE(None,P,P)
Evaluate = C.CFUNCTYPE(None,P,DP,DP,P)


class ExternalFixture:
    def __init__(self,model,processes=(),constraints=(),updates=None):
        # process: (reactant names, net stoichiometry, numeric or named rate)
        # constraint: (row name, coefficients, constant)
        self.model = model
        self.processes,self.constraints,self.updates = list(processes),list(constraints),updates or {}
        pairs = set()
        for reactants,net,_ in self.processes:
            pairs.update((r,c) for r in net for c in reactants)
        for row,coefficients,_ in self.constraints: pairs.update((row,c) for c in coefficients)
        self.pairs = sorted(pairs)
        self.rows = ints([model.ids[r] for r,c in self.pairs]); self.cols = ints([model.ids[c] for r,c in self.pairs])
        self.alg = ints([model.ids[r] for r,_,_ in self.constraints])
        self.bind_cb,self.update_cb = Bind(self.bind),Update(self.update)
        self.force_cb,self.jacobian_cb = Evaluate(self.forcing),Evaluate(self.jacobian)
        self.constraint_cb,self.constraint_jac_cb = Evaluate(self.residual),Evaluate(self.constraint_jacobian)
        self.description = External(len(self.pairs),len(self.constraints),self.rows,self.cols,self.alg,
            C.cast(self.bind_cb,P),C.cast(self.update_cb,P),None,C.cast(self.force_cb,P),C.cast(self.jacobian_cb,P),
            C.cast(self.constraint_cb,P),C.cast(self.constraint_jac_cb,P),None)
        require(model_add_external(model.handle,C.byref(self.description)))
        if not hasattr(model,'callbacks'): model.callbacks = []
        model.callbacks.append(self)

    def bind(self,handle,_):
        self.indices = {n:model_species_index(handle,n.encode()) for n in self.model.ids}
        self.parameters = {model_parameter_name(handle,i).decode():i for i in range(model_parameter_count(handle))}
        self.flat = {(r,c):model_jacobian_index(handle,self.indices[r],self.indices[c]) for r,c in self.pairs}

    def value(self,v,custom,nc,c):
        return custom[self.parameters[v]*nc+c] if isinstance(v,str) else v

    def update(self,state,_):
        nc = state_cells(state); custom = state_parameters(state); conditions = state_conditions(state)
        for name,function in self.updates.items():
            for c in range(nc): custom[self.parameters[name]*nc+c] = function(conditions[c].temperature)

    def forcing(self,state,y,f,_):
        nc = state_cells(state); custom = state_parameters(state)
        for reactants,net,rate in self.processes:
            for c in range(nc):
                v = self.value(rate,custom,nc,c)
                for n in reactants: v *= y[self.indices[n]*nc+c]
                for n,coeff in net.items(): f[self.indices[n]*nc+c] += coeff*v

    def jacobian(self,state,y,j,_):
        nc = state_cells(state); custom = state_parameters(state)
        for reactants,net,rate in self.processes:
            for q,independent in enumerate(reactants):
                for c in range(nc):
                    v = self.value(rate,custom,nc,c)
                    for i,n in enumerate(reactants):
                        if i != q: v *= y[self.indices[n]*nc+c]
                    for n,coeff in net.items(): j[self.flat[n,independent]*nc+c] -= coeff*v

    def residual(self,state,y,f,_):
        nc = state_cells(state); custom = state_parameters(state)
        for row,terms,constant in self.constraints:
            for c in range(nc):
                v = -constant
                for n,coeff in terms.items(): v += self.value(coeff,custom,nc,c)*y[self.indices[n]*nc+c]
                f[self.indices[row]*nc+c] = v

    def constraint_jacobian(self,state,y,j,_):
        nc = state_cells(state); custom = state_parameters(state)
        for row,terms,_ in self.constraints:
            for n,coeff in terms.items():
                for c in range(nc): j[self.flat[row,n]*nc+c] -= self.value(coeff,custom,nc,c)
