"""Residual, Jacobian, metadata, validation, and finite-difference constraint ports."""
import math
from api import *
from arrays import Array, AP
from verify_cases import ForcingFn
from harness import case, equal, near


def sparsity(m):
    n = m.data.nvariables; pattern = (C.c_ubyte*(n*n))()
    equal(lib.chem_constraint_sparsity(m.p,pattern,n),0)
    return {(i,j) for i in range(n) for j in range(n) if pattern[i*n+j]}


def evaluate(m,rows,initial=0.0):
    s = State(m,len(rows)); n, nc = m.data.nvariables,len(rows)
    for c,row in enumerate(rows):
        for i,v in enumerate(row): s.data.variables[i*nc+c] = v
    forcing = arr([initial]*(n*nc)); jacobian = arr([0]*(m.data.nnz*nc))
    lib.chem_constraints(s.p,s.data.variables,forcing,jacobian)
    return s,forcing,jacobian


def jac_value(m,s,jac,row,col,cell=0):
    at = m.data.index[row*m.data.nvariables+col]
    return jac[at*s.data.stride+cell] if at >= 0 else 0.0


def verify_jacobian(s,constraints_only=False,atol=1e-7,rtol=1e-7):
    m = s.model; n, nc = m.data.nvariables,s.data.cells
    v = Array(nc,n)
    for c in range(nc):
        for i in range(n): v[c,i] = s.data.variables[i*s.data.stride+c]
    errors = []
    def force(vp,fp,ctx):
        try:
            vv,ff = Array(pointer=vp),Array(pointer=fp)
            y = arr([vv[c,i] for i in range(n) for c in range(nc)]); f = arr([0]*(n*nc))
            if constraints_only: lib.chem_constraints(s.p,y,f,None)
            else: lib.chem_forcing(s.p,y,f)
            for c in range(nc):
                for i in range(n): ff[c,i] = f[i*nc+c]
        except Exception as error: errors.append(error)
    callback = ForcingFn(force); output = AP()
    equal(lib.chem_finite_difference(v.p,n,1e-8,callback,None,C.byref(output)),0)
    if errors: raise errors[0]
    fd = Array(pointer=output); fd.owned = True
    pattern = [(i,j) for i in range(n) for j in range(n) if m.data.index[i*n+j] >= 0]
    a = Array(nc,pattern=pattern,order=n)
    jac = arr([0]*(m.data.nnz*s.data.stride))
    if constraints_only: lib.chem_constraints(s.p,s.data.variables,None,jac)
    else: lib.chem_jacobian(s.p,s.data.variables,jac)
    try:
        for c in range(nc):
            for i,j in pattern: a[c,i,j] = jac_value(m,s,jac,i,j,c)
        comparison = lib.chem_compare_jacobian(a.p,fd.p,n,atol,rtol)
        equal(comparison.passed,1)
        equal(lib.chem_check_sparsity(a.p,fd.p,n,1e-6).passed,1)
    finally: a.close(); fd.close(); v.close()


LINEAR_DATA = {
    'Construction': ('A_B_conservation',['A','B'],'B',[1,1],1,[[.3,.7]]),
    'AlgebraicSpecies': ('conservation',['A','B','C'],'C',[1,1,1],10,[[2,3,5]]),
    'WeightedTerms': ('weighted_sum',['A','B','C'],'C',[2,3,-1],5,[[1,2,3]]),
    'ZeroConstant': ('A_equals_B',['A','B'],'B',[1,-1],0,[[.5,.5]]),
    'FractionalCoefficients': ('fractional_conservation',['A','B'],'B',[.5,1.5],2,[[1,1]]),
    'ResidualComputationThroughConstraintSet': ('A_B_conservation',['A','B'],'B',[1,1],1,[[.3,.7],[.5,.6]]),
    'JacobianComputationThroughConstraintSet': ('A_B_conservation',['A','B'],'B',[1,1],1,[[.3,.7]]),
    'WeightedSumResidualAndJacobian': ('weighted_sum',['A','B','C'],'C',[2,3,-1],5,[[1,2,3]]),
    'ThreeSpeciesConservationResidual': ('conservation',['A','B','C'],'C',[1,1,1],10,[[2,3,5],[2,3,6]]),
    'ZeroConstantResidual': ('A_equals_B',['A','B'],'B',[1,-1],0,[[.5,.5],[.6,.5]]),
    'FractionalCoefficientsResidualAndJacobian': ('fractional_conservation',['A','B'],'B',[.5,1.5],2,[[1,1]]),
    'JacobianIndependentOfConcentrations': ('A_B_sum',['A','B'],'B',[2,3],1,[[.1,.2],[10,20]]),
    'FiniteDifferenceJacobianSimpleConservation': ('conservation',['A','B'],'B',[1,1],1,[[.3,.7],[.8,.2]]),
    'FiniteDifferenceJacobianWeightedSum': ('weighted',['A','B','C'],'C',[2,3,-1],5,[[1.5,.8,2]]),
}


@case(*['linear:LinearConstraint.'+name for name in LINEAR_DATA])
def linear(identifier):
    name = identifier.split('.')[-1]; label,names,algebraic,coeffs,constant,rows = LINEAR_DATA[name]
    m = Model(names); m.constraint(algebraic,list(zip(names,coeffs)),constant,name=label); m.finalize()
    s,f,j = evaluate(m,rows)
    try:
        q = m.data.constraints[0]; equal(q.name.decode(),label); equal(q.constant,constant); equal(q.nterms,len(names)); equal(q.algebraic,m.var(algebraic))
        equal(list(q.coefficients[:q.nterms]),coeffs); equal(list(q.species[:q.nterms]),list(range(len(names))))
        for c,row in enumerate(rows):
            expected = sum(x*y for x,y in zip(row,coeffs))-constant
            near(f[q.algebraic*s.data.stride+c],expected,1e-10)
            for col,coef in enumerate(coeffs): near(jac_value(m,s,j,q.algebraic,col,c),-coef,1e-10)
        if name.startswith('FiniteDifference'): verify_jacobian(s,True)
    finally: s.close(); m.close()


EQUILIBRIUM_NAMES = ('Construction','AlgebraicSpecies','SingleReactantSingleProduct','MultipleReactantsAndProducts',
    'InvalidEquilibriumConstant','EmptyReactantsThrows','EmptyProductsThrows','InvalidStoichiometryThrows',
    'ResidualComputationThroughConstraintSet','JacobianComputationThroughConstraintSet','ComplexStoichiometryResidual',
    'FiniteDifferenceJacobianSimple','FiniteDifferenceJacobianComplexStoichiometry')


@case(*['equilibrium:EquilibriumConstraint.'+name for name in EQUILIBRIUM_NAMES])
def equilibrium(identifier):
    name = identifier.split('.')[-1]
    if name in ('InvalidEquilibriumConstant','EmptyReactantsThrows','EmptyProductsThrows','InvalidStoichiometryThrows'):
        m = Model(['A','B'])
        try:
            invalid = []
            if name == 'InvalidEquilibriumConstant': invalid = [(k,1,[1,1]) for k in (-1,0)]
            elif name == 'EmptyReactantsThrows': invalid = [(1,0,[1,1])]
            elif name == 'EmptyProductsThrows': invalid = [(1,2,[1,1])]
            else: invalid = [(1,1,coefs) for coefs in ([0,1],[-1,1],[1,0],[1,-2])]
            for k,nr,coeffs in invalid:
                q = Constraint(1,1,2,nr,inds([0,1]),arr(coeffs),k,None,None,None)
                equal(lib.chem_model_constraint(m.p,C.byref(q)),1)
        finally: m.close()
        return
    simple = name == 'SingleReactantSingleProduct'
    complex_ = name in ('MultipleReactantsAndProducts','ComplexStoichiometryResidual','FiniteDifferenceJacobianComplexStoichiometry')
    names = ['A','B'] if simple else ['A','B','C' if complex_ else 'AB']
    coeffs = [1,1] if simple else [2,1,1] if complex_ else [1,1,1]
    nr,k,row = (1,10,1) if simple else (1,100,1) if complex_ else (2,1000,2)
    label = 'A_B_simple' if simple else 'dissociation' if complex_ else 'A_B_equilibrium'
    m = Model(names); m.constraint(row,list(zip(names,coeffs)),k,equilibrium_reactants=nr,name=label); m.finalize()
    rows = [[.1,.5,2]] if complex_ else [[.001,.001,.001],[.02,.03,.05]] if not simple else [[1,10]]
    if name == 'JacobianComputationThroughConstraintSet': rows = [[.01,.02,.05]]
    if name == 'FiniteDifferenceJacobianSimple': rows = [[.02,.03,.05],[.1,.2,5]]
    if name == 'FiniteDifferenceJacobianComplexStoichiometry': rows = [[.15,.8,1.5]]
    s,f,j = evaluate(m,rows)
    try:
        q = m.data.constraints[0]; equal(q.name.decode(),label); equal(q.algebraic,row); equal(q.nreactants,nr); equal(q.nterms,len(names))
        equal(list(q.species[:len(names)]),list(range(len(names)))); equal(list(q.coefficients[:len(names)]),coeffs)
        for c,values in enumerate(rows):
            expected = k*math.prod(max(x,0)**stoich for x,stoich in zip(values[:nr],coeffs[:nr]))-math.prod(max(x,0)**stoich for x,stoich in zip(values[nr:],coeffs[nr:]))
            near(f[row*s.data.stride+c],expected,1e-10)
            for col in range(len(names)):
                lo,hi = (0,nr) if col < nr else (nr,len(names))
                derivative = (k if col < nr else -1)*coeffs[col]*values[col]**(coeffs[col]-1)
                for t in range(lo,hi):
                    if t != col: derivative *= values[t]**coeffs[t]
                near(jac_value(m,s,j,row,col,c),-derivative,1e-8)
        if name.startswith('FiniteDifference'): verify_jacobian(s,True)
    finally: s.close(); m.close()


SET_NAMES = ('Construction','ReplaceStateRowsMapsToAlgebraicSpecies','NonZeroJacobianElements','MultipleConstraints',
    'AddForcingTerms','SubtractJacobianTerms','EmptyConstraintSet','UnknownSpeciesThrows','ThreeDStateOneConstraint',
    'FourDStateTwoConstraints','CoupledConstraintsSharedSpecies','VectorizedMatricesRespectGridCellIndexing','VectorMatrix1','VectorMatrix2','VectorMatrix3')


@case(*['constraint_set:ConstraintSet.'+name for name in SET_NAMES])
def constraint_set(identifier):
    name = identifier.split('.')[-1]
    if name.startswith('VectorMatrix'):
        for subcase in SET_NAMES:
            if not subcase.startswith('VectorMatrix'): constraint_set('constraint_set:ConstraintSet.'+subcase)
        return
    if name == 'UnknownSpeciesThrows':
        m = Model(['A','B'])
        try:
            q = Constraint(1,1,3,2,inds([0,1,2]),arr([1,1,1]),.033,None,None,None)
            equal(lib.chem_model_constraint(m.p,C.byref(q)),3)
        finally: m.close()
        return
    if name == 'EmptyConstraintSet':
        m = Model(['A','B','C']).finalize(); s = State(m)
        try:
            equal(m.data.nconstraints,0); equal(sparsity(m),set())
            f = arr([1]*4); lib.chem_constraints(s.p,s.data.variables,f,None); equal(list(f),[1]*4)
        finally: s.close(); m.close()
        return
    specs = [('AB',[('A',1),('B',1),('AB',1)],2)]
    names,rows = ['A','B','AB'],[[.2,.4,.6],[.1,.3,.7]]
    if name == 'ReplaceStateRowsMapsToAlgebraicSpecies':
        names,rows,specs = ['A','B','C'],[[0,.2,.3]],[('C',[('B',1),('C',1)],1)]
    elif name == 'MultipleConstraints':
        names,rows,specs = ['A','B','AB','C','D'],[[.2,.4,.6,.1,.3]],specs+[('D',[('C',1),('D',1)],1)]
    elif name == 'ThreeDStateOneConstraint':
        names,rows,specs = ['X','Y','Z'],[[10,.2,.5],[10,.33,.3]],[('Y',[('X',1),('Y',1)],1)]
    elif name == 'FourDStateTwoConstraints':
        names,rows,specs = ['A','B','C','D'],[[.33,.01089,10,1],[.33,.01089,5,1],[.33,.02,5,1]],[('B',[('A',1),('B',1)],1),('A',[('C',1),('D',1),('A',1)],2)]
    elif name == 'CoupledConstraintsSharedSpecies':
        names,rows,specs = ['A','B','C'],[[.1,.0033,.0033]],[('B',[('A',1),('B',1)],1),('C',[('A',1),('C',1)],1)]
    elif name == 'VectorizedMatricesRespectGridCellIndexing': rows = [[.01,.02,.05],[.03,.01,.2],[.001,.002,.004]]
    elif name == 'SubtractJacobianTerms': rows = [[.01,.02,.05]]
    m = Model(names)
    for index,(alg,terms,nr) in enumerate(specs): m.constraint(alg,terms,.033,equilibrium_reactants=nr,name=f'eq{index}')
    m.finalize(); s,f,j = evaluate(m,rows)
    try:
        equal(m.data.nconstraints,len(specs))
        expected_pattern = {(m.var(alg),m.var(dep)) for alg,terms,nr in specs for dep,coef in terms}
        equal(sparsity(m),expected_pattern)
        for c,values in enumerate(rows):
            for alg,terms,nr in specs:
                xs = [values[m.var(dep)] for dep,coef in terms]
                expected = .033*math.prod(xs[:nr])-math.prod(xs[nr:])
                near(f[m.var(alg)*s.data.stride+c],expected,1e-10)
                for i,(dep,coef) in enumerate(terms):
                    lo,hi = (0,nr) if i < nr else (nr,len(terms))
                    derivative = .033 if i < nr else -1
                    for t in range(lo,hi):
                        if t != i: derivative *= xs[t]
                    near(jac_value(m,s,j,m.var(alg),m.var(dep),c),-derivative,1e-12)
        for i in range(len(names)):
            if not m.data.algebraic[i]:
                for c in range(len(rows)): equal(f[i*s.data.stride+c],0); equal(jac_value(m,s,j,i,i,c),0)
    finally: s.close(); m.close()
