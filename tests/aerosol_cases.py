"""Two composed aerosol models, with the reference variables and closed forms."""
import math
from api import *
from harness import case, equal, near, real_equal
from external_helpers import custom, transfers
from constraint_cases import jac_value

QUUX = ['QUUX.BAZ','QUUX.QUX']; CORGE = ['CORGE.FO2','CORGE.BAZ','CORGE.QUX']
STUB1 = ['STUB1.MODE1.'+n for n in QUUX]+['STUB1.MODE2.'+n for n in QUUX+CORGE]
STUB2 = ['STUB2.MODE'+str(mode)+'.'+n for mode,names in [(1,['NUMBER']+QUUX),(2,['NUMBER']+CORGE),(3,['NUMBER']+QUUX+CORGE)] for n in names]
PARAM = 'STUB2.PARAM.MODE2.CORGE.FO2_TO_BAZ_RATE_CONSTANT'
TEMP_PARAM = 'STUB2.PARAM.MODE3.QUUX.BAZ_TO_QUX_RATE_CONSTANT'
RXNS = [('FO2','STUB1.MODE2.CORGE.FO2',.001),('STUB1.MODE1.QUUX.BAZ','STUB1.MODE2.QUUX.BAZ',.002)]


def aerosol():
    m = Model(['FO2','BAR']+STUB1+STUB2)
    k1,k2 = custom(m,PARAM),custom(m,TEMP_PARAM)
    model1 = transfers(m,RXNS)
    def update(sp,context):
        s = sp.contents
        for c in range(s.cells): s.custom[k2*s.stride+c] = .005*s.conditions[c].temperature
    transfers(m,[('STUB2.MODE2.CORGE.FO2','STUB2.MODE2.CORGE.BAZ',k1),('STUB2.MODE3.QUUX.BAZ','STUB2.MODE3.QUUX.QUX',k2)],update)
    m.finalize()
    return m,model1


NAMES = ('StateIncludesStubAerosolModel','CanUpdateStateWithStubAerosolModel','CanUpdateMultiCellStateWithStubAerosolModel',
         'CanCalculateSingleGridCellForcingWithStubAerosolModel','CanCalculateSingleGridCellJacobianWithStubAerosolModel',
         'CanSolveSingleGridCellWithStubAerosolModel1','CanSolveSingleGridCellWithTwoStubAerosolModels',
         'CanSolveMultiGridCellWithStubAerosolModel1','CanSolveMultiGridCellWithTwoStubAerosolModels')


@case(*[family+':AerosolModelIntegration.'+name for family in ('aerosol_model_backward_euler','aerosol_model_rosenbrock') for name in NAMES])
def aerosol_cases(identifier):
    name = identifier.split('.')[-1]; backward = identifier.startswith('aerosol_model_backward')
    multi = 'Multi' in name; both = 'Two' in name; cells = 3 if multi else 1
    methods = (0,) if backward else (2,3,4,5,6)
    if 'Solve' not in name: methods = (3,)
    for method in methods:
        for doolittle in ((False,True) if 'Solve' in name else (False,)):
            m,model1 = aerosol(); s = State(m,cells)
            try:
                if name == 'StateIncludesStubAerosolModel':
                    equal(m.data.nvariables,2+len(STUB1)+len(STUB2)); equal(len(STUB1),7); equal(len(STUB2),13); equal(m.data.ncustom,2)
                    for n in ['FO2','BAR']+STUB1+STUB2: assert m.var(n) != MAX
                elif 'Update' in name:
                    values = {'FO2':[1.34,1.35,1.36],'BAR':[2.53,2.54,2.55],'STUB1.MODE1.QUUX.BAZ':[.75,.76,.77],
                      'STUB1.MODE2.CORGE.FO2':[1.23,1.24,1.25],'STUB2.MODE3.CORGE.QUX':[.42,.43,.44],
                      'STUB2.MODE3.CORGE.BAZ':[.33,.34,.35],'STUB2.MODE1.NUMBER':[1000,1001,1002],'STUB2.MODE2.NUMBER':[500,501,502]}
                    for n,x in values.items(): s.set(n,x[:cells])
                    for n,x in values.items():
                        for c in range(cells):
                            real_equal(s.get(n,c),x[c]); real_equal(s.data.variables[m.var(n)*cells+c],x[c])
                else:
                    initial = {'FO2':[1,1.5,2],'BAR':[2,2.5,3],'STUB1.MODE1.QUUX.BAZ':[.5,.7,.9],
                      'STUB1.MODE2.CORGE.FO2':[.8,1,1.2],'STUB2.MODE3.CORGE.QUX':[.3,.4,.5],
                      'STUB2.MODE3.CORGE.BAZ':[.2,.3,.4],'STUB2.MODE1.NUMBER':[1000,1100,1200],'STUB2.MODE2.NUMBER':[500,550,600]}
                    if both:
                        initial.update({'STUB2.MODE2.CORGE.FO2':[.4,.6,.8],'STUB2.MODE2.CORGE.BAZ':[.1,.15,.2],
                          'STUB2.MODE3.QUUX.BAZ':[.2,.3,.4],'STUB2.MODE3.QUUX.QUX':[.3,.4,.5]})
                        initial['STUB2.MODE3.CORGE.BAZ'] = [.2,.25,.3]
                        # The two-model single-cell fixture initializes QUUX.QUX, not CORGE.QUX.
                        if not multi: initial.pop('STUB2.MODE3.CORGE.QUX')
                    for n,x in initial.items(): s.set(n,x[:cells])
                    if name.startswith('CanCalculate'):
                        f,j = arr([0]*(m.data.nvariables*cells)),arr([0]*(m.data.nnz*cells))
                        if 'Forcing' in name: model1[0](m.p,s.p,s.data.variables,f,None)
                        else: model1[1](m.p,s.p,s.data.variables,j,None)
                        for r,p,k in RXNS:
                            ri,pi = m.var(r),m.var(p)
                            if 'Forcing' in name:
                                real_equal(f[ri],-k*s.get(r)); real_equal(f[pi],k*s.get(r))
                            else:
                                real_equal(jac_value(m,s,j,ri,ri),k); real_equal(jac_value(m,s,j,pi,ri),-k)
                    else:
                        if both:
                            s.set(PARAM,[.01,.015,.02][:cells],custom=True)
                            for c in range(cells): s.data.conditions[c].temperature = [275,285,295][c]
                        s.rates(); p = lib.chem_euler_parameters() if backward else lib.chem_solver(method)
                        if doolittle: p.factor = C.cast(lib.chem_lu_doolittle,P)
                        s.advance(10,parameters=p)
                        tolerance = (.1 if both else 2e-4) if backward else (5e-4 if both and method == 2 else 5e-5)
                        expected = {n:x[:cells] for n,x in initial.items()}
                        reactions = list(RXNS)
                        if both: reactions += [('STUB2.MODE2.CORGE.FO2','STUB2.MODE2.CORGE.BAZ',[.01,.015,.02]),('STUB2.MODE3.QUUX.BAZ','STUB2.MODE3.QUUX.QUX',[275*.005,285*.005,295*.005])]
                        for r,product,k in reactions:
                            before = initial[r]; delta = [before[c]*(1-math.exp(-(k[c] if isinstance(k,list) else k)*10)) for c in range(cells)]
                            for c in range(cells):
                                near(s.get(r,c),before[c]-delta[c],tolerance*before[c])
                                near(s.get(product,c),initial.get(product,[0]*cells)[c]+delta[c],tolerance*before[c])
                        for n in ['BAR','STUB2.MODE1.NUMBER','STUB2.MODE2.NUMBER']+([] if both else ['STUB2.MODE3.CORGE.QUX','STUB2.MODE3.CORGE.BAZ']):
                            for c in range(cells): equal(s.get(n,c),initial[n][c])
            finally: s.close(); m.close()
