"""Two external aerosol models, with named modes and temperature-dependent rates."""
from ffi import *
from harness import group
from external import ExternalFixture

PFO = 'STUB2.PARAM.MODE2.CORGE.FO2_TO_BAZ_RATE_CONSTANT'
PBAZ = 'STUB2.PARAM.MODE3.QUUX.BAZ_TO_QUX_RATE_CONSTANT'
NAMES = '''StateIncludesStubAerosolModel CanUpdateStateWithStubAerosolModel CanUpdateMultiCellStateWithStubAerosolModel
CanCalculateSingleGridCellForcingWithStubAerosolModel CanCalculateSingleGridCellJacobianWithStubAerosolModel
CanSolveSingleGridCellWithStubAerosolModel1 CanSolveSingleGridCellWithTwoStubAerosolModels
CanSolveMultiGridCellWithStubAerosolModel1 CanSolveMultiGridCellWithTwoStubAerosolModels'''


def aerosol_model():
    m = Model(['FO2','BAR'])
    phase = {'QUUX':('BAZ','QUX'),'CORGE':('FO2','BAZ','QUX')}
    for stub,modes in [(1,{1:['QUUX'],2:['QUUX','CORGE']}),(2,{1:['QUUX'],2:['CORGE'],3:['QUUX','CORGE']})]:
        names = []
        for mode,phases in modes.items():
            if stub == 2: names.append(f'STUB2.MODE{mode}.NUMBER')
            names.extend(f'STUB{stub}.MODE{mode}.{p}.{s}' for p in phases for s in phase[p])
        for n in sorted(names): m.species(n)
    for p in (PFO,PBAZ): m.parameter(p)
    transfers = [('FO2','STUB1.MODE2.CORGE.FO2',.001),('STUB1.MODE1.QUUX.BAZ','STUB1.MODE2.QUUX.BAZ',.002),
                 ('STUB2.MODE2.CORGE.FO2','STUB2.MODE2.CORGE.BAZ',PFO),('STUB2.MODE3.QUUX.BAZ','STUB2.MODE3.QUUX.QUX',PBAZ)]
    for i in (0,2):
        ExternalFixture(m,[([src],{src:-1,dst:1},k) for src,dst,k in transfers[i:i+2]],
                        updates={PBAZ:lambda t:.005*t} if i else {})
    return m.finalize(reorder=True),transfers


@group('aerosol_model_backward_euler:AerosolModelIntegration.',NAMES)
@group('aerosol_model_rosenbrock:AerosolModelIntegration.',NAMES)
def aerosol(id_):
    name = id_.split('.')[-1]; nc = 3 if 'Multi' in name else 1
    be = id_.startswith('aerosol_model_backward'); two = 'TwoStub' in name
    m,transfers = aerosol_model()
    with m:
        if name == 'StateIncludesStubAerosolModel':
            assert len(m.names) == 22 and model_parameter_count(m.handle) == 2
            assert set(m.names) == set(m.ids)
            return
        if 'Update' in name:
            vals = {'FO2':[1.34,1.35,1.36],'BAR':[2.53,2.54,2.55],
                'STUB1.MODE1.QUUX.BAZ':[.75,.76,.77],'STUB1.MODE2.CORGE.FO2':[1.23,1.24,1.25],
                'STUB2.MODE3.CORGE.QUX':[.42,.43,.44],'STUB2.MODE3.CORGE.BAZ':[.33,.34,.35],
                'STUB2.MODE1.NUMBER':[1000,1001,1002],'STUB2.MODE2.NUMBER':[500,501,502]}
            with State(m,nc) as s:
                for n,v in vals.items(): s.set(n,v[:nc])
                for n,v in vals.items(): near(s.y[m.indices[n]],v[:nc],rtol=0)
            return
        values = {'FO2':[1,1.5,2],'BAR':[2,2.5,3],'STUB1.MODE1.QUUX.BAZ':[.5,.7,.9],
            'STUB1.MODE2.CORGE.FO2':[.8,1,1.2],'STUB2.MODE3.CORGE.QUX':[.3,.4,.5],
            'STUB2.MODE3.CORGE.BAZ':[.2,.3,.4],'STUB2.MODE1.NUMBER':[1000,1100,1200],
            'STUB2.MODE2.NUMBER':[500,550,600]}
        if two:
            values.update({'STUB2.MODE2.CORGE.FO2':[.4,.6,.8],'STUB2.MODE2.CORGE.BAZ':[.1,.15,.2],
                'STUB2.MODE3.QUUX.BAZ':[.2,.3,.4],'STUB2.MODE3.QUUX.QUX':[.3,.4,.5],
                'STUB2.MODE3.CORGE.BAZ':[.2,.25,.3]})
        methods = ('be',) if be else (2,3,4,43,6)
        for method in methods if 'Solve' in name else (3,):
            for lu in (0,1):
                with State(m,nc) as s:
                    require(state_set_lu(s.handle,lu))
                    for n,v in values.items(): s.set(n,v[:nc])
                    if two:
                        s.conditions(np.array([275,285,295])[:nc]); s.parameter(PFO,[.01,.015,.02][:nc])
                    else: s.conditions(0,0,0)
                    s.update(); initial = s.y.copy()
                    if 'Forcing' in name:
                        target = np.zeros_like(s.y)
                        for src,dst,k in transfers[:2]:
                            v = k*s.y[m.indices[src]]; target[m.indices[src]] -= v; target[m.indices[dst]] += v
                        near(s.forcing(),target,rtol=0); continue
                    if 'Jacobian' in name:
                        target = np.zeros((nc,len(m.names),len(m.names)))
                        for src,dst,k in transfers[:2]:
                            target[:,m.indices[src],m.indices[src]] -= k; target[:,m.indices[dst],m.indices[src]] += k
                        near(s.jacobian(),target,rtol=0); continue
                    s.solve(10,method); expected = initial.copy(); tolerances = np.zeros_like(s.y)
                    tol = (.1 if two else 2e-4) if be else (5e-4 if two and method == 2 else 5e-5)
                    for src,dst,k in transfers[:4 if two else 2]:
                        rate = s.custom[model_parameter_index(m.handle,k.encode())] if isinstance(k,str) else k
                        amount = initial[m.indices[src]]*(1-np.exp(-rate*10))
                        expected[m.indices[src]] -= amount; expected[m.indices[dst]] += amount
                        tolerances[m.indices[src]] = tol*initial[m.indices[src]]
                        tolerances[m.indices[dst]] = tol*initial[m.indices[src]]
                    near(s.y,expected,rtol=0,atol=tolerances)
