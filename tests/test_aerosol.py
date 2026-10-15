"""External state registration and aerosol process callback case ports."""
from functools import partial
import json
import math
from pathlib import Path
from builderffi import *
from test_state import named

S1=sorted(['STUB1.MODE1.QUUX.'+s for s in ['BAZ','QUX']]+['STUB1.MODE2.QUUX.'+s for s in ['BAZ','QUX']]+['STUB1.MODE2.CORGE.'+s for s in ['FO2','BAZ','QUX']])
S2=sorted(['STUB2.MODE1.NUMBER','STUB2.MODE2.NUMBER','STUB2.MODE3.NUMBER']+
          ['STUB2.MODE1.QUUX.'+s for s in ['BAZ','QUX']]+['STUB2.MODE2.CORGE.'+s for s in ['FO2','BAZ','QUX']]+
          ['STUB2.MODE3.QUUX.'+s for s in ['BAZ','QUX']]+['STUB2.MODE3.CORGE.'+s for s in ['FO2','BAZ','QUX']])
P2=['STUB2.PARAM.MODE2.CORGE.FO2_TO_BAZ_RATE_CONSTANT','STUB2.PARAM.MODE3.QUUX.BAZ_TO_QUX_RATE_CONSTANT']


class Transfers(Definition):
    def __init__(self,variables,parameters,transfers,temperature_parameter=None):
        self.transfers,self.temperature_parameter=transfers,temperature_parameter
        super().__init__(variables,parameters)
    def configure(self,data,process,constraints):
        jobs=[];coordinates=[];flats={}
        for source,target,rate in self.transfers:
            r,p=lookup(data,source),lookup(data,target)
            parameter=lookup(data,rate,True) if isinstance(rate,str) else None
            jobs.append((r,p,rate,parameter))
            coordinates.extend([(r,r),(p,r)])
        coords=self.array(MatrixIndex,[MatrixIndex(*ij) for ij in sorted(set(coordinates))])
        process.nonzero_count=len(coords);process.nonzeros=coords
        def finish(m,ctx):
            for row,col in coordinates:
                flat=lib.chem_jacobian_index(m,row,col);assert flat!=Z(-1).value
                flats[row,col]=flat
        def forcing(state,y,f,ctx):
            nc=lib.chem_cell_count(state);params=lib.chem_parameters(state)
            for r,p,rate,parameter in jobs:
                for c in range(nc):
                    k=rate if parameter is None else params[parameter*nc+c]
                    value=k*y[r*nc+c];f[r*nc+c]-=value;f[p*nc+c]+=value
        def jacobian(state,y,jac,ctx):
            nc=lib.chem_cell_count(state);params=lib.chem_parameters(state)
            for r,p,rate,parameter in jobs:
                for c in range(nc):
                    k=rate if parameter is None else params[parameter*nc+c]
                    jac[flats[r,r]*nc+c]+=k;jac[flats[p,r]*nc+c]-=k
        parameter=lookup(data,self.temperature_parameter,True) if self.temperature_parameter else None
        def update(state,ctx):
            if parameter is not None:
                nc=lib.chem_cell_count(state);params=lib.chem_parameters(state);temp=lib.chem_temperatures(state)
                for c in range(nc):params[parameter*nc+c]=.005*temp[c]
        self.attach(process,forcing,jacobian,update,finish)


def builder():
    one=Transfers(S1,[],[('FO2','STUB1.MODE2.CORGE.FO2',.001),('STUB1.MODE1.QUUX.BAZ','STUB1.MODE2.QUUX.BAZ',.002)])
    two=Transfers(S2,P2,[('STUB2.MODE2.CORGE.FO2','STUB2.MODE2.CORGE.BAZ',P2[0]),('STUB2.MODE3.QUUX.BAZ','STUB2.MODE3.QUUX.QUX',P2[1])],P2[1])
    return Builder(['FO2','BAR'],models=[one,two])


def metadata():
    b=builder()
    with b.build() as m:
        assert m.n==2+len(S1)+len(S2)==22 and m.data.parameter_count==2
        assert (len(S1),len(S2),len(P2))==(7,13,2)
        for name in ['FO2','BAR']+S1+S2:assert lib.chem_species_index(m.handle,name.encode())!=Z(-1).value
        for name in P2:assert lib.chem_parameter_index(m.handle,name.encode())!=Z(-1).value
        b.check()


def update(multi=False):
    data={'FO2':[1.34,1.35,1.36],'BAR':[2.53,2.54,2.55],'STUB1.MODE1.QUUX.BAZ':[.75,.76,.77],
          'STUB1.MODE2.CORGE.FO2':[1.23,1.24,1.25],'STUB2.MODE3.CORGE.QUX':[.42,.43,.44],
          'STUB2.MODE3.CORGE.BAZ':[.33,.34,.35],'STUB2.MODE1.NUMBER':[1000.,1001.,1002.],
          'STUB2.MODE2.NUMBER':[500.,501.,502.]}
    b=builder()
    with b.build(3 if multi else 1) as m:
        named(m.state,{name:values if multi else values[:1] for name,values in data.items()})
        for name,values in data.items():
            index=lib.chem_species_index(m.handle,name.encode())
            assert index!=Z(-1).value
            expected=values if multi else values[:1]
            assert list(lib.chem_species_column(m.state,name.encode())[:m.nc])==expected
            assert [row[index] for row in m.rows()]==expected
        b.check()


def initial(multi=False,two=False):
    values={'FO2':[1.,1.5,2.],'BAR':[2.,2.5,3.],'STUB1.MODE1.QUUX.BAZ':[.5,.7,.9],
            'STUB1.MODE2.CORGE.FO2':[.8,1.,1.2],'STUB2.MODE3.CORGE.QUX':[.3,.4,.5],
            'STUB2.MODE3.CORGE.BAZ':[.2,.3,.4],'STUB2.MODE1.NUMBER':[1000.,1100.,1200.],
            'STUB2.MODE2.NUMBER':[500.,550.,600.]}
    if two:
        del values['STUB2.MODE3.CORGE.QUX']
        values.update({'STUB2.MODE2.CORGE.FO2':[.4,.6,.8],'STUB2.MODE2.CORGE.BAZ':[.1,.15,.2],
                       'STUB2.MODE3.QUUX.BAZ':[.2,.3,.4],'STUB2.MODE3.QUUX.QUX':[.3,.4,.5],
                       'STUB2.MODE3.CORGE.BAZ':[.2,.25,.3]})
    return {name:array if multi else array[:1] for name,array in values.items()}


def derivative(jacobian=False):
    b=builder()
    with b.build() as m:
        named(m.state,initial())
        process=m.data.external_models[0]
        if jacobian:
            output=(D*lib.chem_jacobian_size(m.handle))()
            MODEL(process.negative_jacobian)(m.state,m.y,output,None)
            for source,target,k in b.models[0].transfers:
                r,p=lookup(C.byref(m.data),source),lookup(C.byref(m.data),target)
                assert output[lib.chem_jacobian_index(m.handle,r,r)]==k
                assert output[lib.chem_jacobian_index(m.handle,p,r)]==-k
        else:
            output=(D*m.n)()
            MODEL(process.forcing)(m.state,m.y,output,None)
            for source,target,k in b.models[0].transfers:
                r,p=lookup(C.byref(m.data),source),lookup(C.byref(m.data),target)
                assert output[r]==-k*m.y[r] and output[p]==k*m.y[r]
        b.check()


def solve(method,group,algorithm,tolerance,multi=False,two=False):
    values=initial(multi,two);nc=3 if multi else 1
    for start in range(0,nc,group or nc):
        count=min(group or nc,nc-start)
        b=builder()
        with b.build(count) as m:
            named(m.state,{name:array[start:start+count] for name,array in values.items()})
            if two:
                for c in range(count):
                    lib.chem_temperatures(m.state)[c]=[275.,285.,295.][start+c]
                named(m.state,{P2[0]:[.01,.015,.02][start:start+count]},parameters=True)
            assert lib.chem_update_rates(m.state)==0
            assert lib.chem_set_linear_solver(m.state,getattr(lib,'chem_linear_'+algorithm)())==0
            if method:status,stats=m.solve(10.,method)
            else:
                p=lib.chem_backward_euler_parameters();stats=Stats()
                status=lib.chem_solve_backward_euler(m.state,C.byref(p),10.,C.byref(stats))
            assert status==0,(method,group,algorithm,status)
            expected={name:list(array[start:start+count]) for name,array in values.items()}
            affected={}
            def transfer(source,target,constants):
                target_values=expected.setdefault(target,[0.]*count)
                for c,k in enumerate(constants):
                    base=values[source][start+c];change=base*(1-math.exp(-k*10.))
                    expected[source][c]-=change;target_values[c]+=change
                affected[source]=affected[target]=values[source][start:start+count]
            transfer('FO2','STUB1.MODE2.CORGE.FO2',[.001]*count)
            transfer('STUB1.MODE1.QUUX.BAZ','STUB1.MODE2.QUUX.BAZ',[.002]*count)
            if two:
                transfer('STUB2.MODE2.CORGE.FO2','STUB2.MODE2.CORGE.BAZ',[.01,.015,.02][start:start+count])
                transfer('STUB2.MODE3.QUUX.BAZ','STUB2.MODE3.QUUX.QUX',[.005*t for t in [275.,285.,295.][start:start+count]])
            for name,array in expected.items():
                actual=list(lib.chem_species_column(m.state,name.encode())[:count])
                for c,(got,want) in enumerate(zip(actual,array)):
                    if name in affected:assert abs(got-want)<=tolerance*affected[name][c],(method,group,algorithm,name,got,want,tolerance)
                    else:assert got==want,(name,got,want)
            b.check()


def check(name,variants):
    if name=='StateIncludesStubAerosolModel':return metadata()
    if name=='CanUpdateStateWithStubAerosolModel':return update()
    if name=='CanUpdateMultiCellStateWithStubAerosolModel':return update(True)
    if 'Forcing' in name:return derivative()
    if 'Jacobian' in name:return derivative(True)
    settings={}
    for row in variants:
        key=(row['method'],row['group'],row['algorithm'])
        settings[key]=min(settings.get(key,math.inf),row['tolerance'])
    for (method,group,algorithm),tol in settings.items():solve(method,group,algorithm,tol,'MultiGridCell' in name,'TwoStub' in name)


CASES={key:partial(check,key.split('.')[-1],variants) for key,variants in
       json.loads((Path(__file__).parent/'data/aerosol_variants.json').read_text()).items()}
