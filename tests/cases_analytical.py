"""Analytical integration fixtures, including all five Rosenbrock formulas.

Values and error bounds follow MICM's analytical_policy.hpp. The oracle for
each rate is written independently of the C rate implementation.
"""
import math
from ffi import *
from harness import group

METHODS = (2,3,4,43,6)
RATE_KINDS = {'Troe':1,'Photolysis':7,'TernaryChemicalActivation':2,
              'Tunneling':4,'Arrhenius':0,'Branched':5}
NAMES = ' '.join(n+s for n in RATE_KINDS for s in ('','SuperStiffButAnalytical'))
RGAS = 1.380649e-23*6.02214076e23


def advance(state,dt,method=3,substeps=1):
    """Resume a solver at its reported time after a maximum-step return."""
    for _ in range(substeps):
        elapsed = 0.
        while elapsed < dt/substeps:
            status,stats = state.solve(dt/substeps-elapsed,method,allow_failure=True)
            assert status in (0,-7), status_string(status)
            assert stats.final_time > 0
            elapsed += stats.final_time


def add_user(m,reactants,products,value,label=None):
    label = label or 'k'+str(model_reaction_count(m.handle))
    rate = rate_default(7); rate.custom_index = m.parameter(label,value)
    m.reaction(reactants,products,rate)


def chain_rate(name,second,stiff=False):
    rate = rate_default(RATE_KINDS[name])
    if name in ('Troe','TernaryChemicalActivation'):
        rate.p[0] = 4e-11 if name == 'Troe' else 4e-5
        if second: rate.p[:8] = (1.2e-3,1.6,3.,136.,5.,24.,.9,.8)
    elif name == 'Tunneling':
        rate.p[:3] = (1.2e-4,167.,1e8) if second else (.004,0.,0.)
    elif name == 'Arrhenius':
        rate.p[:5] = (1.2e-4,1.6 if stiff else 7.,75.,50.,.5) if second else (.004,0.,50.,300.,0.)
    elif name == 'Branched':
        rate.p[:5] = (1e-4,204.3,1e-3,2.,int(second))
    rate_prepare(C.byref(rate))
    return rate


def chain_oracle(name,second,t,p,air,stiff):
    if name == 'Photolysis': return np.full(3,.003 if second else .002)
    rate = chain_rate(name,second,stiff); a = rate.p
    if name in ('Troe','TernaryChemicalActivation'):
        k0 = a[0]*np.exp(a[2]/t)*(t/300.)**a[1]
        ki = a[3]*np.exp(a[5]/t)*(t/300.)**a[4]
        ratio = k0*air/ki
        return k0*(air if name == 'Troe' else 1)/(1+ratio)*a[6]**(a[7]/(a[7]+np.log10(ratio)**2))
    if name == 'Tunneling': return a[0]*np.exp(-a[1]/t+a[2]/t**3)
    if name == 'Arrhenius': return a[0]*np.exp(a[2]/t)*(t/a[3])**a[1]*(1+a[4]*p)
    ar = 2e-22*math.exp(2)*2.45e19; br = .43*(293/298)**-8
    z = ar/(1+ar/br)*.41**(1/(1+math.log10(ar/br)**2))*(1-.001)/.001
    av = 2e-22*math.exp(2)*air*6.02214076e23*1e-6; bv = .43*(t/298)**-8
    act = av/(1+av/bv)*.41**(1/(1+np.log10(av/bv)**2))
    return 1e-4*np.exp(-204.3/t)*(act if second else z)/(z+act)


def chain_tolerance(name,stiff,method):
    if method == 'be':
        return {'Troe':(1e-6,1e-5),'Photolysis':(1e-3,1e-3),
                'TernaryChemicalActivation':(1e-5,1e-2),'Tunneling':(1e-3,1e-3),
                'Arrhenius':(1e-3,1e-3),'Branched':(1e-5,1e-2)}[name][stiff]
    if stiff:
        if name == 'Arrhenius': return (1e-4,2e-5,2e-5,2e-5,1e-5)[METHODS.index(method)]
        return {'Troe':1e-5,'Photolysis':2e-5,'TernaryChemicalActivation':2e-3,
                'Tunneling':1e-4,'Branched':2e-3}[name]
    return {'Troe':1e-10,'Photolysis':2e-6,'TernaryChemicalActivation':1e-8,
            'Tunneling':2e-5 if method == 2 else 1e-8,
            'Arrhenius':4e-6 if method == 2 else 1e-9,
            'Branched':1e-10 if method == 2 else 1e-13}[name]


@group('analytical_rosenbrock_integration:AnalyticalExamples.',NAMES)
@group('analytical_backward_euler:AnalyticalExamples.',NAMES)
def chains(id_):
    name = id_.split('.')[-1]; stiff = name.endswith('SuperStiffButAnalytical')
    name = name.removesuffix('SuperStiffButAnalytical')
    methods = ('be',) if id_.startswith('analytical_backward') else METHODS
    temperature = np.array([272.5,254.7,312.6]); pressure = np.array([101253.3,100672.5,101319.8])
    air = pressure/(RGAS*temperature); initial = 1-np.arange(3)/3
    k1 = chain_oracle(name,False,temperature,pressure,air,stiff)
    k2 = chain_oracle(name,True,temperature,pressure,air,stiff)
    time = np.arange(1000)[:,None]
    a = initial*np.exp(-k1*time)
    b = initial*k1/(k2-k1)*(np.exp(-k1*time)-np.exp(-k2*time))
    c = initial*(1+(k1*np.exp(-k2*time)-k2*np.exp(-k1*time))/(k2-k1))
    expected = np.stack((a,b,c),axis=1)
    with Model(['A1','A2','B','C'] if stiff else ['A','B','C']) as m:
        if stiff:
            for src,dst,k in [('A1','A2',4e10),('A2','A1',.9*4e10)]:
                r = rate_default(0); r.p[0] = k; m.reaction([src],{dst:1},r)
        for src,dst,second in ([('A1','B',False),('A2','B',False),('B','C',True)] if stiff else
                               [('A','B',False),('B','C',True)]):
            r = chain_rate(name,second,stiff)
            if name == 'Photolysis': r.custom_index = m.parameter(src+dst,.003 if second else .002)
            m.reaction([src],{dst:1},r)
        m.finalize(reorder=True)
        for method in methods:
            # Exercise both numerical LU algorithms for each parameter set.
            for lu in (0,1):
                with State(m,3) as s:
                    require(state_set_lu(s.handle,lu)); s.conditions(temperature,pressure)
                    if stiff: s.set('A1',initial/2); s.set('A2',initial/2)
                    else: s.set('A',initial)
                    history = np.zeros_like(expected); history[0] = expected[0]
                    for step in range(1,1000):
                        s.update(); s.solve(1.,method)
                        history[step,0] = s.y[m.indices['A1']]+s.y[m.indices['A2']] if stiff else s.y[m.indices['A']]
                        history[step,1] = s.y[m.indices['B']]; history[step,2] = s.y[m.indices['C']]
                    if not stiff: near(s.rates,np.stack((k1,k2)),rtol=1e-9)
                    near(history,expected,rtol=0,atol=chain_tolerance(name,stiff,method))


@group('analytical_rosenbrock_integration:AnalyticalExamples.','SurfaceRxn')
@group('analytical_backward_euler:AnalyticalExamples.','SurfaceRxn')
def surface(id_):
    methods = ('be',) if id_.startswith('analytical_backward') else METHODS
    diameter,gsd = 1e-6,.1
    radius = diameter/2*math.exp(5*math.log(gsd)**2/2)
    number = 6/(math.pi*diameter**3*math.exp(4.5*math.log(gsd)**2))*(.002/1000+.003/1000)
    k = 4*number*math.pi*radius**2/(radius/.95e-5+4/(math.sqrt(8*RGAS*272.5/(math.pi*.04607))*.02))
    with Model() as m:
        m.species('foo',**{'molecular weight [kg mol-1]':.04607,'diffusion coefficient [m2 s-1]':.95e-5})
        m.species('bar'); m.species('baz')
        r = rate_default(6); r.p[2] = .02; r.custom_index = m.parameter('foo.radius',radius)
        m.parameter('foo.number',number); m.reaction(['foo'],{'bar':1,'baz':.4},r); m.finalize()
        for method in methods:
            with State(m) as s:
                s.conditions(272.5,101253.3); s.set('foo',1.)
                tol = .05 if method == 'be' else (1e-2,1e-5,1e-6,1e-5,1e-7)[METHODS.index(method)]
                for step in range(1,11):
                    s.update(); near(s.rates[0,0],k,rtol=0,atol=1e-8)
                    s.solve(.1/k,method); a = math.exp(-step*.1)
                    near(s.y[:,0],[a,1-a,.4*(1-a)],rtol=tol)


ROBERTSON = np.array([
 [.9664597373330035,.3074626578578675e-4,.3350951640121071e-1],
 [.8413699238414729,.1623390937990473e-4,.1586138422491472],
 [.6172348823960878,.6153591274639123e-5,.3827589640126376],
 [.3368745306607069,.2013702318261393e-5,.6631234556369748],
 [.1073004285378040,.4800166972571660e-6,.8926990914454987],
 [.1786592114209946e-1,.7274751468436319e-7,.9821340061103859],
 [.2031483924973415e-2,.8142277783356159e-8,.9979685079327488],
 [.2076093439016395e-3,.8306077485067610e-9,.9997923898254906],
 [.2082417512179460e-4,.8329841429908955e-10,.9999791757415798],
 [.2083229471647004e-5,.8332935037760723e-11,.9999979167621954],
 [.2083328471883087e-6,.8333315602809495e-12,.9999997916663195],
 [.2083340149701284e-7,.8333360770334744e-13,.9999999791665152]])
E5 = np.array([
 [1.7599259497677897058e-3,1.3846281519376516449e-11,7.6370038530073911180e-13,1.3082581134075777338e-11],
 [1.6180769999072942552e-3,1.3822370304983735443e-10,8.2515735006838336088e-12,1.2997212954915352082e-10],
 [7.4813208224292220114e-6,2.3734781561205975019e-12,2.2123586689581663654e-12,1.6111948716243113653e-13],
 [4.7150333630401632232e-10,1.8188895860807021729e-14,1.8188812376786725407e-14,8.3484020296321693074e-20],
 [3.1317148329356996037e-14,1.4840957952870064294e-16,1.4840957948345691466e-16,4.5243728279782625194e-26],
 [3.8139035189787091771e-49,1.0192582567660293322e-20,1.0192582567660293322e-20,3.7844935507486221171e-65],
 [0,8.8612334976263783420e-23,8.8612334976263783421e-23,0]])
HIRES = np.array([
 [.000737131257332567,.000144248572631618,.000058887297409676,.001175651343283149,.002386356198831330,.006238968252742796,.002849998395185769,.002850001604814231],
 [.000670305503581864,.000130996846986347,.000046862231597733,.001044668020551705,.000594883830951485,.001399628833942774,.001014492757718480,.004685507242281520]])
OREGONATOR = np.array([
 [.1000661467180497e1,.1512778937348249e4,.1035854312767229e5],
 [.1000874625199626e1,.1144336972384497e4,.8372149966624639e2],
 [.1001890368438751e1,.5299926232295553e3,.1662279579042420e1],
 [.1004118022612645e1,.2438326079910346e3,.1008822224048647e1],
 [.1008995416634061e1,.1121664388662539e3,.1007783229065319e1],
 [.1019763472537298e1,.5159761322947535e2,.1016985778956374e1],
 [.1043985088527474e1,.2373442027531524e2,.1037691843544522e1],
 [.1100849071667922e1,.1091533805469020e2,.1085831969810860e1],
 [.1249102130020572e1,.5013945178605446e1,.1208326626237875e1],
 [.1779724751937019e1,.2281852385542403e1,.1613754023671725e1],
 [.1000889326903503e1,.1125438585746596e4,.1641049483777168e5],
 [.1000814870318523e1,.1228178521549889e4,.1320554942846513e3]])


def robertson(m):
    add_user(m,['A'],{'B':1},.04)
    add_user(m,['B','B'],{'B':1,'C':1},3e7)
    add_user(m,['B','C'],{'A':1,'C':1},1e4)


@group('analytical_rosenbrock_integration:AnalyticalExamples.','Robertson E5 Oregonator HIRES')
@group('analytical_backward_euler:AnalyticalExamples.','Oregonator HIRES')
def stiff_examples(id_):
    name = id_.split('.')[-1]; methods = ('be',) if id_.startswith('analytical_backward') else METHODS
    names = {'Robertson':['A','B','C'],'E5':['A'+str(i) for i in range(1,7)],
             'HIRES':['Y'+str(i) for i in range(8)],'Oregonator':['X','Y','Z','P','Q']}[name]
    with Model(names) as m:
        if name == 'Robertson':
            robertson(m); initial = [1,0,0]; rtol,atol = 1e-10,1e-12
            times,expected = 10.**np.arange(12),ROBERTSON
        elif name == 'E5':
            for rs,ps,k in [(['A1'],{'A2':1,'A3':1},7.89e-10),(['A2','A3'],{'A5':1},1.13e9),
                            (['A1','A3'],{'A4':1},1.1e7),(['A4'],{'A3':1,'A6':1},1.13e3)]:
                add_user(m,rs,ps,k)
            initial = [.00176,0,0,0,0,0]; rtol,atol = 1e-13,[1e-7,1e-17,1e-17,1e-17,1e-7,1e-7]
            times,expected = 10.**np.arange(1,14,2),E5
        elif name == 'HIRES':
            for rs,ps,k in [(['Y0'],{'Y1':1},1.71),(['Y1'],{'Y0':.43/8.75,'Y3':8.32/8.75},8.75),
                (['Y2'],{'Y0':8.32/10.03,'Y3':1.71/10.03},10.03),([],{'Y0':1},.0007),
                (['Y3'],{'Y2':.43/1.12,'Y5':.69/1.12},1.12),(['Y4'],{'Y2':.035/1.745,'Y5':1.71/1.745},1.745),
                (['Y5'],{'Y4':1},.43),(['Y6'],{'Y4':.43/1.81,'Y5':.69/1.81,'Y7':1},1.81),
                (['Y5','Y7'],{'Y6':1},280.)]: add_user(m,rs,ps,k)
            initial = [1,0,0,0,0,0,0,.0057]; rtol,atol = 1e-6,1e-8
            times,expected = np.array([321.8122,421.8122]),HIRES
        else:
            sv,qv,wv,k2 = 77.27,8.375e-6,.161,1.6e9
            for rs,ps,k in [(['Y'],{'X':1},1/(sv*wv)),(['X','Y'],{'P':1},k2),
                (['X'],[('Z',1),('X',2)],sv/wv),(['X','X'],{'Q':1},sv*sv*qv*k2/2),(['Z'],{'Y':1},1.)]:
                add_user(m,rs,ps,k)
            alpha,eta,rho = 1/(sv*wv*k2),sv/(wv*k2),1/(wv*wv*k2)
            initial = [alpha,2*eta,3*rho,0,0]; rtol,atol = 1e-9,[alpha*1e-11,eta*1e-11,rho*1e-11,eta,eta]
            times,expected = np.arange(1,13)*30*wv,OREGONATOR*np.array([alpha,eta,rho])
        m.finalize()
        for method in methods:
            with State(m) as s:
                s.conditions(272.5,101253.3,1e6); s.y[:,0] = initial; s.tolerance(rtol,atol); s.update()
                previous = 0.
                for i,t in enumerate(times):
                    advance(s,t-previous,method,18000 if name == 'Oregonator' and method == 'be' else 1)
                    previous = t; actual = s.y[:expected.shape[1],0].copy(); target = expected[i]
                    if name in ('Robertson','E5'):
                        absolute = .0003 if name == 'Robertson' else np.array([1e-6,1e-13,1e-13,1e-13])
                        relative = 1e-6 if name == 'Robertson' else 1e-10
                        assert np.all((abs(actual-target)<absolute)|(abs(actual-target)<relative*abs(target))), (name,method,t,actual,target)
                    elif name == 'HIRES':
                        tol = .1 if method == 'be' else (1e-6,1e-7,1e-7,1e-6,1e-6)[METHODS.index(method)]
                        near(actual,target,rtol=0,atol=tol)
                    else: near(actual,target,rtol=5e-3 if method == 'be' else 1e-6)
