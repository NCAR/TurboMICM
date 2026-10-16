"""Rate tests translated from test/unit/process/rate_constant (double precision)."""
from math import exp, pow, log10, sqrt, pi
from ffi import *
from harness import case, group

T, AIR = 301.24,42.2


def calculate(rate,temperature=T,pressure=0,air=AIR,custom=()):
    return rate_calculate(C.byref(rate),C.byref(Conditions(temperature,pressure,air)),doubles(custom))


@group('arrhenius_rate_constant:ArrheniusRateConstant.', 'CalculateWithSystem CalculateWithPrescribedArguments')
def arrhenius(_):
    r = rate_default(0)
    near(calculate(r),1.)
    for a,c in [(1.,0.),(2.2e-10,0.),(3e-11,-200.),(2.6e-12,-350.)]:
        r.p[0],r.p[2] = a,c
        near(calculate(r),a*exp(c/T))


@group('taylor_series_rate_constant:TaylorSeriesRateConstant.', 'CalculateWithSystem CalculateWithPrescribedArguments')
def taylor(_):
    r = rate_default(8)
    near(calculate(r),1.)
    r.p[0] = 1
    near(calculate(r),1.)
    r.p[0] = 2.2e-10
    near(calculate(r),2.2e-10)
    for a,c,coeff in [(3e-11,-200.,[12.5,1.3e-2,5.2e-4]),(2.6e-12,-350.,[1.,.43,.0073])]:
        r.p[0],r.p[2] = a,c
        r.coefficient_count = 3
        r.p[5:8] = coeff
        near(calculate(r),a*exp(c/T)*(coeff[0]+coeff[1]*T+coeff[2]*T**2))


@group('tunneling_rate_constant:TunnelingRateConstant.', 'CalculateWithMinimalArguments CalculateWithAllArguments')
def tunneling(id_):
    r = rate_default(4)
    if id_.endswith('MinimalArguments'):
        near(calculate(r),1.)
    else:
        r.p[:3] = [1.2,2.3,302.3]
        near(calculate(r),1.2*exp(-2.3/T)*exp(302.3/pow(T,3)))


@group('branched_rate_constant:BranchedRateConstant.', 'CalculateAlkoxyBranchWithAllArguments CalculateNitrateBranchWithAllArguments')
def branched(id_):
    branch = int('Nitrate' in id_)
    r = rate_default(5)
    r.p[:5] = [1.2,204.3,.001,2,branch]
    rate_prepare(C.byref(r))
    def frac(temperature,air):
        a = 2e-22*exp(2)*air
        b = .43*pow(temperature/298.,-8.)
        return a/(1+a/b)*pow(.41,1/(1+pow(log10(a/b),2)))
    z = frac(293.,2.45e19)*.999/.001
    v = frac(T,AIR*6.02214076e23*1e-6)
    near(calculate(r),1.2*exp(-204.3/T)*(v if branch else z)/(z+v))


@group('troe_rate_constant:TroeRateConstant.', 'CalculateWithMinimalArguments CalculateWithAllArguments AnalyticalTroeExampleAB AnalyticalTroeExampleBC')
@group('ternary_chemical_activation_rate_constant:TernaryChemicalActivationRateConstant.', 'CalculateWithMinimalArguments CalculateWithAllArguments')
def falloff(id_):
    troe = id_.startswith('troe_')
    r = rate_default(1 if troe else 2)
    tolerance = 0
    if id_.endswith('AllArguments'):
        r.p[:8] = [1.2,2.3,302.3,2.6,-3.1,402.1,.9,1.2]
    elif id_.endswith('ExampleAB'):
        r.p[0] = 4e-18
    elif id_.endswith('ExampleBC'):
        r.p[:8] = [1.2e-12,167.,3.,136.,5.,24.,.9,.8]
        tolerance = 1e-14
    p = r.p
    k0 = p[0]*exp(p[2]/T)*pow(T/300.,p[1])
    ki = p[3]*exp(p[5]/T)*pow(T/300.,p[4])
    numerator = AIR*k0 if troe and 'Calculate' in id_ else k0*(AIR if troe else 1)
    expected = numerator/(1+AIR*k0/ki)*pow(p[6],1/(1+1/p[7]*pow(log10(AIR*k0/ki),2)))
    near(calculate(r),expected,rtol=tolerance)


@group('ternary_chemical_activation_rate_constant:TernaryChemicalActivationJPL19RateConstant.', 'CalculateWithAllArguments ZeroKintReducesToTroe')
def jpl(id_):
    r = rate_default(3)
    t,air = T,AIR
    if id_.endswith('AllArguments'):
        r.p[:14] = [1.2,2.3,302.3,298.,2.6,-3.1,402.1,310.,.7,.5,120.,290.,.9,1.2]
    else:
        r.p[0],r.p[1],r.p[4],r.p[5],r.p[8] = .003,-2.,1.5,-1.,0.
        t,air = 250.,30.
    p = r.p
    rates = [p[i]*exp(p[i+2]/t)*pow(t/p[i+3],p[i+1]) for i in (0,4,8)]
    k0,ki,km = rates
    ratio = k0*air/ki
    kf = k0*air/(1+ratio)*pow(p[12],p[13]/(p[13]+pow(log10(ratio),2)))
    near(calculate(r,t,air=air),kf+km*(1-kf/ki),rtol=0,atol=1e-30 if km == 0 else 0)


@group('reversible_rate_constant:ReversibleRateConstant.', 'DefaultConstructor CalculateWithSystem CalculateWithPrescribedArguments TemperatureDependence ZeroReverseRate EquilibriumConstantOnly')
def reversible(id_):
    name = id_.split('.')[-1]
    r = rate_default(9)
    if name == 'DefaultConstructor':
        near(calculate(r,298.15),0.,atol=1e-13)
        return
    p,temperatures = {
        'CalculateWithSystem': ([.0114,2300.,.32],[T]),
        'CalculateWithPrescribedArguments': ([.0114,2300.,.32],[298.15]),
        'TemperatureDependence': ([1.5,2000.,.5],[273.15,298.15,323.15]),
        'ZeroReverseRate': ([5.,1000.,0.],[298.15]),
        'EquilibriumConstantOnly': ([3.5,1800.,1.],[298.15]),
    }[name]
    r.p[:3] = p
    for t in temperatures:
        near(calculate(r,t),p[0]*exp(p[1]/t)*p[2],atol=1e-13 if p[2] == 0 else 0)


@group('user_defined_rate_constant:UserDefinedRateConstant.', 'CalculateWithSystem ConstructorWithRate ConstructorWithRateAndName CustomScalingFactor')
def user(id_):
    r = rate_default(7)
    r.custom_index = 0
    value = .5 if id_.endswith('System') else 1.1
    if id_.endswith('CustomScalingFactor'):
        r.p[0],value = 2.,1.2
    near(calculate(r,custom=[value]),value*r.p[0],rtol=0)


@case('lambda_rate_constant:LambdaRateConstant.CalculateWithSystem')
def lambda_rate(_):
    callback = C.CFUNCTYPE(D,C.POINTER(Conditions),DP,P)(lambda cond,custom,context: 2*cond.contents.temperature)
    r = rate_default(10)
    r.function = C.cast(callback,P)
    near(calculate(r,300.),600.)


@group('surface_rate_constant:SurfaceRateConstant.', 'CalculateDefaultProbability CalculateSpecifiedProbability DiffusionCoefficientIsMissing MolecularWeightIsMissing')
def surface(id_):
    if id_.endswith('IsMissing'):
        with Model(['foo']) as model:
            key = 'molecular weight [kg mol-1]' if 'Diffusion' in id_ else 'diffusion coefficient [m2 s-1]'
            require(model_set_species_property(model.handle,0,key.encode(),.025 if 'Diffusion' in id_ else 230.))
            model.parameter('radius'); model.parameter('number')
            r = rate_default(6); r.custom_index = 0
            assert model_add_reaction(model.handle,1,ints([0]),0,None,None,C.byref(r)) == -3
        return
    probability = .74 if 'Specified' in id_ else 1.
    r = rate_default(6)
    r.p[:3] = [230.,8*(1.380649e-23*6.02214076e23)/(pi*.025),probability]
    r.custom_index = 0
    expected = 4*2.5e6*pi*pow(1e-7,2)/(1e-7/230.+4/(sqrt(8*(1.380649e-23*6.02214076e23)*273.65/(pi*.025))*probability))
    near(calculate(r,273.65,custom=[1e-7,2.5e6]),expected)
