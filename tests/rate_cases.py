"""Ports of test/unit/process/rate_constant: same data and double tolerances.

Original assertions: Copyright (C) 2023-2026 University Corporation for Atmospheric Research.
SPDX-License-Identifier: Apache-2.0
"""
import ctypes as C
import math

from api import *
from harness import case, equal, near


@case('arrhenius_rate_constant:ArrheniusRateConstant.CalculateWithSystem',
      'arrhenius_rate_constant:ArrheniusRateConstant.CalculateWithPrescribedArguments')
def arrhenius(_):
    r = rate()
    near(value(r),1.0,1e-13)
    for a,c in [(1,0),(2.2e-10,0),(3e-11,-200),(2.6e-12,-350)]:
        r.A,r.C = a,c
        expected = a*math.exp(c/301.24)
        near(value(r),expected,1e-13*expected)


@case('branched_rate_constant:BranchedRateConstant.CalculateAlkoxyBranchWithAllArguments',
      'branched_rate_constant:BranchedRateConstant.CalculateNitrateBranchWithAllArguments')
def branched(identifier):
    branch = int('Nitrate' in identifier)
    r = rate(BRANCHED,branch=branch,X=1.2,Y=204.3,a0=1e-3,n=2)
    a = 2e-22*math.exp(2)*2.45e19
    b = .43*(293.0/298.0)**-8.0
    z = a/(1+a/b)*.41**(1/(1+math.log10(a/b)**2))*(1-1e-3)/1e-3
    a = 2e-22*math.exp(2)*(42.2*6.02214076e23*1e-6)
    b = .43*(301.24/298.0)**-8.0
    capital_a = a/(1+a/b)*.41**(1/(1+math.log10(a/b)**2))
    expected = 1.2*math.exp(-204.3/301.24)*((capital_a if branch else z)/(z+capital_a))
    near(value(r),expected,1e-13*expected)


@case(*['reversible_rate_constant:ReversibleRateConstant.'+name for name in (
    'DefaultConstructor','CalculateWithSystem','CalculateWithPrescribedArguments','TemperatureDependence','ZeroReverseRate','EquilibriumConstantOnly')])
def reversible(identifier):
    name = identifier.split('.')[-1]
    if name == 'DefaultConstructor':
        near(value(rate(REVERSIBLE),temperature=298.15),0.0,1e-13)
    elif name == 'TemperatureDependence':
        for t in (273.15,298.15,323.15):
            expected = 1.5*math.exp(2000/t)*.5
            near(value(rate(REVERSIBLE,A=1.5,C=2000,B=.5),temperature=t),expected,1e-13*expected)
    elif name == 'ZeroReverseRate':
        near(value(rate(REVERSIBLE,A=5,C=1000,B=0),temperature=298.15),0.0,1e-13)
    elif name == 'EquilibriumConstantOnly':
        expected = 3.5*math.exp(1800/298.15)
        near(value(rate(REVERSIBLE,A=3.5,C=1800,B=1),temperature=298.15),expected,1e-13*expected)
    else:
        t = 301.24 if name == 'CalculateWithSystem' else 298.15
        expected = 1.14e-2*math.exp(2300/t)*.32
        near(value(rate(REVERSIBLE,A=1.14e-2,C=2300,B=.32),temperature=t),expected,1e-13*expected)


@case('tunneling_rate_constant:TunnelingRateConstant.CalculateWithMinimalArguments',
      'tunneling_rate_constant:TunnelingRateConstant.CalculateWithAllArguments')
def tunneling(identifier):
    if 'Minimal' in identifier:
        near(value(rate(TUNNELING)),1,1e-13)
    else:
        expected = 1.2*math.exp(-2.3/301.24)*math.exp(302.3/301.24**3)
        near(value(rate(TUNNELING,A=1.2,B=2.3,C=302.3)),expected,1e-13*expected)


@case('taylor_series_rate_constant:TaylorSeriesRateConstant.CalculateWithSystem',
      'taylor_series_rate_constant:TaylorSeriesRateConstant.CalculateWithPrescribedArguments')
def taylor(_):
    r = rate(TAYLOR)
    near(value(r),1.0,1e-13)
    for a,c,coeffs in [(1,0,[1]),(2.2e-10,0,[1]),(3e-11,-200,[12.5,1.3e-2,5.2e-4]),(2.6e-12,-350,[1,.43,7.3e-3])]:
        r = rate(TAYLOR,A=a,C=c,coefficients=coeffs)
        expected = a*math.exp(c/301.24)*sum(v*301.24**j for j,v in enumerate(coeffs))
        near(value(r),expected,1e-13*expected)


@case(*['troe_rate_constant:TroeRateConstant.'+name for name in (
    'CalculateWithMinimalArguments','CalculateWithAllArguments','AnalyticalTroeExampleAB','AnalyticalTroeExampleBC')],
      *['ternary_chemical_activation_rate_constant:TernaryChemicalActivationRateConstant.'+name for name in (
          'CalculateWithMinimalArguments','CalculateWithAllArguments')])
def falloff(identifier):
    kind = TROE if identifier.startswith('troe') else TERNARY
    tolerance = 0.0
    if 'AllArguments' in identifier:
        r = rate(kind,k0_A=1.2,k0_B=2.3,k0_C=302.3,kinf_A=2.6,kinf_B=-3.1,kinf_C=402.1,Fc=.9,N=1.2)
    elif 'ExampleAB' in identifier:
        r = rate(kind,k0_A=4e-18)
    elif 'ExampleBC' in identifier:
        r = rate(kind,k0_A=1.2e-12,k0_B=167,k0_C=3,kinf_A=136,kinf_B=5,kinf_C=24,Fc=.9,N=.8)
        tolerance = 1e-14
    else:
        r = rate(kind)
    k0 = r.k0_A*math.exp(r.k0_C/301.24)*(301.24/300)**r.k0_B
    ki = r.kinf_A*math.exp(r.kinf_C/301.24)*(301.24/300)**r.kinf_B
    numerator = 42.2 if kind == TROE else 1
    expected = numerator*k0/(1+42.2*k0/ki)*r.Fc**(1/(1+1/r.N*math.log10(42.2*k0/ki)**2))
    near(value(r),expected,tolerance*abs(expected))


@case('ternary_chemical_activation_rate_constant:TernaryChemicalActivationJPL19RateConstant.CalculateWithAllArguments',
      'ternary_chemical_activation_rate_constant:TernaryChemicalActivationJPL19RateConstant.ZeroKintReducesToTroe')
def jpl19(identifier):
    if 'ZeroKint' in identifier:
        t,density = 250,30
        r = rate(JPL19,k0_A=3e-3,k0_B=-2,kinf_A=1.5,kinf_B=-1,kint_A=0)
    else:
        t,density = 301.24,42.2
        r = rate(JPL19,k0_A=1.2,k0_B=2.3,k0_C=302.3,k0_D=298,kinf_A=2.6,kinf_B=-3.1,kinf_C=402.1,kinf_D=310,
                 kint_A=.7,kint_B=.5,kint_C=120,kint_D=290,Fc=.9,N=1.2)
    k0 = r.k0_A*math.exp(r.k0_C/t)*(t/r.k0_D)**r.k0_B
    ki = r.kinf_A*math.exp(r.kinf_C/t)*(t/r.kinf_D)**r.kinf_B
    kint = r.kint_A*math.exp(r.kint_C/t)*(t/r.kint_D)**r.kint_B
    ratio = k0*density/ki
    kf = k0*density/(1+ratio)*r.Fc**(r.N/(r.N+math.log10(ratio)**2))
    expected = kf+kint*(1-kf/ki)
    near(value(r,temperature=t,air_density=density),expected,1e-30 if 'ZeroKint' in identifier else 0)


@case(*['user_defined_rate_constant:UserDefinedRateConstant.'+name for name in (
    'CalculateWithSystem','ConstructorWithRate','ConstructorWithRateAndName','CustomScalingFactor')])
def user_defined(identifier):
    if 'CalculateWithSystem' in identifier:
        v,scale = .5,1
    elif 'CustomScaling' in identifier:
        v,scale = 1.2,2
    else:
        v,scale = 1.1,1
    equal(value(rate(USER,A=scale,custom_index=0),custom=[v]),v*scale)


@case('lambda_rate_constant:LambdaRateConstant.CalculateWithSystem')
def lambda_rate(_):
    callback_type = C.CFUNCTYPE(D,C.POINTER(Conditions),DP,P)
    callback = callback_type(lambda c,u,ctx: 2*c.contents.temperature)
    r = rate(LAMBDA,function=C.cast(callback,P))
    near(value(r,temperature=300),600,1e-13*600)


@case('surface_rate_constant:SurfaceRateConstant.CalculateDefaultProbability',
      'surface_rate_constant:SurfaceRateConstant.CalculateSpecifiedProbability')
def surface(identifier):
    probability = .74 if 'Specified' in identifier else 1.0
    r = rate(SURFACE,diffusion=2.3e2,molecular_weight=.025,probability=probability,custom_index=0)
    gas_constant = 1.380649e-23*6.02214076e23
    expected = 4*2.5e6*math.pi*(1e-7)**2/(1e-7/2.3e2+4/(math.sqrt(8*gas_constant*273.65/(math.pi*.025))*probability))
    near(value(r,temperature=273.65,custom=[1e-7,2.5e6]),expected,1e-13*expected)


@case('surface_rate_constant:SurfaceRateConstant.DiffusionCoefficientIsMissing',
      'surface_rate_constant:SurfaceRateConstant.MolecularWeightIsMissing',
      'reaction_rate_store:ReactionRateConstantStore.SurfaceMissingDiffusionCoefficientThrows',
      'reaction_rate_store:ReactionRateConstantStore.SurfaceMissingMolecularWeightThrows')
def surface_missing(identifier):
    m = Model(['foo'])
    try:
        r = rate(SURFACE,label='foo')
        if 'Diffusion' in identifier:
            r.molecular_weight = .025
        else:
            r.diffusion = 2.3e2
        equal(lib.chem_model_reaction(m.p,inds([0]),1,inds([]),arr([]),0,r),1)
    finally:
        m.close()
