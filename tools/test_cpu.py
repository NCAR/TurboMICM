"""Run ported cases and report unported cases as failures, never as passes."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
CASE_LIST = ROOT / '../reference/micm_cases.txt'
# Each listed executable must check every assertion in each mapped case.
CASES = {
    **{case: ('build/test_rates', 'rate calculation') for case in [
        'arrhenius_rate_constant:ArrheniusRateConstant.CalculateWithSystem',
        'arrhenius_rate_constant:ArrheniusRateConstant.CalculateWithPrescribedArguments',
        'tunneling_rate_constant:TunnelingRateConstant.CalculateWithMinimalArguments',
        'tunneling_rate_constant:TunnelingRateConstant.CalculateWithAllArguments',
        'user_defined_rate_constant:UserDefinedRateConstant.CalculateWithSystem',
        'user_defined_rate_constant:UserDefinedRateConstant.ConstructorWithRate',
        'user_defined_rate_constant:UserDefinedRateConstant.ConstructorWithRateAndName',
        'user_defined_rate_constant:UserDefinedRateConstant.CustomScalingFactor',
    ]},
    **{case: ('build/test_rate_extensions', 'rate calculation') for case in [
        'reversible_rate_constant:ReversibleRateConstant.DefaultConstructor',
        'reversible_rate_constant:ReversibleRateConstant.CalculateWithSystem',
        'reversible_rate_constant:ReversibleRateConstant.CalculateWithPrescribedArguments',
        'reversible_rate_constant:ReversibleRateConstant.TemperatureDependence',
        'reversible_rate_constant:ReversibleRateConstant.ZeroReverseRate',
        'reversible_rate_constant:ReversibleRateConstant.EquilibriumConstantOnly',
        'taylor_series_rate_constant:TaylorSeriesRateConstant.CalculateWithSystem',
        'taylor_series_rate_constant:TaylorSeriesRateConstant.CalculateWithPrescribedArguments',
        'ternary_chemical_activation_rate_constant:TernaryChemicalActivationJPL19RateConstant.CalculateWithAllArguments',
        'ternary_chemical_activation_rate_constant:TernaryChemicalActivationJPL19RateConstant.ZeroKintReducesToTroe',
    ]},
    **{case: ('build/test_falloff', 'rate calculation') for case in [
        'branched_rate_constant:BranchedRateConstant.CalculateAlkoxyBranchWithAllArguments',
        'branched_rate_constant:BranchedRateConstant.CalculateNitrateBranchWithAllArguments',
        'surface_rate_constant:SurfaceRateConstant.CalculateDefaultProbability',
        'surface_rate_constant:SurfaceRateConstant.CalculateSpecifiedProbability',
        'ternary_chemical_activation_rate_constant:TernaryChemicalActivationRateConstant.CalculateWithMinimalArguments',
        'ternary_chemical_activation_rate_constant:TernaryChemicalActivationRateConstant.CalculateWithAllArguments',
        'troe_rate_constant:TroeRateConstant.CalculateWithMinimalArguments',
        'troe_rate_constant:TroeRateConstant.CalculateWithAllArguments',
        'troe_rate_constant:TroeRateConstant.AnalyticalTroeExampleAB',
        'troe_rate_constant:TroeRateConstant.AnalyticalTroeExampleBC',
    ]},
    **{f'{prefix}:{suite}.{case}': ('build/test_doolittle', 'LU reconstruction')
       for prefix, suite in [
           ('lu_decomposition_doolittle', 'LuDecompositionDoolittle'),
           ('lu_decomposition_doolittle_in_place', 'LuDecompositionDoolittleInPlace'),
       ] for case in [
           'DenseMatrixStandardOrdering', 'RandomMatrixStandardOrdering',
           'DiagonalMatrixStandardOrdering', 'AgnosticToInitialValueStandardOrdering',
           'DenseMatrixVectorOrdering', 'RandomMatrixVectorOrdering',
           'DiagonalMatrixVectorOrdering', 'AgnosticToInitialValueVectorOrdering',
       ]},
    'analytical_rosenbrock_integration:AnalyticalExamples.Photolysis': ('build/test_analytical_methods', 'analytical trajectories'),
    'analytical_rosenbrock_integration:AnalyticalExamples.PhotolysisSuperStiffButAnalytical': ('build/test_analytical_methods', 'stiff analytical trajectories'),
    **{f'linear:LinearConstraint.{case}': ('build/test_constraints', 'linear constraint') for case in [
        'Construction', 'AlgebraicSpecies', 'WeightedTerms', 'ZeroConstant', 'FractionalCoefficients',
        'ResidualComputationThroughConstraintSet', 'JacobianComputationThroughConstraintSet',
        'WeightedSumResidualAndJacobian', 'ThreeSpeciesConservationResidual', 'ZeroConstantResidual',
        'FractionalCoefficientsResidualAndJacobian', 'JacobianIndependentOfConcentrations',
        'FiniteDifferenceJacobianSimpleConservation', 'FiniteDifferenceJacobianWeightedSum',
    ]},
    **{f'equilibrium:EquilibriumConstraint.{case}': ('build/test_constraints', 'equilibrium constraint') for case in [
        'Construction', 'AlgebraicSpecies', 'SingleReactantSingleProduct', 'MultipleReactantsAndProducts',
        'InvalidEquilibriumConstant', 'EmptyReactantsThrows', 'EmptyProductsThrows', 'InvalidStoichiometryThrows',
        'ResidualComputationThroughConstraintSet', 'JacobianComputationThroughConstraintSet',
        'ComplexStoichiometryResidual', 'FiniteDifferenceJacobianSimple', 'FiniteDifferenceJacobianComplexStoichiometry',
    ]},
    **{f'constraint_set:ConstraintSet.{case}': ('build/test_constraint_set', 'constraint set') for case in [
        'Construction', 'ReplaceStateRowsMapsToAlgebraicSpecies', 'NonZeroJacobianElements', 'MultipleConstraints',
        'AddForcingTerms', 'SubtractJacobianTerms', 'EmptyConstraintSet', 'UnknownSpeciesThrows',
        'ThreeDStateOneConstraint', 'FourDStateTwoConstraints', 'CoupledConstraintsSharedSpecies',
        'VectorizedMatricesRespectGridCellIndexing', 'VectorMatrix1', 'VectorMatrix2', 'VectorMatrix3',
    ]},
    'rosenbrock:RosenbrockSolver.StandardNormalizedError': ('build/test_normalized_error', 'normalized error'),
    'rosenbrock:RosenbrockSolver.VectorNormalizedError': ('build/test_normalized_error', 'normalized error'),
    'backward_euler:BackwardEuler.CanCallSolve': ('build/test_backward_euler', 'chain solve'),
    'backward_euler:BackwardEuler.IsConverged': ('build/test_backward_euler', 'convergence thresholds'),
}


def main():
    results = {}
    failures = 0
    for line in CASE_LIST.read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        case, _kind = line.split()
        if case not in CASES:
            status, reason = 'FAIL', 'not yet ported'
        else:
            executable, reason = CASES[case]
            if executable not in results:
                results[executable] = subprocess.run([str(ROOT / executable)],
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True).returncode
            status = 'PASS' if results[executable] == 0 else 'FAIL'
            if status == 'FAIL': reason = 'ported check failed'
        failures += status == 'FAIL'
        print(case, status, reason)
    return bool(failures)

if __name__ == '__main__':
    sys.exit(main())
