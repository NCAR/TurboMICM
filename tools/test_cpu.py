"""Run ported cases and report unported cases as failures, never as passes."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
CASE_LIST = ROOT / '../reference/micm_cases.txt'
# Each listed executable must check every assertion in each mapped case.
CASES = {
    **{f"{prefix}:{suite}.{case}": ("build/test_matrix", "matrix operations")
       for prefix, suite, cases in [
           ("matrix", "Matrix", ["SmallMatrix", "SmallConstMatrix", "InitializeMatrix", "InitializeConstMatrix", "LoopOverMatrix", "LoopOverConstMatrix"]),
           ("vector_matrix", "VectorMatrix", ["SmallVectorMatrix", "SmallConstVectorMatrix", "InitializeVectorMatrix", "InitializeConstVectorMatrix", "LoopOverVectorMatrix", "LoopOverConstVectorMatrix"]),
       ] for case in cases},
    **{f"{prefix}:{suite}.{case}": ("build/test_matrix", "matrix operations")
       for prefix, suite in [("matrix", "Matrix"), ("vector_matrix", "VectorMatrix")]
       for case in ["Strides", "ConversionToVector", "ConstConversionToVector", "ConversionFromVector", "AssignmentFromVector",
                    "Axpy", "ForEach", "SetScaler", "Max", "Min", "Print", "TestFill", "TestCopy"]},
    **{f"external_model_constraints:ExternalModelConstraints.{case}": ("build/test_external_composition", "composed constraints")
       for case in ["KineticVsConservativeConstraintConvergence", "SimpleConstraintPreservesRatioNotMass",
                    "BuiltInVsExternalModelConstraintStepByStep", "MultiEquilibriumKineticVsComposedConstraints",
                    "TemperatureDependentConstraintParameter"]},
    **{f"external_model_constraints:ExternalModelFiniteDifferenceJacobian.{case}": ("build/test_external_constraints", "external model Jacobian")
       for case in ["ProcessForcingJacobian", "ConstraintResidualJacobian", "EquilibriumConstraintModelJacobian"]},
    **{f"external_model_constraints:ExternalModelConstraints.{case}": ("build/test_external_constraints", "external model constraints")
       for case in ["DAESolveEnforcesConservation", "AlgebraicRowHoldsOnlyConstraintTerms", "AddExternalModelConstraintsOnly", "MultiGridCell"]},
    **{f"regression_test_rosenbrock_p_force:RegressionRosenbrock.{case}": ("build/test_nonlinear", "Chapman rates and forcing")
       for case in ["RateConstants", "VectorRateConstants", "Forcing", "VectorForcing"]},
    **{f"regression_test_rosenbrock_dforce_dy:RegressionRosenbrock.{case}": ("build/test_nonlinear", "Chapman Jacobian")
       for case in ["Jacobian", "VectorJacobian"]},
    "chapman_integration:ChapmanIntegration.CanBuildChapmanSystem": ("build/test_nonlinear", "Chapman conservation"),
    "integrated_reaction_rates:ChapmanIntegration.CanBuildChapmanSystem": ("build/test_nonlinear", "integrated reaction rates"),
    **{f"terminator:RosenbrockSolver.{name}": ("build/test_nonlinear", "terminator trajectories") for name in ["Terminator", "VectorTerminator"]},
    "analytical_rosenbrock_integration:AnalyticalExamples.SurfaceRxn": ("build/test_nonlinear", "surface reaction trajectories"),
    "analytical_backward_euler:AnalyticalExamples.SurfaceRxn": ("build/test_nonlinear", "surface reaction trajectories"),
    "analytical_rosenbrock_integration:AnalyticalExamples.E5": ("build/test_nonlinear", "E5 trajectories"),
    "analytical_rosenbrock_integration:AnalyticalExamples.Oregonator": ("build/test_nonlinear", "Oregonator trajectories"),
    "analytical_backward_euler:AnalyticalExamples.Oregonator": ("build/test_nonlinear", "Oregonator trajectories"),
    "analytical_rosenbrock_integration:AnalyticalExamples.Robertson": ("build/test_nonlinear", "Robertson trajectories"),
    "analytical_rosenbrock_integration:AnalyticalExamples.HIRES": ("build/test_nonlinear", "HIRES trajectories"),
    "analytical_backward_euler:AnalyticalExamples.HIRES": ("build/test_nonlinear", "HIRES trajectories"),
    **{f"process_set:ProcessSet.{case}": ("build/test_process_set", "process forcing and Jacobian") for case in ["Matrix", "VectorMatrix"]},
    "process_set:RandomProcessSet.Matrix": ("build/test_process_random", "random process systems"),
    **{f"{prefix}:AnalyticalExamples.{family}{suffix}": ("build/test_chain_rates", "analytical reaction trajectories")
       for prefix in ["analytical_rosenbrock_integration", "analytical_backward_euler"]
       for family in ["Troe", "TernaryChemicalActivation", "Tunneling", "Arrhenius", "Branched", "Photolysis"]
       for suffix in ["", "SuperStiffButAnalytical"]},
    **{f"state:State.{case}": ("build/test_state", "state ownership and named updates") for case in [
        "DefaultConstructor", "Constructor", "CopyConstructor", "CopyAssignmentOperator", "MoveConstructor", "MoveAssignmentOperator",
        "SettingSingleConcentrationWithInvalidArgumentsThowsException", "SetSingleConcentration", "SetConcentrationByElementSingleValue",
        "SetConcentrationByElementVector", "SettingConcentrationsWithInvalidArguementsThrowsException", "SetConcentrations",
        "SetStateWithSquareBracketOperator", "SetMultiCellStateWithSquareBrackets", "SettingCustomRateParameterWithInvalidVectorSizeThrowsException",
        "SettingCustomRateParameterWithInvalidLabelThrowsException", "SetCustomRateParameter", "SetCustomRateParameters",
        "UnsafelySetCustomRateParameterOneCell", "UnsafelySetCustomRateParameterMultiCell",
        "UnsafelySetCustomRateParameterCatchesTooFewGridCells", "UnsafelySetCustomRateParameterCatchesTooParameters",
    ]},
    "process_set:ProcessSetAlgebraicVariables.CudaMatrix": ("build/test_process_set", "algebraic masking"),
    "process_set:ProcessSetFiniteDifferenceJacobian.Matrix": ("build/test_process_set", "process Jacobian"),
    **{f"jacobian_verification:JacobianVerification.{case}": ("build/test_jacobian_check", "Jacobian verification") for case in [
        "FiniteDifferenceMatchesAnalytical", "MultiBlockFiniteDifference", "ComparePassesForCorrectJacobian",
        "CompareFailsForWrongJacobian", "SparsityCompletenessPassesWhenComplete",
        "SparsityCompletenessFailsWhenMissingEntry", "NearZeroVariableHandled",
    ]},
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
           ('lu_decomposition_mozart', 'LuDecompositionMozart'),
           ('lu_decomposition_mozart_in_place', 'LuDecompositionMozartInPlace'),
       ] for case in [
           'DenseMatrixStandardOrdering', 'RandomMatrixStandardOrdering',
           'DiagonalMatrixStandardOrdering', 'AgnosticToInitialValueStandardOrdering',
           'DenseMatrixVectorOrdering', 'RandomMatrixVectorOrdering',
           'DiagonalMatrixVectorOrdering', 'AgnosticToInitialValueVectorOrdering',
       ]},
    **{f'{prefix}:{suite}.{case}': ('build/test_linear_solve', 'linear solve')
       for prefix, suite in [('linear_solver', 'LinearSolver'), ('linear_solver_in_place', 'LinearSolverInPlace')]
       for case in ['DenseMatrixStandardOrdering', 'RandomMatrixStandardOrdering', 'DiagonalMatrixStandardOrdering',
                    'StandardOrderingAgnosticToInitialValue', 'DenseMatrixVectorOrdering', 'RandomMatrixVectorOrdering',
                    'VectorOrderingAgnosticToInitialValue', 'DiagonalMatrixVectorOrdering']},
    **{f'{prefix}:{suite}.{case}': ('tests/test_markowitz.py', 'Markowitz fill reduction')
       for prefix, suite in [('linear_solver', 'LinearSolver'), ('linear_solver_in_place', 'LinearSolverInPlace')]
       for case in ['DiagonalMarkowitzReorder', 'VectorDiagonalMarkowitzReordering']},
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
    **{f'constraint_initialization:ConstraintInitialization.{case}': ('build/test_constraint_init', 'constraint projection') for case in [
        'ConsistentICsUnchanged', 'MildlyInconsistentICsCorrected', 'SeverelyInconsistentICsConverge',
        'ScaledConstraintRowsSolveIdenticallyThroughSolve', 'PureODESystemUnaffected', 'MultiCellSystems',
        'SubsequentSolveCallsReinitialize', 'SolverStateToStringNewValue',
        'FailedInitializationRestoresCallerState', 'SuccessfulInitializationKeepsTheCorrection',
        'WeightedCorrectionIsInvariantToConstraintRowScaling', 'WeightedCorrectionUsesStateTolerances',
    ]},
    **{f'equilibrium_constraint:EquilibriumIntegration.{case}': ('build/test_dae', 'DAE integration') for case in [
        'DAESolveWithConstraint', 'DAESolveWithTwoCoupledConstraints', 'DAESolveWithNonUnitStoichiometry',
    ]},
    **{case: ('build/test_parameterized', 'parameterized reaction rates') for case in [
        'process:Process.Matrix', 'process:Process.VectorMatrix',
        'reaction_rate_store:ReactionRateConstantStore.ParameterizedMultipliers',
    ]},
    **{case: ('build/test_callback', 'callback rate') for case in [
        'lambda_rate_constant:LambdaRateConstant.CalculateWithSystem',
        'reaction_rate_store:ReactionRateConstantStore.LambdaEntriesRcIndex',
        'solver_lambda_rate_constant:Solver.GetLambdaRateConstantByNameCanOverrideLambda',
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
                command = [str(ROOT / executable)]
                if executable.endswith('.py'):
                    command.insert(0, str(ROOT / '.conda/bin/python'))
                results[executable] = subprocess.run(command,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True).returncode
            status = 'PASS' if results[executable] == 0 else 'FAIL'
            if status == 'FAIL': reason = 'ported check failed'
        failures += status == 'FAIL'
        print(case, status, reason)
    return bool(failures)

if __name__ == '__main__':
    sys.exit(main())
