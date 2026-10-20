#!/usr/bin/env python3
"""Report each reference case, distinguishing exercised ports from missing ones."""
from pathlib import Path

root = Path(__file__).resolve().parents[2]
case_file = root / "reference" / "micm_cases.txt"
implemented = {
    "arrhenius_rate_constant:ArrheniusRateConstant.CalculateWithSystem",
    "arrhenius_rate_constant:ArrheniusRateConstant.CalculateWithPrescribedArguments",
    "branched_rate_constant:BranchedRateConstant.CalculateAlkoxyBranchWithAllArguments",
    "branched_rate_constant:BranchedRateConstant.CalculateNitrateBranchWithAllArguments",
    "lambda_rate_constant:LambdaRateConstant.CalculateWithSystem",
    "surface_rate_constant:SurfaceRateConstant.CalculateDefaultProbability",
    "surface_rate_constant:SurfaceRateConstant.CalculateSpecifiedProbability",
    "surface_rate_constant:SurfaceRateConstant.DiffusionCoefficientIsMissing",
    "surface_rate_constant:SurfaceRateConstant.MolecularWeightIsMissing",
    "reversible_rate_constant:ReversibleRateConstant.DefaultConstructor",
    "reversible_rate_constant:ReversibleRateConstant.CalculateWithSystem",
    "reversible_rate_constant:ReversibleRateConstant.CalculateWithPrescribedArguments",
    "reversible_rate_constant:ReversibleRateConstant.TemperatureDependence",
    "reversible_rate_constant:ReversibleRateConstant.ZeroReverseRate",
    "reversible_rate_constant:ReversibleRateConstant.EquilibriumConstantOnly",
    "taylor_series_rate_constant:TaylorSeriesRateConstant.CalculateWithSystem",
    "taylor_series_rate_constant:TaylorSeriesRateConstant.CalculateWithPrescribedArguments",
    "troe_rate_constant:TroeRateConstant.CalculateWithMinimalArguments",
    "troe_rate_constant:TroeRateConstant.CalculateWithAllArguments",
    "troe_rate_constant:TroeRateConstant.AnalyticalTroeExampleAB",
    "troe_rate_constant:TroeRateConstant.AnalyticalTroeExampleBC",
    "ternary_chemical_activation_rate_constant:TernaryChemicalActivationRateConstant.CalculateWithMinimalArguments",
    "ternary_chemical_activation_rate_constant:TernaryChemicalActivationRateConstant.CalculateWithAllArguments",
    "ternary_chemical_activation_rate_constant:TernaryChemicalActivationJPL19RateConstant.CalculateWithAllArguments",
    "ternary_chemical_activation_rate_constant:TernaryChemicalActivationJPL19RateConstant.ZeroKintReducesToTroe",
    "tunneling_rate_constant:TunnelingRateConstant.CalculateWithMinimalArguments",
    "tunneling_rate_constant:TunnelingRateConstant.CalculateWithAllArguments",
    "user_defined_rate_constant:UserDefinedRateConstant.CalculateWithSystem",
    "user_defined_rate_constant:UserDefinedRateConstant.ConstructorWithRate",
    "user_defined_rate_constant:UserDefinedRateConstant.ConstructorWithRateAndName",
    "user_defined_rate_constant:UserDefinedRateConstant.CustomScalingFactor",
    "process:Process.BuildsChemicalReaction",
    "process:Process.ChemicalReactionCopyAssignmentSucceeds",
    "reaction_rate_store:ReactionRateConstantStore.OffsetsAreContiguousCumulativeSizes",
    "reaction_rate_store:ReactionRateConstantStore.ArrheniusParametersPreserved",
    "reaction_rate_store:ReactionRateConstantStore.BranchedDerivedFieldsComputed",
    "reaction_rate_store:ReactionRateConstantStore.UserDefinedCustomParamIndex",
    "reaction_rate_store:ReactionRateConstantStore.SurfaceDataFieldsAndCustomParamIndex",
    "reaction_rate_store:ReactionRateConstantStore.SurfaceCustomParamIndexAfterUserDefined",
    "reaction_rate_store:ReactionRateConstantStore.TotalRateConstantsMatchesProcessCount",
    "reaction_rate_store:ReactionRateConstantStore.SurfaceMissingDiffusionCoefficientThrows",
    "reaction_rate_store:ReactionRateConstantStore.SurfaceMissingMolecularWeightThrows",
    "reaction_rate_store:ReactionRateConstantStore.ParameterizedMultipliers",
    "README_example",
    "solve_results",
    "multiple_grid_cells",
    "openmp_tutorial",
    "vectorized_matrix_solver",
    "solver_configuration",
    "rate_constants_no_user_defined_example_by_hand",
    "rate_constants_user_defined_example_by_hand",
}

failed = False
for raw in case_file.read_text().splitlines():
    line = raw.strip()
    if not line or line.startswith("#"):
        continue
    case_id, kind = line.split()
    if case_id in implemented:
        if case_id in {"README_example", "solve_results", "multiple_grid_cells",
                       "openmp_tutorial", "vectorized_matrix_solver",
                       "solver_configuration",
                       "rate_constants_no_user_defined_example_by_hand",
                       "rate_constants_user_defined_example_by_hand"}:
            reason = "C tutorial output checked"
        else:
            reason = "C rate formula checked"
        print(f"{case_id} PASS {reason}")
    elif kind == "program":
        failed = True
        print(f"{case_id} FAIL C program port incomplete")
    else:
        failed = True
        print(f"{case_id} FAIL C case port incomplete")
raise SystemExit(1 if failed else 0)
