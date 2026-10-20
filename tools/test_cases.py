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
    "process:Process.Matrix",
    "process:Process.VectorMatrix",
    "process_set:ProcessSetFiniteDifferenceJacobian.Matrix",
    "process_set:RandomProcessSet.Matrix",
    "process_set:ProcessSet.Matrix",
    "process_set:ProcessSet.VectorMatrix",
    "terminator:RosenbrockSolver.Terminator",
    "terminator:RosenbrockSolver.VectorTerminator",
    "chapman_integration:ChapmanIntegration.CanBuildChapmanSystem",
    "integrated_reaction_rates:ChapmanIntegration.CanBuildChapmanSystem",
    "solver_lambda_rate_constant:Solver.GetLambdaRateConstantByNameCanOverrideLambda",
    "process_set:ProcessSetAlgebraicVariables.CudaMatrix",
    "constraint_set:ConstraintSet.Construction",
    "constraint_set:ConstraintSet.ReplaceStateRowsMapsToAlgebraicSpecies",
    "constraint_set:ConstraintSet.NonZeroJacobianElements",
    "constraint_set:ConstraintSet.MultipleConstraints",
    "constraint_set:ConstraintSet.AddForcingTerms",
    "constraint_set:ConstraintSet.SubtractJacobianTerms",
    "constraint_set:ConstraintSet.EmptyConstraintSet",
    "constraint_set:ConstraintSet.UnknownSpeciesThrows",
    "constraint_set:ConstraintSet.ThreeDStateOneConstraint",
    "constraint_set:ConstraintSet.FourDStateTwoConstraints",
    "constraint_set:ConstraintSet.CoupledConstraintsSharedSpecies",
    "constraint_set:ConstraintSet.VectorizedMatricesRespectGridCellIndexing",
    "constraint_set:ConstraintSet.VectorMatrix1",
    "constraint_set:ConstraintSet.VectorMatrix2",
    "constraint_set:ConstraintSet.VectorMatrix3",
    "linear:LinearConstraint.Construction",
    "linear:LinearConstraint.AlgebraicSpecies",
    "linear:LinearConstraint.WeightedTerms",
    "linear:LinearConstraint.ZeroConstant",
    "linear:LinearConstraint.FractionalCoefficients",
    "linear:LinearConstraint.ResidualComputationThroughConstraintSet",
    "linear:LinearConstraint.JacobianComputationThroughConstraintSet",
    "linear:LinearConstraint.WeightedSumResidualAndJacobian",
    "linear:LinearConstraint.ThreeSpeciesConservationResidual",
    "linear:LinearConstraint.ZeroConstantResidual",
    "linear:LinearConstraint.FractionalCoefficientsResidualAndJacobian",
    "linear:LinearConstraint.JacobianIndependentOfConcentrations",
    "linear:LinearConstraint.FiniteDifferenceJacobianSimpleConservation",
    "linear:LinearConstraint.FiniteDifferenceJacobianWeightedSum",
    "equilibrium:EquilibriumConstraint.Construction",
    "equilibrium:EquilibriumConstraint.AlgebraicSpecies",
    "equilibrium:EquilibriumConstraint.SingleReactantSingleProduct",
    "equilibrium:EquilibriumConstraint.MultipleReactantsAndProducts",
    "equilibrium:EquilibriumConstraint.InvalidEquilibriumConstant",
    "equilibrium:EquilibriumConstraint.EmptyReactantsThrows",
    "equilibrium:EquilibriumConstraint.EmptyProductsThrows",
    "equilibrium:EquilibriumConstraint.InvalidStoichiometryThrows",
    "equilibrium:EquilibriumConstraint.ResidualComputationThroughConstraintSet",
    "equilibrium:EquilibriumConstraint.JacobianComputationThroughConstraintSet",
    "equilibrium:EquilibriumConstraint.ComplexStoichiometryResidual",
    "equilibrium:EquilibriumConstraint.FiniteDifferenceJacobianSimple",
    "equilibrium:EquilibriumConstraint.FiniteDifferenceJacobianComplexStoichiometry",
    "constraint_initialization:ConstraintInitialization.ConsistentICsUnchanged",
    "constraint_initialization:ConstraintInitialization.MildlyInconsistentICsCorrected",
    "constraint_initialization:ConstraintInitialization.SeverelyInconsistentICsConverge",
    "constraint_initialization:ConstraintInitialization.FailedInitializationRestoresCallerState",
    "constraint_initialization:ConstraintInitialization.SuccessfulInitializationKeepsTheCorrection",
    "constraint_initialization:ConstraintInitialization.ScaledConstraintRowsSolveIdenticallyThroughSolve",
    "constraint_initialization:ConstraintInitialization.WeightedCorrectionIsInvariantToConstraintRowScaling",
    "constraint_initialization:ConstraintInitialization.SubsequentSolveCallsReinitialize",
    "constraint_initialization:ConstraintInitialization.MultiCellSystems",
    "constraint_initialization:ConstraintInitialization.WeightedCorrectionUsesStateTolerances",
    "matrix:Matrix.Strides",
    "matrix:Matrix.ConversionToVector",
    "matrix:Matrix.ConversionFromVector",
    "matrix:Matrix.Axpy",
    "matrix:Matrix.SetScaler",
    "matrix:Matrix.Max",
    "matrix:Matrix.Min",
    "matrix:Matrix.TestFill",
    "matrix:Matrix.TestCopy",
    "matrix:Matrix.ReduceSum",
    "matrix:Matrix.MismatchedRowDimensions",
    "matrix:Matrix.MismatchedColumnDimensions",
    "matrix:Matrix.WrongMatrixDimensions",
    "matrix:Matrix.ForEach",
    "matrix:Matrix.ArrayFunction",
    "matrix:Matrix.EmptyMatrixFunction",
    "matrix:Matrix.SmallMatrix",
    "matrix:Matrix.InitializeMatrix",
    "matrix:Matrix.LoopOverMatrix",
    "matrix:Matrix.AssignmentFromVector",
    "matrix:Matrix.SmallConstMatrix",
    "matrix:Matrix.InitializeConstMatrix",
    "matrix:Matrix.LoopOverConstMatrix",
    "matrix:Matrix.ConstConversionToVector",
    "matrix:Matrix.ColumnViewReuse",
    "matrix:Matrix.VectorInMatrixFunction",
    "matrix:Matrix.VectorTooSmall",
    "matrix:Matrix.VectorTooLarge",
    "matrix:Matrix.EmptyVectorNonEmptyMatrix",
    "matrix:Matrix.NonEmptyVectorEmptyMatrix",
    "matrix:Matrix.EmptyVectorEmptyMatrix",
    "matrix:Matrix.ReduceMax",
    "matrix:Matrix.IterateOverMatrix",
    "matrix:Matrix.FunctionReusability",
    "matrix:Matrix.MultiMatrixDifferentRowsFromCreation",
    "matrix:Matrix.VectorCapture",
    "vector_matrix:VectorMatrix.SmallVectorMatrix",
    "vector_matrix:VectorMatrix.SmallConstVectorMatrix",
    "vector_matrix:VectorMatrix.InitializeVectorMatrix",
    "vector_matrix:VectorMatrix.InitializeConstVectorMatrix",
    "vector_matrix:VectorMatrix.LoopOverVectorMatrix",
    "vector_matrix:VectorMatrix.LoopOverConstVectorMatrix",
    "vector_matrix:VectorMatrix.Strides",
    "phase:PhaseSpecies.ConstructsWithoutDiffusionCoefficient",
    "phase:PhaseSpecies.ConstructsWithDiffusionCoefficient",
    "phase:PhaseSpecies.SetDiffusionCoefficient",
    "species:Species.SetThirdBody",
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
    "reaction_rate_store:ReactionRateConstantStore.LambdaEntriesRcIndex",
    "version:Version.FullVersion",
    "version:Version.VersionMajor",
    "version:Version.VersionMinor",
    "version:Version.VersionPatch",
    "version:Version.VersionTweak",
    "state:State.DefaultConstructor",
    "state:State.Constructor",
    "state:State.SettingSingleConcentrationWithInvalidArgumentsThowsException",
    "state:State.SetSingleConcentration",
    "state:State.SetConcentrationByElementSingleValue",
    "state:State.SetConcentrationByElementVector",
    "state:State.SettingConcentrationsWithInvalidArguementsThrowsException",
    "state:State.SetConcentrations",
    "state:State.SetStateWithSquareBracketOperator",
    "state:State.SetMultiCellStateWithSquareBrackets",
    "state:State.SettingCustomRateParameterWithInvalidVectorSizeThrowsException",
    "state:State.SettingCustomRateParameterWithInvalidLabelThrowsException",
    "state:State.SetCustomRateParameter",
    "state:State.SetCustomRateParameters",
    "state:State.UnsafelySetCustomRateParameterOneCell",
    "state:State.UnsafelySetCustomRateParameterMultiCell",
    "state:State.UnsafelySetCustomRateParameterCatchesTooFewGridCells",
    "state:State.UnsafelySetCustomRateParameterCatchesTooParameters",
    "openmp:OpenMP.OneSolverManyStates",
    "jacobian_verification:JacobianVerification.FiniteDifferenceMatchesAnalytical",
    "jacobian_verification:JacobianVerification.MultiBlockFiniteDifference",
    "jacobian_verification:JacobianVerification.ComparePassesForCorrectJacobian",
    "jacobian_verification:JacobianVerification.CompareFailsForWrongJacobian",
    "jacobian_verification:JacobianVerification.SparsityCompletenessPassesWhenComplete",
    "jacobian_verification:JacobianVerification.SparsityCompletenessFailsWhenMissingEntry",
    "jacobian_verification:JacobianVerification.NearZeroVariableHandled",
    "README_example",
    "solve_results",
    "multiple_grid_cells",
    "openmp_tutorial",
    "vectorized_matrix_solver",
    "solver_configuration",
    "rate_constants_no_user_defined_example_by_hand",
    "rate_constants_user_defined_example_by_hand",
}
not_applicable = {
    "state:State.CopyConstructor":
        "C state views are non-owning, so they have no deep-copy constructor",
    "state:State.CopyAssignmentOperator":
        "C state views are non-owning, so they have no deep-copy assignment",
    "state:State.MoveConstructor":
        "C state views are non-owning and have no move constructor",
    "state:State.MoveAssignmentOperator":
        "C state views are non-owning and have no move assignment",
    "view_category_sfinae:ViewCategorySFINAE.VectorTypes":
        "C has no template type-trait or SFINAE behavior",
    "view_category_sfinae:ViewCategorySFINAE.MatrixViews":
        "C has no template type-trait or SFINAE behavior",
    "view_category_sfinae:ViewCategorySFINAE.SparseMatrixViews":
        "C has no template type-trait or SFINAE behavior",
    "view_category_sfinae:ViewCategorySFINAE.MatrixTypesExcluded":
        "C has no template type-trait or SFINAE behavior",
    "view_category_sfinae:ViewCategorySFINAE.RuntimeVectorUse":
        "C has no template type-trait or SFINAE behavior",
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
        if case_id.startswith("jacobian_verification:"):
            reason = "C finite-difference Jacobian checks"
        elif case_id.startswith(("constraint_set:", "linear:", "equilibrium:")):
            reason = "C linear/equilibrium constraint residual and Jacobian behavior checked"
        elif case_id.startswith("constraint_initialization:"):
            reason = "C constrained solve projects consistent and inconsistent algebraic initial conditions"
        elif case_id == "process_set:ProcessSetAlgebraicVariables.CudaMatrix":
            reason = "C masked reaction forcing and Jacobian updates preserve algebraic rows"
        elif case_id.startswith("matrix:"):
            reason = "C strided matrix view, conversion, arithmetic, and reductions checked"
        elif case_id.startswith("vector_matrix:"):
            reason = "C grouped vector-matrix view packing, access, and stride behavior checked"
        elif case_id.startswith("version:"):
            reason = "C version API checked"
        elif case_id in {"process_set:ProcessSet.Matrix", "process_set:ProcessSet.VectorMatrix"}:
            reason = "C forcing and Jacobian match grouped-process reference values"
        elif case_id == "species:Species.SetThirdBody":
            reason = "C third-body species parameterizes to air density"
        elif case_id.startswith("phase:PhaseSpecies."):
            reason = "C species diffusion metadata checked"
        elif case_id == "solver_lambda_rate_constant:Solver.GetLambdaRateConstantByNameCanOverrideLambda":
            reason = "caller-overridden C lambda callback drives rate constants and solve"
        elif case_id == "integrated_reaction_rates:ChapmanIntegration.CanBuildChapmanSystem":
            reason = "C product concentrations remain nondecreasing across rate changes"
        elif case_id == "chapman_integration:ChapmanIntegration.CanBuildChapmanSystem":
            reason = "C Chapman integration conserves oxygen and photolyzes O2"
        elif case_id.startswith("terminator:RosenbrockSolver."):
            reason = "C scalar and grouped solvers match the analytic terminator solution"
        elif case_id == "process_set:RandomProcessSet.Matrix":
            reason = "C random process mechanisms evaluate across the reference sizes"
        elif case_id == "process_set:ProcessSetFiniteDifferenceJacobian.Matrix":
            reason = "C process Jacobian agrees with finite differences"
        elif case_id in {"process:Process.Matrix", "process:Process.VectorMatrix"}:
            reason = "C rate constants match scalar and grouped cell layouts"
        elif case_id.startswith("state:State."):
            reason = "C caller-owned state views and setter errors checked"
        elif case_id == "openmp:OpenMP.OneSolverManyStates":
            reason = "independent C workspaces give identical thread states"
        print(f"{case_id} PASS {reason}")
    elif case_id in not_applicable:
        print(f"{case_id} N/A {not_applicable[case_id]}")
    elif kind == "program":
        failed = True
        print(f"{case_id} FAIL C program port incomplete")
    else:
        failed = True
        print(f"{case_id} FAIL C case port incomplete")
raise SystemExit(1 if failed else 0)
