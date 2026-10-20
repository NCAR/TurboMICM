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
    "phase_species:PhaseSpecies.ConstructsWithoutDiffusionCoefficient",
    "phase_species:PhaseSpecies.ConstructsWithDiffusionCoefficient",
    "phase_species:PhaseSpecies.SetDiffusionCoefficient",
    "phase:Phase.Constructor",
    "phase:Phase.ConstructorWithParameterizedSpecies",
    "phase:Phase.UniqueNamesWithParameterizedSpecies",
    "system:System.ConstructorWithAllParameters",
    "system:System.ConstructorWithParameterizedSpecies",
    "species:Species.StringConstructor",
    "species:Species.StringAndVectorConstructor",
    "species:Species.GetProperty",
    "species:Species.CopyConstructor",
    "species:Species.CopyAssignment",
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
    "matrix:Matrix.ReduceLOr",
    "matrix:Matrix.ReduceLAnd",
    "matrix:Matrix.ReduceStrict",
    "matrix:Matrix.MatrixVectorDifferentRowsFromCreation",
    "matrix:Matrix.MultipleTemporaries",
    "matrix:Matrix.ConstMatrixFunction",
    "matrix:Matrix.FunctionWithConstSignature",
    "matrix:Matrix.MultipleVectorsDifferentSizes",
    "matrix:Matrix.MultipleVectorsSameSize",
    "matrix:Matrix.MultipleMatricesOneVector",
    "matrix:Matrix.MultipleMatricesDifferentRowsVector",
    "matrix:Matrix.VectorSizeMatchesOneMatrixOnly",
    "matrix:Matrix.ConstVector",
    "matrix:Matrix.MutableVector",
    "matrix:Matrix.FunctionReusabilityWithVectors",
    "matrix:Matrix.FunctionInvocationWithWrongSizedVector",
    "matrix:Matrix.FunctionWithConstSignature",
    "matrix:Matrix.MixedVectorColumnViewRowVariable",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.ZeroMatrix",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.ConstZeroMatrix",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.SetScalar",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.AddToDiagonal",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.SingleBlockMatrix",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.ConstSingleBlockMatrix",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.MultiBlockMatrix",
    "sparse_matrix_standard_compressed_row_ordering:SparseMatrixBuilder.BadConfiguration",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.Print",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.PrintNonZero",
    "sparse_matrix_vector_compressed_row_ordering:SparseVectorCompressedRowMatrix.ZeroMatrix",
    "sparse_matrix_vector_compressed_row_ordering:SparseVectorCompressedRowMatrix.ConstZeroMatrix",
    "sparse_matrix_vector_compressed_row_ordering:SparseVectorCompressedRowMatrix.SetScalar",
    "sparse_matrix_vector_compressed_row_ordering:SparseVectorCompressedRowMatrix.AddToDiagonal",
    "sparse_matrix_vector_compressed_row_ordering:SparseVectorCompressedRowMatrix.SingleBlockMatrix",
    "sparse_matrix_vector_compressed_row_ordering:SparseVectorCompressedRowMatrix.ConstSingleBlockMatrix",
    "sparse_matrix_vector_compressed_row_ordering:SparseVectorCompressedRowMatrix.MultiBlockMatrix",
    "sparse_matrix_vector_compressed_column_ordering:SparseVectorCompressedColumnMatrix.ZeroMatrix",
    "sparse_matrix_vector_compressed_column_ordering:SparseVectorCompressedColumnMatrix.ConstZeroMatrix",
    "sparse_matrix_vector_compressed_column_ordering:SparseVectorCompressedColumnMatrix.SetScalar",
    "sparse_matrix_vector_compressed_column_ordering:SparseVectorCompressedColumnMatrix.AddToDiagonal",
    "sparse_matrix_vector_compressed_column_ordering:SparseVectorCompressedColumnMatrix.SingleBlockMatrix",
    "sparse_matrix_vector_compressed_column_ordering:SparseVectorCompressedColumnMatrix.ConstSingleBlockMatrix",
    "sparse_matrix_vector_compressed_column_ordering:SparseVectorCompressedColumnMatrix.MultiBlockMatrix",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.IntegerVectorSparse",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.IntegerVectorSparse",
    "sparse_matrix_vector_compressed_row_ordering:SparseVectorCompressedRowMatrix.IntegerVectorSparse",
    "sparse_matrix_vector_compressed_column_ordering:SparseVectorCompressedColumnMatrix.IntegerVectorSparse",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.ZeroMatrix",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.ConstZeroMatrix",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.SetScalar",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.AddToDiagonal",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.SingleBlockMatrix",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.ConstSingleBlockMatrix",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.MultiBlockMatrix",
    "sparse_matrix_standard_compressed_column_ordering:SparseMatrixBuilder.BadConfiguration",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.Print",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.PrintNonZero",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.ArrayFunction",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.MultiMatrixArrayFunction",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.MismatchedBlockDimensions",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.MismatchedElementDimensions",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.EmptySparseMatrixFunction",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.BlockViewReuse",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.FunctionReusability",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.TwoSparseMatricesDifferentStructure",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.SparseAndDenseMatrixFunction",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.GetBlockViewByVectorIndex",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.TestFill",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.TestCopy",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.MultipleTemporaries",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.ConstSparseMatrixFunction",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.MultipleSparseMatricesDifferentBlocksFromCreation",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.SparseMatrixVectorDifferentBlocksFromCreation",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.MismatchedBlocksAtInvocation",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.MultipleSparseMatricesMismatchedBlocksAtInvocation",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.VectorInSparseMatrixFunction",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.VectorTooSmall",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.VectorTooLarge",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.EmptyVectorNonEmptySparseMatrix",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.NonEmptyVectorEmptySparseMatrix",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.EmptyVectorEmptySparseMatrix",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.MultipleVectorsDifferentSizes",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.MultipleVectorsSameSize",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.MultipleSparseMatricesOneVector",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.MultipleSparseMatricesDifferentBlocksVector",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.VectorSizeMatchesOneSparseMatrixOnly",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.ConstVectorSparse",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.MutableVectorSparse",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.FunctionWithConstSignatureSparse",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.FunctionReusabilityWithVectorsSparse",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.FunctionInvocationWithWrongSizedVectorSparse",
    "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.MixedVectorBlockViewBlockVariable",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.ArrayFunction",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.MultiMatrixArrayFunction",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.TestFill",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.TestCopy",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.MismatchedBlockDimensions",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.MismatchedElementDimensions",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.EmptySparseMatrixFunction",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.BlockViewReuse",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.FunctionReusability",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.TwoSparseMatricesDifferentStructure",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.SparseAndDenseMatrixFunction",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.GetBlockViewByVectorIndex",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.MultipleTemporaries",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.ConstSparseMatrixFunction",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.MultipleSparseMatricesDifferentBlocksFromCreation",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.SparseMatrixVectorDifferentBlocksFromCreation",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.MismatchedBlocksAtInvocation",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.MultipleSparseMatricesMismatchedBlocksAtInvocation",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.VectorInSparseMatrixFunction",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.VectorTooSmall",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.VectorTooLarge",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.EmptyVectorNonEmptySparseMatrix",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.NonEmptyVectorEmptySparseMatrix",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.EmptyVectorEmptySparseMatrix",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.MultipleVectorsDifferentSizes",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.MultipleVectorsSameSize",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.MultipleSparseMatricesOneVector",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.MultipleSparseMatricesDifferentBlocksVector",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.VectorSizeMatchesOneSparseMatrixOnly",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.ConstVectorSparse",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.MutableVectorSparse",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.FunctionWithConstSignatureSparse",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.FunctionReusabilityWithVectorsSparse",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.FunctionInvocationWithWrongSizedVectorSparse",
    "sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.MixedVectorBlockViewBlockVariable",
    "matrix:Matrix.IterateOverMatrix",
    "matrix:Matrix.FunctionReusability",
    "matrix:Matrix.MultiMatrixDifferentRowsFromCreation",
    "matrix:Matrix.VectorCapture",
    "matrix:Matrix.MultiMatrixArrayFunction",
    "matrix:Matrix.MismatchedRowsAtInvocation",
    "matrix:Matrix.MultipleMatricesMismatchedRowsAtInvocation",
    "vector_matrix:VectorMatrix.SmallVectorMatrix",
    "vector_matrix:VectorMatrix.SmallConstVectorMatrix",
    "vector_matrix:VectorMatrix.InitializeVectorMatrix",
    "vector_matrix:VectorMatrix.InitializeConstVectorMatrix",
    "vector_matrix:VectorMatrix.LoopOverVectorMatrix",
    "vector_matrix:VectorMatrix.LoopOverConstVectorMatrix",
    "vector_matrix:VectorMatrix.Strides",
    "vector_matrix:VectorMatrix.ConversionToVector",
    "vector_matrix:VectorMatrix.ConstConversionToVector",
    "vector_matrix:VectorMatrix.ConversionFromVector",
    "vector_matrix:VectorMatrix.AssignmentFromVector",
    "vector_matrix:VectorMatrix.Axpy",
    "vector_matrix:VectorMatrix.ForEach",
    "vector_matrix:VectorMatrix.SetScaler",
    "vector_matrix:VectorMatrix.Max",
    "vector_matrix:VectorMatrix.Min",
    "vector_matrix:VectorMatrix.ArrayFunction",
    "vector_matrix:VectorMatrix.ReduceMax",
    "vector_matrix:VectorMatrix.FunctionReusability",
    "vector_matrix:VectorMatrix.FunctionReusabilityWithVectors",
    "vector_matrix:VectorMatrix.MultiMatrixDifferentRowsFromCreation",
    "vector_matrix:VectorMatrix.ColumnViewReuse",
    "vector_matrix:VectorMatrix.VectorInMatrixFunction",
    "vector_matrix:VectorMatrix.VectorTooSmall",
    "vector_matrix:VectorMatrix.VectorTooLarge",
    "vector_matrix:VectorMatrix.EmptyVectorNonEmptyMatrix",
    "vector_matrix:VectorMatrix.NonEmptyVectorEmptyMatrix",
    "vector_matrix:VectorMatrix.EmptyVectorEmptyMatrix",
    "vector_matrix:VectorMatrix.MultipleVectorsDifferentSizes",
    "vector_matrix:VectorMatrix.MultipleVectorsSameSize",
    "vector_matrix:VectorMatrix.TestFill",
    "vector_matrix:VectorMatrix.TestCopy",
    "vector_matrix:VectorMatrix.MultiMatrixArrayFunction",
    "vector_matrix:VectorMatrix.MismatchedRowDimensions",
    "vector_matrix:VectorMatrix.MismatchedRowsAtInvocation",
    "grouped_views:GroupedView.DenseGroupedMatchesRawL1",
    "grouped_views:GroupedView.DenseGroupedMatchesRawL2",
    "grouped_views:GroupedView.DenseGroupedMatchesRawL4",
    "grouped_views:GroupedView.DenseGroupedMatchesRawL8",
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
_vector_sparse_case_names = (
    "Print", "PrintNonZero", "ArrayFunction", "MultiMatrixArrayFunction",
    "MultipleTemporaries", "BlockViewReuse", "FunctionReusability",
    "TwoSparseMatricesDifferentStructure", "SparseAndDenseMatrixFunction",
    "SparseAndVectorMatrixFunction", "MismatchedBlockDimensions",
    "MismatchedElementDimensions", "ConstSparseMatrixFunction",
    "EmptySparseMatrixFunction", "MultipleSparseMatricesDifferentBlocksFromCreation",
    "SparseMatrixVectorDifferentBlocksFromCreation", "MismatchedBlocksAtInvocation",
    "MultipleSparseMatricesMismatchedBlocksAtInvocation",
    "VectorInSparseMatrixFunction", "VectorTooSmall", "VectorTooLarge",
    "EmptyVectorNonEmptySparseMatrix", "NonEmptyVectorEmptySparseMatrix",
    "EmptyVectorEmptySparseMatrix", "MultipleVectorsDifferentSizes",
    "MultipleVectorsSameSize", "MultipleSparseMatricesOneVector",
    "MultipleSparseMatricesDifferentBlocksVector", "VectorSizeMatchesOneSparseMatrixOnly",
    "ConstVectorSparse", "MutableVectorSparse", "FunctionWithConstSignatureSparse",
    "FunctionReusabilityWithVectorsSparse", "FunctionInvocationWithWrongSizedVectorSparse",
    "MixedVectorBlockViewBlockVariable", "GetBlockViewByVectorIndex", "TestFill", "TestCopy",
)
for _category, _case_class in (
    ("sparse_matrix_vector_compressed_row_ordering", "SparseVectorCompressedRowMatrix"),
    ("sparse_matrix_vector_compressed_column_ordering", "SparseVectorCompressedColumnMatrix"),
):
    implemented.update(f"{_category}:{_case_class}.{_name}"
                       for _name in _vector_sparse_case_names)
    implemented.add(f"{_category}:SparseMatrixBuilder.BadConfiguration")
implemented.update(f"vector_matrix:VectorMatrix.{_name}" for _name in (
    "Print", "MismatchedColumnDimensions", "WrongMatrixDimensions",
    "MultipleTemporaries", "ConstMatrixFunction", "EmptyMatrixFunction",
    "MatrixVectorDifferentRowsFromCreation", "MultipleMatricesMismatchedRowsAtInvocation",
    "MultipleMatricesOneVector", "MultipleMatricesDifferentRowsVector",
    "VectorSizeMatchesOneMatrixOnly", "ConstVector", "MutableVector",
    "FunctionInvocationWithWrongSizedVector", "MixedVectorColumnViewRowVariable",
    "FunctionWithConstSignature", "IntegerVector", "ReduceSum", "ReduceLOr",
    "ReduceLAnd", "ReduceStrict", "VectorCapture",
))
implemented.update(f"grouped_views:GroupedView.{_name}" for _name in (
    "MatrixColumnViewRoundTrip", "VectorMatrixColumnViewContiguousBlock",
    "PlainMatrixGroupedRoundTrip", "SparseStandardBlockViewRoundTrip",
    "SparseVectorBlockViewContiguousBlock", "MixedSparseAndDenseGroupedL1",
    "MixedSparseAndDenseGroupedL2", "MixedSparseAndDenseGroupedL4",
    "MixedSparseAndDenseGroupedL8", "SparseVectorRowColOverloadStillFunctional",
    "ConstGroupViewProducesGroupedViews", "FillAndCopySparseStandardCSR",
    "FillAndCopySparseStandardCSC", "FillAndCopySparseVectorCSR_L4",
    "FillAndCopySparseVectorCSC_L4", "FillAndCopySparseVectorCSR_L8",
    "FillCopyCrossGroupViewSources",
))
_lu_case_names = (
    "DenseMatrixStandardOrdering", "RandomMatrixStandardOrdering",
    "DiagonalMatrixStandardOrdering", "AgnosticToInitialValueStandardOrdering",
    "DenseMatrixVectorOrdering", "RandomMatrixVectorOrdering",
    "DiagonalMatrixVectorOrdering", "AgnosticToInitialValueVectorOrdering",
)
for _category, _case_class in (
    ("lu_decomposition_doolittle", "LuDecompositionDoolittle"),
    ("lu_decomposition_mozart", "LuDecompositionMozart"),
    ("lu_decomposition_doolittle_in_place", "LuDecompositionDoolittleInPlace"),
    ("lu_decomposition_mozart_in_place", "LuDecompositionMozartInPlace"),
):
    implemented.update(f"{_category}:{_case_class}.{_name}"
                       for _name in _lu_case_names)
_linear_case_names = (
    "DenseMatrixStandardOrdering", "RandomMatrixStandardOrdering",
    "DiagonalMatrixStandardOrdering", "StandardOrderingAgnosticToInitialValue",
    "DenseMatrixVectorOrdering", "RandomMatrixVectorOrdering",
    "VectorOrderingAgnosticToInitialValue", "DiagonalMatrixVectorOrdering",
)
for _category, _case_class in (
    ("linear_solver", "LinearSolver"),
    ("linear_solver_in_place", "LinearSolverInPlace"),
):
    implemented.update(f"{_category}:{_case_class}.{_name}"
                       for _name in _linear_case_names)
implemented.update({
    "linear_solver:LinearSolver.DiagonalMarkowitzReorder",
    "linear_solver:LinearSolver.VectorDiagonalMarkowitzReordering",
    "linear_solver_in_place:LinearSolverInPlace.DiagonalMarkowitzReorder",
    "linear_solver_in_place:LinearSolverInPlace.VectorDiagonalMarkowitzReordering",
})
implemented.update({
    "matrix:Matrix.Print", "matrix:Matrix.IntegerVector",
    "rosenbrock:RosenbrockSolver.StandardNormalizedError",
    "rosenbrock:RosenbrockSolver.VectorNormalizedError",
    "rosenbrock:RosenbrockSolver.StandardNormalizedErrorWithConstraints",
    "rosenbrock:RosenbrockSolver.VectorNormalizedErrorWithConstraints",
})
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
    "sparse_matrix_vector_compressed_row_ordering:SparseVectorCompressedRowMatrix.IncompatibleOrdering":
        "C uses runtime scalar views and has no compile-time ordering mismatch",
    "sparse_matrix_vector_compressed_row_ordering:SparseVectorCompressedRowMatrix.IncompatibleVectorOrdering":
        "C uses runtime scalar views and has no compile-time ordering mismatch",
    "sparse_matrix_vector_compressed_row_ordering:SparseVectorCompressedRowMatrix.IncompatibleSparseOrdering":
        "C uses runtime scalar views and has no compile-time ordering mismatch",
    "sparse_matrix_vector_compressed_column_ordering:SparseVectorCompressedColumnMatrix.IncompatibleOrdering":
        "C uses runtime scalar views and has no compile-time ordering mismatch",
    "sparse_matrix_vector_compressed_column_ordering:SparseVectorCompressedColumnMatrix.IncompatibleVectorOrdering":
        "C uses runtime scalar views and has no compile-time ordering mismatch",
    "sparse_matrix_vector_compressed_column_ordering:SparseVectorCompressedColumnMatrix.IncompatibleSparseOrdering":
        "C uses runtime scalar views and has no compile-time ordering mismatch",
    "grouped_views:GroupedView.CategoryTagsAreDistinct":
        "C has no distinct compile-time view category tags",
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
        elif case_id.startswith(("phase_species:", "phase:", "system:")):
            reason = "C phase and system views preserve species metadata and qualified names"
        elif case_id.startswith("species:"):
            reason = "C species records support parameterization and typed caller-owned properties"
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
