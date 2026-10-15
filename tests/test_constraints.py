"""Ports of the equilibrium, linear, and constraint-set reference cases."""
# Copyright (C) 2023-2026 University Corporation for Atmospheric Research
# SPDX-License-Identifier: Apache-2.0
from chemffi import *

CASES = {}

def case(suite, *names):
    def decorate(fn):
        for name in names:
            CASES[f'{suite}.{name}'] = fn
        return fn
    return decorate


def near(got, expected, tolerance=1e-10):
    assert abs(got - expected) <= tolerance, (got, expected, tolerance)


def pattern(m):
    return {(r, c) for r in range(m.n) for c in range(m.n)
            if lib.chem_constraint_nonzero(m.handle, r, c)}


def descriptor(m):
    return lib.chem_mechanism_data(m.handle).contents.constraints[0]


def finite_difference(m):
    rows = m.rows()
    jac = m.jacobian(True)
    for col in range(m.n):
        for cell in range(m.nc):
            x = rows[cell][col]
            h = 1e-8 * max(1, abs(x))
            m.y[col * m.nc + cell] = x + h
            plus = m.force(True)[cell]
            m.y[col * m.nc + cell] = x - h if x - h >= 0 else x
            minus = m.force(True)[cell]
            m.y[col * m.nc + cell] = x
            for row in range(m.n):
                fd = (plus[row] - minus[row]) / (2*h if x-h >= 0 else h)
                analytical = -jac[cell][row][col]
                near(analytical, fd, 1e-7 + 1e-7 * max(abs(fd), abs(analytical)))
                if (row, col) not in pattern(m):
                    near(fd, 0, 1e-6)


@case('LinearConstraint', 'Construction')
def linear_construction():
    c = linear(1, [(0, 1), (1, 1)], 1)
    c[0].name = b'A_B_conservation'
    with Mechanism(['A', 'B'], constraints=[c]) as m:
        d = descriptor(m)
        assert d.name == b'A_B_conservation' and d.constant == 1
        assert d.reactant_count == 2
        assert [d.reactants[i].species for i in range(2)] == [0, 1]


@case('LinearConstraint', 'AlgebraicSpecies')
def linear_algebraic():
    with Mechanism(['A', 'B', 'C'], constraints=[linear(2, [(0,1),(1,1),(2,1)],10)]) as m:
        assert descriptor(m).algebraic_species == 2
        assert lib.chem_is_algebraic(m.handle,2)


@case('LinearConstraint', 'WeightedTerms')
def linear_weighted_terms():
    c = linear(2, [(0,2),(1,3),(2,-1)],5)
    c[0].name = b'weighted_sum'
    with Mechanism(['A','B','C'], constraints=[c]) as m:
        d = descriptor(m)
        assert d.name == b'weighted_sum' and d.constant == 5 and d.reactant_count == 3
        assert [d.reactants[i].yield_ for i in range(3)] == [2,3,-1]


@case('LinearConstraint', 'ZeroConstant')
def linear_zero():
    c = linear(1, [(0,1),(1,-1)],0)
    c[0].name = b'A_equals_B'
    with Mechanism(['A','B'], constraints=[c]) as m:
        d = descriptor(m)
        assert d.name == b'A_equals_B' and d.constant == 0 and d.reactant_count == 2


@case('LinearConstraint', 'FractionalCoefficients')
def linear_fractional():
    c = linear(1, [(0,.5),(1,1.5)],2)
    c[0].name = b'fractional_conservation'
    with Mechanism(['A','B'], constraints=[c]) as m:
        d = descriptor(m)
        assert d.name == b'fractional_conservation' and d.constant == 2
        assert [d.reactants[i].yield_ for i in range(2)] == [.5,1.5]


@case('LinearConstraint', 'ResidualComputationThroughConstraintSet')
def linear_residual():
    with Mechanism(['A','B'], constraints=[linear(1,[(0,1),(1,1)],1)]) as m:
        for row, expected in [([.3,.7],0),([.5,.6],.1)]:
            m.set_y([row]); near(m.force(True)[0][1],expected)


@case('LinearConstraint', 'JacobianComputationThroughConstraintSet')
def linear_jacobian():
    with Mechanism(['A','B'], constraints=[linear(1,[(0,1),(1,1)],1)]) as m:
        m.set_y([[.3,.7]])
        for value in m.jacobian(True)[0][1]: near(value,-1)


@case('LinearConstraint', 'WeightedSumResidualAndJacobian')
def linear_weighted():
    with Mechanism(['A','B','C'], constraints=[linear(2,[(0,2),(1,3),(2,-1)],5)]) as m:
        m.set_y([[1,2,3]]); near(m.force(True)[0][2],0)
        for got, want in zip(m.jacobian(True)[0][2],[-2,-3,1]): near(got,want)


@case('LinearConstraint', 'ThreeSpeciesConservationResidual')
def linear_three():
    with Mechanism(['A','B','C'], constraints=[linear(2,[(0,1),(1,1),(2,1)],10)]) as m:
        for row, expected in [([2,3,5],0),([2,3,6],1)]:
            m.set_y([row]); near(m.force(True)[0][2],expected)


@case('LinearConstraint', 'ZeroConstantResidual')
def linear_zero_residual():
    with Mechanism(['A','B'], constraints=[linear(1,[(0,1),(1,-1)],0)]) as m:
        for row, expected in [([.5,.5],0),([.6,.5],.1)]:
            m.set_y([row]); near(m.force(True)[0][1],expected)


@case('LinearConstraint', 'FractionalCoefficientsResidualAndJacobian')
def linear_fractional_residual():
    with Mechanism(['A','B'], constraints=[linear(1,[(0,.5),(1,1.5)],2)]) as m:
        m.set_y([[1,1]]); near(m.force(True)[0][1],0)
        for got, want in zip(m.jacobian(True)[0][1],[-.5,-1.5]): near(got,want)


@case('LinearConstraint', 'JacobianIndependentOfConcentrations')
def linear_constant_jacobian():
    with Mechanism(['A','B'], constraints=[linear(1,[(0,2),(1,3)],1)]) as m:
        for row in [[.1,.2],[10,20]]:
            m.set_y([row])
            for got, want in zip(m.jacobian(True)[0][1],[-2,-3]): near(got,want)


@case('LinearConstraint', 'FiniteDifferenceJacobianSimpleConservation')
def linear_fd_simple():
    with Mechanism(['A','B'], constraints=[linear(1,[(0,1),(1,1)],1)], cells=2) as m:
        m.set_y([[.3,.7],[.8,.2]]); finite_difference(m)


@case('LinearConstraint', 'FiniteDifferenceJacobianWeightedSum')
def linear_fd_weighted():
    with Mechanism(['A','B','C'], constraints=[linear(2,[(0,2),(1,3),(2,-1)],5)]) as m:
        m.set_y([[1.5,.8,2]]); finite_difference(m)


def eq_metadata(reactants, products, row, name, n=3):
    c = equilibrium(row, reactants, products, 1000, -2400)
    c[0].name = name
    with Mechanism(['A','B','AB'][:n], constraints=[c]) as m:
        d = descriptor(m)
        assert d.name == name
        assert d.algebraic_species == row
        assert d.reactant_count == len(reactants) and d.product_count == len(products)
        assert [(d.reactants[i].species,d.reactants[i].yield_) for i in range(d.reactant_count)] == reactants
        assert [(d.products[i].species,d.products[i].yield_) for i in range(d.product_count)] == products


@case('EquilibriumConstraint', 'Construction', 'AlgebraicSpecies')
def eq_construction():
    eq_metadata([(0,1),(1,1)],[(2,1)],2,b'A_B_equilibrium')


@case('EquilibriumConstraint', 'SingleReactantSingleProduct')
def eq_single():
    eq_metadata([(0,1)],[(1,1)],1,b'A_B_simple',2)


@case('EquilibriumConstraint', 'MultipleReactantsAndProducts')
def eq_multiple():
    eq_metadata([(0,2)],[(1,1),(2,1)],1,b'dissociation')


@case('EquilibriumConstraint', 'InvalidEquilibriumConstant')
def eq_invalid_constant():
    for k in [-1,0]:
        with Mechanism(['A','B'], constraints=[equilibrium(1,[(0,1)],[(1,1)],k,-2400)],expected_status=1): pass


@case('EquilibriumConstraint', 'EmptyReactantsThrows')
def eq_empty_reactants():
    with Mechanism(['A','B'], constraints=[equilibrium(1,[],[(1,1)])],expected_status=1): pass


@case('EquilibriumConstraint', 'EmptyProductsThrows')
def eq_empty_products():
    with Mechanism(['A','B'], constraints=[equilibrium(1,[(0,1)],[])],expected_status=1): pass


@case('EquilibriumConstraint', 'InvalidStoichiometryThrows')
def eq_invalid_stoich():
    for rr, pp in [([(0,0)],[(1,1)]), ([(0,-1)],[(1,1)]), ([(0,1)],[(1,0)]), ([(0,1)],[(1,-2)])]:
        with Mechanism(['A','B'],constraints=[equilibrium(1,rr,pp)],expected_status=1): pass


@case('EquilibriumConstraint', 'ResidualComputationThroughConstraintSet')
def eq_residual():
    with Mechanism(['A','B','AB'],constraints=[equilibrium(2,[(0,1),(1,1)],[(2,1)],1000)]) as m:
        for row, expected in [([.001,.001,.001],0),([.02,.03,.05],.55)]:
            m.set_y([row]); near(m.force(True)[0][2],expected)


@case('EquilibriumConstraint', 'JacobianComputationThroughConstraintSet')
def eq_jacobian():
    with Mechanism(['A','B','AB'],constraints=[equilibrium(2,[(0,1),(1,1)],[(2,1)],1000)]) as m:
        m.set_y([[.01,.02,.05]])
        for got, want in zip(m.jacobian(True)[0][2],[-20,-10,1]): near(got,want,1e-8)


@case('EquilibriumConstraint', 'ComplexStoichiometryResidual')
def eq_complex():
    with Mechanism(['A','B','C'],constraints=[equilibrium(1,[(0,2)],[(1,1),(2,1)],100)]) as m:
        m.set_y([[.1,.5,2]]); near(m.force(True)[0][1],0)


@case('EquilibriumConstraint', 'FiniteDifferenceJacobianSimple')
def eq_fd_simple():
    with Mechanism(['A','B','AB'],constraints=[equilibrium(2,[(0,1),(1,1)],[(2,1)],1000)],cells=2) as m:
        m.set_y([[.02,.03,.05],[.1,.2,5]]); finite_difference(m)


@case('EquilibriumConstraint', 'FiniteDifferenceJacobianComplexStoichiometry')
def eq_fd_complex():
    with Mechanism(['A','B','C'],constraints=[equilibrium(1,[(0,2)],[(1,1),(2,1)],100)]) as m:
        m.set_y([[.15,.8,1.5]]); finite_difference(m)

@case('ConstraintSet', 'Construction')
def set_construction():
    with Mechanism(['A','B','AB'],constraints=[equilibrium(2,[(0,1),(1,1)],[(2,1)],.033,-24000)]) as m:
        assert lib.chem_mechanism_data(m.handle).contents.constraint_count == 1


@case('ConstraintSet', 'ReplaceStateRowsMapsToAlgebraicSpecies')
def set_rows():
    with Mechanism(['A','B','C'],constraints=[equilibrium(2,[(1,1)],[(2,1)],.033,-24000)]) as m:
        assert [lib.chem_is_algebraic(m.handle,i) for i in range(3)] == [0,0,1]
        assert pattern(m) == {(2,1),(2,2)}


@case('ConstraintSet', 'NonZeroJacobianElements')
def set_nonzero():
    with Mechanism(['A','B','AB'],constraints=[equilibrium(2,[(0,1),(1,1)],[(2,1)])]) as m:
        assert pattern(m) == {(2,0),(2,1),(2,2)}


@case('ConstraintSet', 'MultipleConstraints')
def set_multiple():
    cs = [equilibrium(2,[(0,1),(1,1)],[(2,1)]),equilibrium(4,[(3,1)],[(4,1)])]
    with Mechanism(['A','B','AB','C','D'],constraints=cs) as m:
        assert lib.chem_mechanism_data(m.handle).contents.constraint_count == 2
        assert pattern(m) == {(2,0),(2,1),(2,2),(4,3),(4,4)}


@case('ConstraintSet', 'AddForcingTerms')
def set_forcing():
    with Mechanism(['A','B','AB'],constraints=[equilibrium(2,[(0,1),(1,1)],[(2,1)])],cells=2) as m:
        m.set_y([[.2,.4,.6],[.1,.3,.7]])
        for row, expected in zip(m.force(True),[-.59736,-.69901]): near(row[2],expected,1e-5)


@case('ConstraintSet', 'SubtractJacobianTerms')
def set_jacobian():
    with Mechanism(['A','B','AB'],constraints=[equilibrium(2,[(0,1),(1,1)],[(2,1)])]) as m:
        m.set_y([[.01,.02,.05]])
        for got,want in zip(m.jacobian(True)[0][2],[-.00066,-.00033,1]): near(got,want)


@case('ConstraintSet', 'EmptyConstraintSet')
def set_empty():
    with Mechanism(['A','B','C']) as m:
        assert lib.chem_mechanism_data(m.handle).contents.constraint_count == 0
        assert pattern(m) == set()
        m.set_y([[.1,.2,.3]])
        forcing=(D*4)(1,1,1,1)
        lib.chem_constraint_forcing(m.state,m.y,forcing)
        assert list(forcing)==[1,1,1,1]


@case('ConstraintSet', 'UnknownSpeciesThrows')
def set_unknown():
    with Mechanism(['A','B'],constraints=[equilibrium(4,[(2,1),(3,1)],[(4,1)])],expected_status=1): pass


@case('ConstraintSet', 'ThreeDStateOneConstraint')
def set_three():
    with Mechanism(['X','Y','Z'],constraints=[equilibrium(1,[(0,1)],[(1,1)])],cells=2) as m:
        assert pattern(m) == {(1,0),(1,1)}
        m.set_y([[10,.2,.5],[10,.33,.3]])
        for cell,expected in enumerate([.13,0]): near(m.force(True)[cell][1],expected)
        for matrix in m.jacobian(True):
            near(matrix[1][0],-.033); near(matrix[1][1],1); near(matrix[2][2],0)


@case('ConstraintSet', 'FourDStateTwoConstraints')
def set_four():
    cs=[equilibrium(1,[(0,1)],[(1,1)]),equilibrium(0,[(2,1),(3,1)],[(0,1)])]
    with Mechanism(['A','B','C','D'],constraints=cs,cells=3) as m:
        assert pattern(m)=={(1,0),(1,1),(0,2),(0,3),(0,0)}
        rows=[[.33,.01089,10,1],[.33,.01089,5,1],[.33,.02,5,1]]
        m.set_y(rows)
        for i,(g0,g1) in enumerate([(0,0),(-.165,0),(-.165,-.00911)]):
            f=m.force(True)[i];near(f[0],g0);near(f[1],g1,1e-5)
        for i,j in enumerate(m.jacobian(True)):
            near(j[1][0],-.033);near(j[1][1],1)
            near(j[0][2],-.033*rows[i][3]);near(j[0][3],-.033*rows[i][2]);near(j[0][0],1)


@case('ConstraintSet', 'CoupledConstraintsSharedSpecies')
def set_coupled():
    cs=[equilibrium(1,[(0,1)],[(1,1)]),equilibrium(2,[(0,1)],[(2,1)])]
    with Mechanism(['A','B','C'],constraints=cs) as m:
        assert pattern(m)=={(1,0),(1,1),(2,0),(2,2)}
        m.set_y([[.1,.0033,.0033]])
        near(m.force(True)[0][1],0);near(m.force(True)[0][2],0)
        j=m.jacobian(True)[0]
        for row in [1,2]: near(j[row][0],-.033);near(j[row][row],1)


@case('ConstraintSet', 'VectorizedMatricesRespectGridCellIndexing')
def set_multicell():
    rows=[[.01,.02,.05],[.03,.01,.2],[.001,.002,.004]]
    with Mechanism(['A','B','AB'],constraints=[equilibrium(2,[(0,1),(1,1)],[(2,1)])],cells=3) as m:
        m.set_y(rows)
        for cell,row in enumerate(rows):
            near(m.force(True)[cell][2],.033*row[0]*row[1]-row[2],1e-9)
            j=m.jacobian(True)[cell]
            near(j[2][0],-.033*row[1],1e-12);near(j[2][1],-.033*row[0],1e-12);near(j[2][2],1,1e-12)


@case('ConstraintSet', 'VectorMatrix1','VectorMatrix2','VectorMatrix3')
def set_vector():
    # Each reference vector-layout case invokes these four policy checks.
    set_construction();set_nonzero();set_forcing();set_jacobian()
