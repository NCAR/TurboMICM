"""Structural checks complement the generated-versus-generic numerical tests."""
import copy
import unittest

import configure
import specialize


class GeneratorTests(unittest.TestCase):
    def test_species_and_fill(self):
        mechanism = configure.read('mechanisms/ts1.json')
        plan = specialize.bind(mechanism)
        self.assertEqual((plan.n, len(plan.reactions), plan.nnz), (209, 547, 2381))
        for row in range(plan.n):
            self.assertIn((row, row), plan.indices)
        for diagonal, divides, updates in plan.pivots:
            for index in [diagonal] + divides + [i for op in updates for i in op]:
                self.assertLess(index, plan.nnz)

    def test_edit_reaction_and_parameter(self):
        mechanism = configure.read('mechanisms/ts1.json')
        edited = copy.deepcopy(mechanism)
        reaction = copy.deepcopy(edited.reactions[0])
        reaction['reactants'] = [0, 1, 0]
        reaction['products'] = [(2, 0.75), (2, 0.1)]
        reaction['rate']['A'] = 1.23456789e-12
        edited.reactions.append(reaction)
        plan = specialize.bind(edited)
        self.assertEqual(len(plan.reactions), 548)
        emitted = specialize.emit(edited, 'edited', width=0)
        self.assertIn('1.23456789e-12', emitted)
        self.assertIn('.nreactions=548', emitted)
        self.assertNotEqual(plan.indices, specialize.bind(mechanism).indices)

    def test_determinism_and_parameter_extremes(self):
        mechanism = configure.read('mechanisms/tutorial/robertson/config.json')
        for batch, split, algebra, width in [(1, 1, 1, 0), (64, 7, 64, 3)]:
            first = specialize.emit(mechanism, 'fixture', batch, split, algebra, width)
            self.assertEqual(first, specialize.emit(mechanism, 'fixture', batch, split, algebra, width))
            self.assertIn('const ck_specialized fixture', first)

    def test_isolated_species_cost(self):
        self.assertEqual(specialize.markowitz([set(), {1}, {1, 2}]), [1, 2, 0])


if __name__ == '__main__':
    unittest.main()
