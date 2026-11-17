"""Check that the generic diagonal Markowitz permutation reduces LU fill."""
from pathlib import Path
import random
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
from generate import markowitz_pattern


def fill(pattern):
    edges = [set(row) for row in pattern]
    n = len(edges)
    for pivot in range(n):
        above = {j for j in edges[pivot] if j > pivot}
        for row in range(pivot+1,n):
            if pivot in edges[row]:
                edges[row].update(above)
    return sum(len(row) for row in edges)


class ReorderTest(unittest.TestCase):
    def test_reference_size_and_fill(self):
        rng = random.Random(4723)
        n = 50
        pattern = [{j for j in range(n) if i == j or rng.randrange(2)} for i in range(n)]
        order = markowitz_pattern(pattern)
        self.assertEqual(sorted(order), list(range(n)))
        inverse = {old: new for new, old in enumerate(order)}
        reordered = [{inverse[j] for j in pattern[old]} for old in order]
        self.assertLess(fill(reordered), fill(pattern))

if __name__ == '__main__':
    unittest.main()
