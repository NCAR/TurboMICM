"""Mechanism-independent symbolic ordering, fill, and arithmetic schedule.

All lists preserve the generic library's arithmetic order. The unsigned-cost
rule matches MICM's Markowitz ordering even for isolated species.
"""
from dataclasses import dataclass


@dataclass
class Plan:
    variables: list
    variable: dict
    pattern: set
    flat: dict
    jobs: list

    @classmethod
    def create(cls, config):
        variables = [i for i, sp in enumerate(config['species'])
                     if sp.get('parameterization') != 'air density']
        if not variables:
            raise ValueError('A solver needs at least one variable species')
        variable = {sp: i for i, sp in enumerate(variables)}
        graph = set()
        for rx in config['reactions']:
            rs = [variable[s] for s in rx['reactants'] if s in variable]
            ps = [variable[s] for s, _ in rx['products'] if s in variable]
            graph.update((row, col) for col in rs for row in rs + ps)
        alive = set(range(len(variables)))
        permutation = []
        while alive:
            def cost(v):
                row = sum((v, j) in graph for j in alive)
                col = sum((i, v) in graph for i in alive)
                # Only relative ordering of underflowed costs matters.
                return ((row - 1) * (col - 1)) & ((1 << 64) - 1)
            pivot = min(alive, key=lambda v: (cost(v), v))
            permutation.append(pivot)
            alive.remove(pivot)
            rows = [i for i in alive if (i, pivot) in graph]
            cols = [j for j in alive if (pivot, j) in graph]
            graph.update((i, j) for i in rows for j in cols)
        variables = [variables[i] for i in permutation]
        variable = {sp: i for i, sp in enumerate(variables)}
        pattern = {(i, i) for i in range(len(variables))}
        jobs = []
        for r, rx in enumerate(config['reactions']):
            rs = [variable[s] for s in rx['reactants'] if s in variable]
            ps = [(variable[s], y) for s, y in rx['products'] if s in variable]
            for col in rs:
                pattern.update((row, col) for row in rs + [s for s, _ in ps])
                deps = rs.copy()
                deps.remove(col)
                jobs.append((col, r, deps, ps))
        for i in range(len(variables)):
            rows = [j for j in range(i + 1, len(variables)) if (j, i) in pattern]
            cols = [k for k in range(i + 1, len(variables)) if (i, k) in pattern]
            pattern.update((j, k) for k in cols for j in rows)
        return cls(variables, variable, pattern,
                   {pos: i for i, pos in enumerate(sorted(pattern))},
                   sorted(jobs, key=lambda j: (j[0], j[1])))

    def elimination(self):
        """Yield each pivot and its ordered multiplier/update operations."""
        n = len(self.variables)
        for i in range(n):
            rows = [j for j in range(i + 1, n) if (j, i) in self.pattern]
            cols = [k for k in range(i + 1, n) if (i, k) in self.pattern]
            yield self.flat[i, i], [self.flat[j, i] for j in rows], [
                (self.flat[j, k], self.flat[j, i], self.flat[i, k])
                for k in cols for j in rows]
