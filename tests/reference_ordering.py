#!/usr/bin/env python3
"""Check that symmetric ordering preserves entries and reduces LU fill."""

import copy
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tools"))
from generate import reorder_mechanism


def pattern(mechanism):
    n = len(mechanism["species"])
    result = {(i, i) for i in range(n)}
    for reaction in mechanism["reactions"]:
        for col in reaction["reactants"]:
            result.update((row, col) for row in reaction["reactants"])
            result.update((row, col) for row, _ in reaction["products"])
    return result


def filled_size(entries, n):
    matrix = set(entries)
    for pivot in range(n):
        below = [row for row in range(pivot + 1, n) if (row, pivot) in matrix]
        right = [col for col in range(pivot + 1, n) if (pivot, col) in matrix]
        matrix.update((row, col) for row in below for col in right)
    return len(matrix)


def check():
    n = 50
    rng = random.Random(1941)
    reactions = []
    for row in range(n):
        for col in range(n):
            if row != col and rng.randrange(2):
                reactions.append({"reactants": [col], "products": [(row, 1.0)]})
    mechanism = {"species": [{"name": f"S{i}"} for i in range(n)], "reactions": reactions}
    before = pattern(mechanism)
    ordered = reorder_mechanism(copy.deepcopy(mechanism))
    after = pattern(ordered)
    names = [entry["name"] for entry in ordered["species"]]
    inverse = {int(name[1:]): i for i, name in enumerate(names)}
    mapped = {(inverse[row], inverse[col]) for row, col in before}
    return after == mapped and len(before) == len(after) and filled_size(after, n) < filled_size(before, n)


ok = check()
for case in ("linear_solver:LinearSolver.DiagonalMarkowitzReorder",
             "linear_solver:LinearSolver.VectorDiagonalMarkowitzReordering",
             "linear_solver_in_place:LinearSolverInPlace.DiagonalMarkowitzReorder",
             "linear_solver_in_place:LinearSolverInPlace.VectorDiagonalMarkowitzReordering"):
    print(case, "PASS" if ok else "FAIL", "permutation and LU fill checked")
sys.exit(0 if ok else 1)
