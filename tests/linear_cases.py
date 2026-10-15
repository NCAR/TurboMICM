"""LU reconstruction, solve residuals, fill, ordering, and scratch initialization.

The reference uses implementation-defined C++ random distributions. The port uses
Python's fixed seed 1 with the same lognormal parameters, dimensions, and bounds.
Original assertions: Copyright (C) 2023-2026 University Corporation for Atmospheric Research.
SPDX-License-Identifier: Apache-2.0
"""
import math
import random
from api import *
from harness import case, equal, near


def matrix_model(n,pattern):
    m = Model([f'x{i}' for i in range(n)])
    flat = (C.c_ubyte*(n*n))(*(int((i,j) in pattern) for i in range(n) for j in range(n)))
    equal(lib.chem_model_sparsity(m.p,flat,n),0)
    return m.finalize()


def random_matrix(n,cells,diagonal=False,repeat=False):
    pattern_rng,values_rng = random.Random(1),random.Random(1)
    pattern = {(i,j) for i in range(n) for j in range(n) if i == j or (not diagonal and pattern_rng.randrange(2))}
    matrices = [[[0.0]*n for i in range(n)] for c in range(cells)]
    for i,j in sorted(pattern):
        values = [values_rng.lognormvariate(-2,4 if diagonal else 2) for _ in range(1 if repeat else cells)]
        for c in range(cells): matrices[c][i][j] = values[0 if repeat else c]
    return pattern,matrices,values_rng


def flatten(m,matrices,stride,initial=0):
    a = arr([initial]*(m.data.nnz*stride)); n = m.data.nvariables
    for i in range(n):
        for j in range(n):
            at = m.data.index[i*n+j]
            if at >= 0:
                for c,matrix in enumerate(matrices): a[at*stride+c] = matrix[i][j]
    return a


def check_lu(m,matrices,stride,lower,upper,tolerance,combined=False):
    n = m.data.nvariables
    def val(buffer,i,j,c):
        at = m.data.index[i*n+j]
        return buffer[at*stride+c] if at >= 0 else 0.0
    for c,a in enumerate(matrices):
        for i in range(n):
            for j in range(n):
                result = 0.0
                for k in range(n):
                    l = (1.0 if k == i else val(lower,i,k,c)) if combined else val(lower,i,k,c)
                    if k > i: l = 0.0
                    u = val(upper,k,j,c) if k <= j else 0.0
                    result += l*u
                near(result,a[i][j],tolerance)
                if not combined:
                    if j > i: equal(val(lower,i,j,c),0)
                    if i > j: equal(val(upper,i,j,c),0)


def factor_case(identifier):
    family,name = identifier.split(':'); inplace = '_in_place' in family
    doolittle = 'doolittle' in family
    if 'DenseMatrix' in name:
        matrices = [[[2,-1,-2],[-4,6,3],[-4,-2,8]]]; n,cells,tolerance = 3,1,1e-10
        pattern = {(i,j) for i in range(n) for j in range(n)}
        initials = [0]
    else:
        diagonal = 'DiagonalMatrix' in name
        extreme = 'Agnostic' in name
        n,cells = (6,5) if diagonal else ((5 if inplace else 10),5) if extreme else (10,5)
        pattern,matrices,_ = random_matrix(n,cells,diagonal,repeat=inplace)
        tolerance = 1e-10 if diagonal else (1e-12 if inplace else 1e-9)
        initials = [-math.inf,-1,0,1,math.inf] if extreme else [0]
    m = matrix_model(n,pattern)
    try:
        for group in ((1,2,3,4) if 'Vector' in name else (1,)):
            stride = ((cells+group-1)//group)*group
            for initial in initials:
                a = flatten(m,matrices,stride,initial)
                if inplace:
                    factor = lib.chem_lu_doolittle if doolittle else lib.chem_lu_mozart
                    factor(m.p,cells,stride,a)
                    check_lu(m,matrices,stride,a,a,tolerance,True)
                else:
                    lower,upper = arr([initial]*(m.data.nnz*stride)),arr([initial]*(m.data.nnz*stride))
                    equal(lib.chem_lu_separate(m.p,cells,stride,a,lower,upper,doolittle),0)
                    check_lu(m,matrices,stride,lower,upper,tolerance)
    finally: m.close()


def markowitz_case():
    n = 50; rng = random.Random(1)
    pattern = {(i,j) for i in range(n) for j in range(n) if i == j or rng.randrange(2)}
    flat = (C.c_ubyte*(n*n))(*(int((i,j) in pattern) for i in range(n) for j in range(n)))
    perm = (Z*n)(); equal(lib.chem_markowitz(n,flat,perm),0)
    equal(sorted(perm),list(range(n)))
    reordered = {(i,j) for i in range(n) for j in range(n) if (perm[i],perm[j]) in pattern}
    original,moved = matrix_model(n,pattern),matrix_model(n,reordered)
    try: assert original.data.nnz > moved.data.nnz
    finally: original.close(); moved.close()


def linear_case(identifier):
    family,name = identifier.split(':'); inplace = '_in_place' in family
    if 'Markowitz' in name: return markowitz_case()
    dense = 'DenseMatrix' in name; diagonal = 'DiagonalMatrix' in name; extreme = 'Agnostic' in name
    n,cells = (3,1) if dense else ((6,5) if diagonal else ((30,5) if extreme else (10,5)))
    if dense:
        matrices = [[[2,-1,-2],[-4,6,3],[-4,-2,8]]]; pattern = {(i,j) for i in range(3) for j in range(3)}
        rhs = [[23,42,9]]
    else:
        pattern,matrices,rng = random_matrix(n,cells,diagonal)
        rhs = [[0.0]*n for c in range(cells)]
        for i in range(n):
            for c in range(cells): rhs[c][i] = rng.lognormvariate(-2,4 if diagonal else 2)
    tolerance = 1e-5 if dense or diagonal else (2e-6 if extreme else 1e-6)
    initials = [-math.inf,-1,0,1,math.inf] if extreme else [0]
    m = matrix_model(n,pattern)
    try:
        for group in ((1,2,3,4) if 'Vector' in name else (1,)):
            stride = ((cells+group-1)//group)*group
            for initial in initials:
                a = flatten(m,matrices,stride)
                x = arr([initial]*(n*stride))
                for i in range(n):
                    for c in range(cells): x[i*stride+c] = rhs[c][i]
                if inplace:
                    lib.chem_lu_mozart(m.p,cells,stride,a); lib.chem_linear_solve(m.p,cells,stride,a,x)
                else:
                    lower,upper = arr([initial]*(m.data.nnz*stride)),arr([initial]*(m.data.nnz*stride))
                    equal(lib.chem_lu_separate(m.p,cells,stride,a,lower,upper,1),0)
                    lib.chem_linear_solve_separate(m.p,cells,stride,lower,upper,x)
                for c in range(cells):
                    for i in range(n): near(sum(matrices[c][i][j]*x[j*stride+c] for j in range(n)),rhs[c][i],tolerance)
    finally: m.close()


for _line in (ROOT/'case-list.txt').read_text().splitlines():
    if _line.startswith('lu_decomposition_'): case(_line.split()[0])(factor_case)
    elif _line.startswith(('linear_solver:','linear_solver_in_place:')): case(_line.split()[0])(linear_case)
