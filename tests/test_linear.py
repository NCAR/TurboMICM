"""Reference LU reconstruction and linear-system residual assertions.

The C chemistry symbolic pattern accepts arbitrary external-model coordinates;
these tests use that generic interface to exercise the identical numerical
kernels used by both integrators and by the benchmark.
"""
from functools import partial
import math
from chemffi import *
from reference_random import Minstd, MT19937, booleans, bernoulli, lognormals

FACTOR = C.CFUNCTYPE(I, P, Z, DP)
SOLVE = C.CFUNCTYPE(None, P, Z, DP, DP)
bind('chem_factor_separate', I, P, Z, DP, DP, DP, Linear)
bind('chem_linear_solve_separate', None, P, Z, DP, DP, DP)
bind('chem_diagonal_markowitz', I, Z, Z, C.POINTER(MatrixIndex), C.POINTER(Z))
MISSING = Z(-1).value


class Pattern:
    def __init__(self, n, entries):
        self.n = n
        self.names = (C.c_char_p * n)(*[str(i).encode() for i in range(n)])
        self.coords = (MatrixIndex * len(entries))(*[MatrixIndex(i,j) for i,j in entries])
        self.model = ExternalModel(nonzero_count=len(entries), nonzeros=self.coords)
        self.data = Data(species_count=n, species_names=self.names,
                         external_model_count=1, external_models=C.pointer(self.model))
        self.handle = P()
        assert lib.chem_mechanism_create(C.byref(self.data), C.byref(self.handle)) == 0
        self.nnz = lib.chem_jacobian_size(self.handle)
        self.indices = [[lib.chem_jacobian_index(self.handle,i,j) for j in range(n)] for i in range(n)]

    def __enter__(self): return self
    def __exit__(self,*args): lib.chem_mechanism_free(self.handle)

    def pack(self, matrices, initial):
        nc = len(matrices)
        result = (D * (self.nnz*nc))(*([initial]*(self.nnz*nc)))
        for i in range(self.n):
            for j in range(self.n):
                k = self.indices[i][j]
                if k != MISSING:
                    for c, matrix in enumerate(matrices): result[k*nc+c] = matrix[i][j]
        return result

    def unpack(self, values, nc):
        return [[[0 if k == MISSING else values[k*nc+c] for k in row] for row in self.indices] for c in range(nc)]


def problem(kind, blocks, inplace, linear=False):
    if kind == 'Dense':
        return [(i,j) for i in range(3) for j in range(3)], [[[2,-1,-2],[-4,6,3],[-4,-2,8]]], [[23,42,9]]
    diagonal = kind == 'Diagonal'
    extreme = kind == 'Extreme'
    n = 6 if diagonal else (30 if linear and extreme else (5 if inplace and extreme else 10))
    bits = bernoulli() if linear and extreme else booleans()
    numbers = lognormals(4 if diagonal else 2, MT19937() if linear and extreme else Minstd())
    entries = [(i,j) for i in range(n) for j in range(n) if i == j or (not diagonal and next(bits))]
    matrices = [[[0.]*n for _ in range(n)] for _ in range(blocks)]
    for i,j in entries:
        if inplace and not linear and not diagonal:
            value = next(numbers)
            for matrix in matrices: matrix[i][j] = value
        else:
            for matrix in matrices: matrix[i][j] = next(numbers)
    rhs = [[0.]*n for _ in range(blocks)]
    for i in range(n):
        for c in range(blocks): rhs[c][i] = next(numbers)
    return entries, matrices, rhs


def near(a,b,tolerance,context):
    assert abs(a-b) <= tolerance, (context,a,b,tolerance)


def execute(kind, blocks, group, inplace, algorithm, linear, initial=0):
    entries, matrices, rhs = problem(kind,blocks,inplace,linear)
    n, blocks = len(matrices[0]), len(matrices)
    solver = getattr(lib,'chem_linear_'+algorithm)()
    with Pattern(n, entries) as pattern:
        # Runtime group sizes replace compile-time vector policies. Each group
        # has independent scratch, just as it does in the benchmark.
        for start in range(0,blocks,group or blocks):
            subset = matrices[start:start+(group or blocks)]
            nc = len(subset)
            matrix = pattern.pack(subset,initial)
            saved = list(matrix)
            if inplace:
                assert FACTOR(solver.factor)(pattern.handle,nc,matrix) == 0
                combined = pattern.unpack(matrix,nc)
                lower = [[[1. if i==j else (combined[c][i][j] if i>j else 0.) for j in range(n)] for i in range(n)] for c in range(nc)]
                upper = [[[combined[c][i][j] if i<=j else 0. for j in range(n)] for i in range(n)] for c in range(nc)]
            else:
                l = (D*(pattern.nnz*nc))(*([initial]*(pattern.nnz*nc)))
                u = (D*(pattern.nnz*nc))(*([initial]*(pattern.nnz*nc)))
                assert lib.chem_factor_separate(pattern.handle,nc,matrix,l,u,solver) == 0
                assert list(matrix) == saved
                lower,upper = pattern.unpack(l,nc),pattern.unpack(u,nc)
                for c in range(nc):
                    for i in range(n):
                        assert lower[c][i][i] == 1
                        for j in range(n):
                            assert i>=j or lower[c][i][j] == 0
                            assert j>=i or upper[c][i][j] == 0
            if linear:
                b = rhs[start:start+nc]
                x = (D*(n*nc))(*[b[c][i] for i in range(n) for c in range(nc)])
                if inplace: SOLVE(solver.solve)(pattern.handle,nc,matrix,x)
                else: lib.chem_linear_solve_separate(pattern.handle,nc,l,u,x)
                for c in range(nc):
                    for i in range(n):
                        result = sum(subset[c][i][j] * x[j*nc+c] for j in range(n))
                        if kind == 'Diagonal' and inplace:
                            # Reference EXPECT_REAL_EQ is the double four-ULP matcher.
                            tol = 4*max(math.ulp(result),math.ulp(b[c][i]))
                        else: tol = {'Dense':1e-5,'Random':1e-6,'Extreme':2e-6,'Diagonal':1e-5}[kind]
                        near(result,b[c][i],tol,(start+c,i))
            else:
                tol = 1e-10 if kind in ('Dense','Diagonal') else (1e-12 if inplace else 1e-9)
                for c in range(nc):
                    for i in range(n):
                        for j in range(n):
                            result = sum(lower[c][i][k]*upper[c][k][j] for k in range(n))
                            near(result,subset[c][i][j],tol,(start+c,i,j))


def check(kind, vector, inplace, algorithm='mozart', linear=False):
    for group in ([1,2,3,4] if vector else [0]):
        counts = [1] if kind == 'Dense' else ([1,5] if kind == 'Random' and not vector and not linear else [5])
        if linear and vector and kind == 'Extreme': counts = [1 if group==1 else (2 if group==2 else 5)]
        for blocks in counts:
            for initial in ([-math.inf,-1.,0.,1.,math.inf] if kind == 'Extreme' else [0.]):
                execute(kind,blocks,group,inplace,algorithm,linear,initial)


def markowitz(vector=False,inplace=False):
    n = 50
    bits = booleans()
    entries = [(i,j) for i in range(n) for j in range(n) if i == j or next(bits)]
    coords = (MatrixIndex*len(entries))(*[MatrixIndex(i,j) for i,j in entries])
    permutation = (Z*n)()
    assert lib.chem_diagonal_markowitz(n,len(entries),coords,permutation) == 0
    assert sorted(permutation) == list(range(n))
    inverse = {old:new for new,old in enumerate(permutation)}
    reordered = [(inverse[i],inverse[j]) for i,j in entries]
    with Pattern(n,entries) as original, Pattern(n,reordered) as ordered:
        # Separate factors include the duplicated unit diagonal.
        assert original.nnz + (0 if inplace else n) > ordered.nnz + (0 if inplace else n)
    # The output and fill quality must agree with the configuration reader's
    # independent implementation of the source's stable minimum-degree ordering.
    import sys
    from pathlib import Path
    sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'tools'))
    from configure import reorder
    mech = {'species':[{'name':str(i)} for i in range(n)],
            'reactions':[{'reactants':[j], 'products':[(i,1.)]} for i,j in entries]}
    expected = [int(s['name']) for s in reorder(mech)['species']]
    assert list(permutation) == expected


CASES = {}
for inplace in [False,True]:
    for algorithm in ['doolittle','mozart']:
        suite = 'LuDecomposition'+algorithm.title()+('InPlace' if inplace else '')
        for vector in [False,True]:
            ordering = 'VectorOrdering' if vector else 'StandardOrdering'
            for name,kind in [('DenseMatrix','Dense'),('RandomMatrix','Random'),('DiagonalMatrix','Diagonal'),('AgnosticToInitialValue','Extreme')]:
                CASES[suite+'.'+name+ordering] = partial(check,kind,vector,inplace,algorithm)
    suite = 'LinearSolver'+('InPlace' if inplace else '')
    for vector in [False,True]:
        ordering = 'VectorOrdering' if vector else 'StandardOrdering'
        for name,kind in [('DenseMatrix','Dense'),('RandomMatrix','Random'),('DiagonalMatrix','Diagonal')]:
            CASES[suite+'.'+name+ordering] = partial(check,kind,vector,inplace,linear=True)
        CASES[suite+'.'+('VectorOrdering' if vector else 'StandardOrdering')+'AgnosticToInitialValue'] = partial(check,'Extreme',vector,inplace,linear=True)
        CASES[suite+'.'+('VectorDiagonalMarkowitzReordering' if vector else 'DiagonalMarkowitzReorder')] = partial(markowitz,vector,inplace)
