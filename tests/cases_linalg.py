"""Sparse factorization and linear solve cases, with independent reconstruction."""
from ffi import *
from harness import group

LU_NAMES = '''DenseMatrixStandardOrdering RandomMatrixStandardOrdering DiagonalMatrixStandardOrdering
AgnosticToInitialValueStandardOrdering DenseMatrixVectorOrdering RandomMatrixVectorOrdering
DiagonalMatrixVectorOrdering AgnosticToInitialValueVectorOrdering'''
SOLVE_NAMES = '''DenseMatrixStandardOrdering RandomMatrixStandardOrdering DiagonalMatrixStandardOrdering
DiagonalMarkowitzReorder StandardOrderingAgnosticToInitialValue DenseMatrixVectorOrdering
RandomMatrixVectorOrdering VectorOrderingAgnosticToInitialValue DiagonalMatrixVectorOrdering
VectorDiagonalMarkowitzReordering'''


def sparse(pattern):
    data = np.ascontiguousarray(pattern,dtype=np.uint8)
    handle = sparse_create(len(data),data.ctypes.data_as(BP))
    assert handle
    return handle


def pack(handle, dense, initial=0):
    nc,n,_ = dense.shape
    result = np.full((sparse_size(handle),nc),initial,dtype=float)
    for r in range(n):
        for c in range(n):
            i = sparse_index(handle,r,c)
            if i >= 0:
                result[i,:] = dense[:,r,c]
    return result


def unpack(handle, values, n):
    nc = values.shape[1]
    result = np.zeros((nc,n,n))
    for r in range(n):
        for c in range(n):
            i = sparse_index(handle,r,c)
            if i >= 0:
                result[:,r,c] = values[i,:]
    return result


def markowitz_case():
    rng = np.random.default_rng(42)
    pattern = rng.integers(0,2,(50,50),dtype=np.uint8)
    np.fill_diagonal(pattern,1)
    perm = ints([0]*50)
    require(markowitz(50,pattern.ctypes.data_as(BP),perm))
    order = list(perm)
    assert sorted(order) == list(range(50))
    # Independent dense elimination with stable original-index tie breaking.
    live,expected,p = list(range(50)),[],pattern.copy()
    while live:
        costs = [(int(p[v,live].sum())-1)*(int(p[live,v].sum())-1) for v in live]
        pivot = live[costs.index(min(costs))]
        expected.append(pivot); live.remove(pivot)
        for i in live:
            if p[i,pivot]:
                for j in live:
                    if p[pivot,j]: p[i,j] = 1
    assert order == expected
    before,after = sparse(pattern),sparse(pattern[np.ix_(order,order)])
    try:
        assert sparse_size(after) <= sparse_size(before)
    finally:
        sparse_free(before); sparse_free(after)


@group('linear_solver:LinearSolver.', SOLVE_NAMES)
@group('linear_solver_in_place:LinearSolverInPlace.', SOLVE_NAMES)
@group('lu_decomposition_doolittle:LuDecompositionDoolittle.', LU_NAMES)
@group('lu_decomposition_mozart:LuDecompositionMozart.', LU_NAMES)
@group('lu_decomposition_doolittle_in_place:LuDecompositionDoolittleInPlace.', LU_NAMES)
@group('lu_decomposition_mozart_in_place:LuDecompositionMozartInPlace.', LU_NAMES)
def linalg(id_):
    if 'Markowitz' in id_:
        markowitz_case()
        return
    linear = id_.startswith('linear_')
    inplace = '_in_place:' in id_
    algorithm = 0 if 'mozart' in id_ else 1
    rng = np.random.default_rng(1)
    cell_counts = (1,2,3,4,5) if 'Vector' in id_ else (5,)
    initials = (-np.inf,-1,0,1,np.inf) if 'InitialValue' in id_ else (0,)
    for nc in cell_counts:
        if '.Dense' in id_:
            base = np.array([[2.,-1.,-2.],[-4.,6.,3.],[-4.,-2.,8.]])
            matrices = np.broadcast_to(base,(nc,3,3)).copy()
            rhs = np.broadcast_to([23.,42.,9.],(nc,3)).T.copy()
            pattern = np.ones((3,3),dtype=np.uint8)
        elif '.Diagonal' in id_:
            pattern = np.eye(6,dtype=np.uint8)
            matrices = np.zeros((nc,6,6))
            for c in range(nc): np.fill_diagonal(matrices[c],rng.lognormal(-2,4,6))
            rhs = rng.lognormal(-2,4,(6,nc))
        else:
            n = 5 if 'InitialValue' in id_ else 10
            pattern = rng.integers(0,2,(n,n),dtype=np.uint8)
            np.fill_diagonal(pattern,1)
            base = rng.lognormal(-2,2,(n,n))*pattern
            matrices = np.broadcast_to(base,(nc,n,n)).copy()
            rhs = rng.lognormal(-2,2,(n,nc))
        n = len(pattern)
        handle = sparse(pattern)
        try:
            for initial in initials:
                values = pack(handle,matrices,initial)
                if inplace:
                    require(lu_factor(handle,nc,ptr(values),algorithm))
                else:
                    original = values.copy()
                    require(lu_factor_copy(handle,nc,ptr(original),ptr(values),algorithm))
                    near(original,pack(handle,matrices),rtol=0)
                lu = unpack(handle,values,n)
                for c in range(nc):
                    lower,upper = np.tril(lu[c],-1)+np.eye(n),np.triu(lu[c])
                    # Same L*U quantities and absolute bounds as the MICM cases.
                    near(lower@upper,matrices[c],rtol=0,atol=1e-12 if '.Dense' not in id_ else 1e-10)
                if linear:
                    x = rhs.copy()
                    lu_solve(handle,nc,ptr(values),ptr(x))
                    for c in range(nc):
                        got = matrices[c]@x[:,c]
                        if '.Diagonal' in id_:
                            near(got,rhs[:,c],rtol=4*np.finfo(float).eps)
                        else:
                            near(got,rhs[:,c],rtol=0,atol=2e-6 if 'InitialValue' in id_ else 1e-5)
        finally:
            sparse_free(handle)
