from arrays import *
from harness import case, equal, near
from sparse_cases import SparsePool


class Comparison(C.Structure):
    _fields_ = [('passed',I), ('max_abs_error',D), ('max_rel_error',D)] + [
        (name,Z) for name in ('worst_block','worst_row','worst_col')] + [('worst_analytical',D), ('worst_fd',D)]


ForcingFn = C.CFUNCTYPE(None,AP,AP,P)
bind('chem_finite_difference',I,AP,Z,D,ForcingFn,P,C.POINTER(AP))
bind('chem_compare_jacobian',Comparison,AP,AP,Z,D,D)
bind('chem_check_sparsity',Comparison,AP,AP,Z,D)


@case(*['jacobian_verification:JacobianVerification.'+name for name in (
    'FiniteDifferenceMatchesAnalytical','MultiBlockFiniteDifference','ComparePassesForCorrectJacobian',
    'CompareFailsForWrongJacobian','SparsityCompletenessPassesWhenComplete','SparsityCompletenessFailsWhenMissingEntry','NearZeroVariableHandled')])
def verification(identifier):
    name = identifier.split('.')[-1]
    inputs = [(1,2),(.5,4),(10,.1)] if name == 'MultiBlockFiniteDifference' else [(0,1)] if name == 'NearZeroVariableHandled' else [(3,5)]
    with SparsePool(1,False) as p:
        v = p.matrix(len(inputs),2)
        for b,(x,y) in enumerate(inputs): v[b,0] = x; v[b,1] = y
        def forcing(vp,fp,context):
            v,f = Array(pointer=vp),Array(pointer=fp)
            for b in range(v.data.rows):
                x,y = v[b,0],v[b,1]; f[b,0] += x*y; f[b,1] += x*x-y
        callback = ForcingFn(forcing); pointer = AP()
        equal(lib.chem_finite_difference(v.p,2,1e-8,callback,None,C.byref(pointer)),0)
        fd = Array(pointer=pointer); fd.owned = True; p.resources.append(fd)
        for b,(x,y) in enumerate(inputs):
            for j,expected in enumerate([y,x,2*x,-1]): near(fd[b,j],expected,1e-5)
        a = p.sparse(len(inputs),2,[(0,0),(0,1),(1,0),(1,1)])
        for b,(x,y) in enumerate(inputs):
            for pair,expected in zip([(0,0),(0,1),(1,0),(1,1)],[-y,-x,-2*x,1]): a[(b,)+pair] = expected
        result = lib.chem_compare_jacobian(a.p,fd.p,2,1e-7,1e-7); equal(result.passed,1); assert result.max_abs_error < 1e-5
        if 'WrongJacobian' in name:
            a[0,0,0] = -999; result = lib.chem_compare_jacobian(a.p,fd.p,2,1e-7,1e-7)
            equal(result.passed,0); equal(result.worst_row,0); equal(result.worst_col,0)
        equal(lib.chem_check_sparsity(a.p,fd.p,2,1e-6).passed,1)
        if 'MissingEntry' in name:
            incomplete = p.sparse(len(inputs),2,[(0,0),(0,1),(1,1)])
            result = lib.chem_check_sparsity(incomplete.p,fd.p,2,1e-6)
            equal(result.passed,0); equal(result.worst_row,1); equal(result.worst_col,0)
