"""Finite-difference perturbations, comparisons, and sparsity diagnostics."""
from functools import partial
from chemffi import DP
from sparseffi import *

FORCING = C.CFUNCTYPE(None,Z,Z,DP,DP,P)


class Comparison(C.Structure):
    _fields_ = [('passed',I),('max_absolute_error',D),('max_relative_error',D),
                ('worst_block',Z),('worst_row',Z),('worst_column',Z),
                ('worst_analytical',D),('worst_finite_difference',D)]


bind('chem_finite_difference',I,Z,Z,FORCING,P,DP,D,DP)
bind('chem_compare_jacobian',I,P,DP,D,D,C.POINTER(Comparison))
bind('chem_check_jacobian_sparsity',I,P,DP,D,C.POINTER(Comparison))


def finite_difference(rows,forcing,perturbation=1e-8):
    nc,n = len(rows),len(rows[0])
    values = (D*(nc*n))(*[rows[c][i] for i in range(n) for c in range(nc)])
    out = (D*(nc*n*n))()
    errors = []
    def wrapped(cells,species,y,f,ctx):
        try:forcing(cells,species,y,f)
        except Exception as error:errors.append(error)
    callback = FORCING(wrapped)
    assert lib.chem_finite_difference(nc,n,callback,None,values,perturbation,out) == 0
    if errors:raise errors[0]
    return out


def simple(nc,n,y,f):
    for c in range(nc):
        x,v = y[c],y[nc+c]
        f[c] += x*v
        f[nc+c] += x*x-v


def difference(mode=0):
    rows = [[3.,5.]] if mode==0 else ([[1.,2.],[.5,4.],[10.,.1]] if mode==1 else [[0.,1.]])
    fd = finite_difference(rows,simple)
    for c,(x,y) in enumerate(rows):
        for j,value in enumerate([y,x,2*x,-1.]):assert abs(fd[j*len(rows)+c]-value) <= 1e-5


def comparison(wrong=False):
    fd = finite_difference([[3.,5.]],simple)
    with Sparse(2,entries=[(0,0),(0,1),(1,0),(1,1)]) as a:
        for (i,j),value in zip([(0,0),(0,1),(1,0),(1,1)],[-999. if wrong else -5.,-3.,-6.,1.]):a[0,i,j] = value
        out = Comparison()
        assert lib.chem_compare_jacobian(a.handle,fd,1e-7,1e-7,C.byref(out)) == 0
        assert bool(out.passed) == (not wrong)
        if wrong:assert (out.worst_row,out.worst_column) == (0,0)
        else:assert out.max_absolute_error < 1e-5


def completeness(missing=False):
    fd = finite_difference([[3.,5.]],simple)
    entries = [(0,0),(0,1),(1,1)] if missing else [(0,0),(0,1),(1,0),(1,1)]
    with Sparse(2,entries=entries) as a:
        out = Comparison()
        assert lib.chem_check_jacobian_sparsity(a.handle,fd,1e-6,C.byref(out)) == 0
        assert bool(out.passed) == (not missing)
        if missing:assert (out.worst_row,out.worst_column) == (1,0)


CASES = {
    'JacobianVerification.FiniteDifferenceMatchesAnalytical':difference,
    'JacobianVerification.MultiBlockFiniteDifference':partial(difference,1),
    'JacobianVerification.ComparePassesForCorrectJacobian':comparison,
    'JacobianVerification.CompareFailsForWrongJacobian':partial(comparison,True),
    'JacobianVerification.SparsityCompletenessPassesWhenComplete':completeness,
    'JacobianVerification.SparsityCompletenessFailsWhenMissingEntry':partial(completeness,True),
    'JacobianVerification.NearZeroVariableHandled':partial(difference,2),
}
