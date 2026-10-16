from ffi import *
from harness import group


@group('version:Version.','FullVersion VersionMajor VersionMinor VersionPatch VersionTweak')
def version_api(_):
    assert version().decode()=='.'.join(str(version_component(i)) for i in range(4))
    assert version_component(0)>=1 and version_component(4)==-1


@group('jacobian_verification:JacobianVerification.',
       'FiniteDifferenceMatchesAnalytical MultiBlockFiniteDifference ComparePassesForCorrectJacobian CompareFailsForWrongJacobian SparsityCompletenessPassesWhenComplete SparsityCompletenessFailsWhenMissingEntry NearZeroVariableHandled')
def jacobian_verification(_):
    @ForcingFunction
    def function(n,c,y,f,ctx):
        for i in range(c):
            x,v=y[i],y[c+i]
            f[i]+=x*v;f[c+i]+=x*x-v
    for values in ([[3,5]],[[1,2],[.5,4],[10,.1]],[[0,1]]):
        y=np.ascontiguousarray(np.array(values,dtype=float).T);n,c=y.shape
        fd=np.zeros((2,2,c));require(finite_difference(n,c,ptr(y),1e-8,function,None,ptr(fd)))
        analytical=np.array([[y[1],y[0]],[2*y[0],np.full(c,-1.)]])
        near(fd,analytical,rtol=0,atol=1e-5)
        pat=np.ones((2,2),np.uint8);negative=-analytical
        result=compare_jacobian(2,c,pat.ctypes.data_as(BP),ptr(negative),ptr(fd),1e-7,1e-7)
        assert result.passed and result.max_absolute_error<1e-5
        negative[0,0]=-999
        result=compare_jacobian(2,c,pat.ctypes.data_as(BP),ptr(negative),ptr(fd),1e-7,1e-7)
        assert not result.passed and result.worst_row==0 and result.worst_column==0
        assert check_sparsity(2,c,pat.ctypes.data_as(BP),ptr(fd),1e-6).passed
        if y[0,0]:
            pat[1,0]=0;result=check_sparsity(2,c,pat.ctypes.data_as(BP),ptr(fd),1e-6)
            assert not result.passed and result.worst_row==1 and result.worst_column==0
