"""Keep callback ownership explicit when exercising the C extension interface."""
from api import *
from harness import equal


def custom(m,label,initial=0):
    index = Z()
    equal(lib.chem_model_custom(m.p,label.encode(),initial,C.byref(index)),0)
    return index.value


def external(m,forcing=None,jacobian=None,update=None,setup=None):
    functions = [kind(fn) if fn else None for kind,fn in
                 [(ExternalKernel,forcing),(ExternalKernel,jacobian),(ExternalUpdate,update),(ExternalSetup,setup)]]
    m.callbacks.extend(fn for fn in functions if fn)
    item = External(*(C.cast(fn,P) if fn else None for fn in functions),None)
    equal(lib.chem_model_external(m.p,item),0)
    return functions


def transfers(m,reactions,update=None):
    """An external process model; cache all index lookups during finalization."""
    cache = []
    def setup(mp,context):
        n = mp.contents.nvariables
        for r,p,k in reactions:
            ri,pi = m.var(r),m.var(p)
            cache.append((ri,pi,k,mp.contents.index[ri*n+ri],mp.contents.index[pi*n+ri]))
    def forcing(mp,sp,y,f,context):
        s = sp.contents; stride = s.stride
        for r,p,k,_,_ in cache:
            for c in range(s.cells):
                v = (s.custom[k*stride+c] if isinstance(k,int) else k)*y[r*stride+c]
                f[r*stride+c] -= v; f[p*stride+c] += v
    def jacobian(mp,sp,y,j,context):
        s = sp.contents; stride = s.stride
        for r,p,k,rr,pr in cache:
            for c in range(s.cells):
                v = s.custom[k*stride+c] if isinstance(k,int) else k
                j[rr*stride+c] += v; j[pr*stride+c] -= v
    return external(m,forcing,jacobian,update,setup)
