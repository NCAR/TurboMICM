"""Ports of test_matrix_policy.hpp, shared by row-major and grouped arrays.

Callbacks use the C views and dispatch; Python only supplies the test operations
and independently checks the resulting values, as a caller of the C API would.
"""
from contextlib import ExitStack
import math
from arrays import *
from harness import case

CHECKS = {}
def check(*names):
    def register(fun):
        for name in names: CHECKS[name] = fun
        return fun
    return register


@check('SmallMatrix','SmallConstMatrix')
def small(make):
    a = make(3,5)
    a[1,3],a[0,0],a[2,4] = 64.7,41.2,102.3
    near(a.values(),[[41.2,0,0,0,0],[0,0,0,64.7,0],[0,0,0,0,102.3]],rtol=0)
    g = a.group
    assert len(a.raw) == ((3+g-1)//g*g if g else 3)*5
    for r,c,v in ((1,3,64.7),(0,0,41.2),(2,4,102.3)):
        index = (r//g)*g*5+c*g+r%g if g else r*5+c
        assert a.raw[index] == v
    if g:
        assert array_groups(a.handle) == (3+g-1)//g
        assert array_group_size(a.handle) == g


@check('InitializeMatrix','InitializeConstMatrix')
def initialize(make):
    a = make(2,3,initial=12.4)
    near(a.values(),np.full((2,3),12.4),rtol=0)


@check('LoopOverMatrix','LoopOverConstMatrix')
def loop(make):
    a = make(3,4,integer=True)
    expected = [[r*100+c for c in range(4)] for r in range(3)]
    a.initialize(expected)
    near(a.values(),expected,rtol=0)


@check('IterateOverMatrix')
def row_iteration(make):
    a = make(3,4,integer=True)
    require(array_set_row(a.handle,1,doubles([42,43,44,45]),4))
    near(a.values(),[[0,0,0,0],[42,43,44,45],[0,0,0,0]],rtol=0)


@check('Strides')
def strides(make):
    a = make(3,4)
    assert array_row_stride(a.handle) == (1 if a.group else 4)
    assert array_column_stride(a.handle) == (a.group or 1)


@check('ConversionToVector','ConstConversionToVector')
def to_vector(make):
    a = make(2,3)
    a.initialize([[0,0,0],[13.2,31.2,314.2]])
    near(a.values()[1],[13.2,31.2,314.2],rtol=0)


@check('ConversionFromVector')
def from_vector(make):
    a = make(2,3)
    expected = [[412.3,32.4,41.3],[5.33,-.3,31.2]]
    a.initialize(expected)
    near(a.values(),expected,rtol=0)
    assert array_set_row(a.handle,1,doubles([1,2]),2) == -4
    near(a.values(),expected,rtol=0)


@check('AssignmentFromVector')
def assign_vector(make):
    a = make(4,3)
    require(array_set_row(a.handle,2,doubles([12.3,15.1,24.3]),3))
    near(a.values(),[[0,0,0],[0,0,0],[12.3,15.1,24.3],[0,0,0]],rtol=0)
    require(array_set_row(a.handle,2,doubles([14.3,52.3,65.7,8]),4))
    near(a.values()[2],[14.3,52.3,65.7],rtol=0)
    assert array_set_row(a.handle,2,doubles([1,2]),2) == -4


@check('Axpy')
def axpy(make):
    a,b = make(4,3),make(4,3)
    y,x = np.indices((4,3)),np.indices((4,3))
    y = y[0]*10.3+y[1]*100.5
    x = x[0]*1.7+x[1]*10.2
    a.initialize(y); b.initialize(x)
    require(array_axpy(a.handle,1.39,b.handle))
    near(a.values(),y+1.39*x,rtol=0,atol=1e-5)
    assert abs(a.values().sum()-(y+1.39*x).sum()) <= 1e-5


@check('ForEach')
def for_each(make):
    a,b,c = (make(4,3) for _ in range(3))
    r,j = np.indices((4,3))
    y,x,z = r*10.3+j*100.5,r*1.7+j*10.2,r*19.5+j*32.2
    a.initialize(y); b.initialize(x); c.initialize(z)
    @Element
    def add(v,n,_): v[0] += v[1]
    require(array_map(2,(P*2)(a.handle,b.handle),(C.c_ubyte*2)(0,1),add,None))
    near(a.values(),y+x,rtol=0,atol=1e-5)
    a.initialize(y)
    @Element
    def subtract(v,n,_): v[0] = v[0]+v[1]-v[2]
    require(array_map(3,(P*3)(a.handle,b.handle,c.handle),(C.c_ubyte*3)(0,1,1),subtract,None))
    near(a.values(),y+x-z,rtol=0,atol=1e-5)


@check('SetScaler')
def scalar(make):
    a = make(2,3)
    require(array_fill(a.handle,2))
    near(a.raw,np.full(len(a.raw),2),rtol=0)


@check('Max','Min')
def bounds(make):
    a = make(2,3)
    require(array_max(a.handle,2))
    near(a.raw,np.full(len(a.raw),2),rtol=0)
    require(array_fill(a.handle,0)); require(array_min(a.handle,2))
    near(a.raw,np.zeros(len(a.raw)),rtol=0)
    require(array_fill(a.handle,1)); a[1,1] = 3
    require(array_max(a.handle,2))
    near(a.values(),[[2,2,2],[2,3,2]],rtol=0)
    require(array_fill(a.handle,1)); a[1,1] = 3
    require(array_min(a.handle,2))
    near(a.values(),[[1,1,1],[1,2,1]],rtol=0)


@check('Print')
def printing(make):
    a = make(2,3); a[1,1] = 3
    assert a.text() == '0,0,0\n0,3,0\n'


def sum_op(groups,factor=4,vector=False):
    g = groups[0]
    tmp = temporary(g)
    assign(tmp,[column(g,j) for j in range(3)],lambda a,b,c:a+b+c)
    if vector:
        assign(column(g,2),[column(groups[1],0),tmp],lambda v,t:v*t)
    else:
        assign(column(g,2),[tmp],lambda t:factor*t)


@check('ArrayFunction','VectorInMatrixFunction')
def array_function(make):
    a = make(5,3)
    values = np.array([[i-2+10*j for j in range(3)] for i in range(5)])
    a.initialize(values)
    with Function([a],sum_op) as f:
        f(a)
        expected = values.copy(); expected[:,2] = 4*values.sum(axis=1)
        near(a.values(),expected,rtol=0)
        b = make(3,3,initial=-1); f(b)
        near(b.values(),[[-1,-1,-12]]*3,rtol=0)
    a.initialize(values)
    v = make(5,vector=True).initialize([100,200,300,400,500])
    with Function([a,v],lambda g:sum_op(g,vector=True),[0,1]) as f:
        f(a,v)
    expected = values.copy(); expected[:,2] = np.arange(1,6)*100*values.sum(axis=1)
    near(a.values(),expected,rtol=0)


@check('MultiMatrixArrayFunction')
def multi_array(make):
    a,b = make(3,2),make(3,3)
    av,bv = np.array([[i,i+10] for i in range(3)]),np.array([[2*i,2*i+20,4*i] for i in range(3)])
    a.initialize(av); b.initialize(bv)
    def op(gs):
        tmp = temporary(gs[0])
        assign(tmp,[column(gs[0],0),column(gs[1],2)],lambda x,y:x+y)
        assign(column(gs[0],1),[tmp],lambda x:x)
    with Function([a,b],op,[0,1]) as f: f(a,b)
    near(a.values(),[[0,0],[1,5],[2,10]],rtol=0)
    near(b.values(),bv,rtol=0)


@check('MultipleTemporaries')
def temporaries(make):
    a = make(4,5)
    for i in range(4): a[i,0],a[i,1] = i+1,(i+1)*10
    def op(gs):
        g=gs[0]; t,u=temporary(g),temporary(g)
        assign(t,[column(g,0),column(g,1)],lambda x,y:x*y)
        assign(u,[column(g,0),column(g,1)],lambda x,y:x+y)
        for j,expr in ((2,lambda x,y:x+y),(3,lambda x,y:x-y),(4,lambda x,y:x*y)):
            assign(column(g,j),[t,u],expr)
    with Function([a],op) as f: f(a)
    near(a.values(),[[i,10*i,10*i*i+11*i,10*i*i-11*i,110*i**3] for i in range(1,5)],rtol=0)


@check('ColumnViewReuse')
def reuse_columns(make):
    a = make(3,4)
    for i in range(3): a[i,0] = i+1
    def op(gs):
        v=[column(gs[0],i) for i in range(4)]
        assign(v[1],v[:1],lambda x:2*x)
        assign(v[2],v[:2],lambda x,y:x+y)
        assign(v[3],[v[2],v[1]],lambda x,y:x*y)
    with Function([a],op) as f: f(a)
    near(a.values(),[[1,2,3,6],[2,4,6,24],[3,6,9,54]],rtol=0)


@check('FunctionReusability')
def reusability(make):
    a=make(2,3).initialize([[0,1,2],[1,2,3]])
    with Function([a],lambda g:sum_op(g,2)) as f:
        f(a); near(a.values(),[[0,1,6],[1,2,12]],rtol=0)
        b=make(2,3,initial=5); f(b); near(b.values(),[[5,5,30]]*2,rtol=0)
        c=make(2,3).initialize([[0,0,0],[10,0,0]]); f(c)
        near(c.values(),[[0,0,0],[10,0,20]],rtol=0)


@check('ConstMatrixFunction')
def const_matrix(make):
    values=np.array([[r*10+c for c in range(4)] for r in range(3)])
    a=make(3,4).initialize(values)
    total=D()
    def op(gs):
        v=column(gs[0],2)
        assert view_set(C.byref(v),0,99) == -5
        reduce_view(v,total,lambda x,acc:acc+x,True)
        t=temporary(gs[0]); assign(t,[column(gs[0],0),column(gs[0],1)],lambda x,y:x+y)
    with Function([a],op,[1]) as f:f(a)
    assert total.value == 36
    near(a.values(),values,rtol=0)


@check('EmptyMatrixFunction','EmptyVectorEmptyMatrix')
def empty(make):
    a,v=make(0,3),make(0,vector=True)
    def never(gs): raise AssertionError('empty dispatch invoked callback')
    with Function([a,v],never) as f:f(a,v)
    b=make(3,0)
    with Function([b],lambda gs:None) as f:f(b)


@check('MultiMatrixDifferentRowsFromCreation','MatrixVectorDifferentRowsFromCreation')
def flexible(make):
    a,b,v=make(3,2),make(3,3),make(3,vector=True)
    def op(gs,factor,column_index):
        t=temporary(gs[0]); assign(t,[column(gs[0],0),column(gs[1],column_index)],lambda x,y:x+y)
        assign(column(gs[0],1),[t],lambda x:factor*x)
    with Function([a,b],lambda g:op(g,2,2),[0,1]) as f:
        c,d=make(5,2),make(5,3)
        for i in range(5): c[i,0],d[i,2]=i+1,i*10
        f(c,d); near(c.values()[:,1],[2,24,46,68,90],rtol=0)
        e,h=make(2,2),make(2,3)
        e[0,0],e[1,0],h[0,2],h[1,2]=10,20,5,15
        f(e,h); near(e.values()[:,1],[30,70],rtol=0)
    with Function([a,v],lambda g:op(g,3,0),[0,1]) as f:
        c,d=make(5,2),make(5,vector=True)
        for i in range(5): c[i,0],d[i,0]=i+1,i*10
        f(c,d);near(c.values()[:,1],[3,36,69,102,135],rtol=0)


@check('MismatchedRowDimensions','MismatchedRowsAtInvocation','MultipleMatricesMismatchedRowsAtInvocation',
       'VectorTooSmall','VectorTooLarge','EmptyVectorNonEmptyMatrix','NonEmptyVectorEmptyMatrix',
       'MultipleVectorsDifferentSizes','MultipleMatricesDifferentRowsVector','VectorSizeMatchesOneMatrixOnly',
       'FunctionInvocationWithWrongSizedVector')
def bad_rows(make):
    # Each source error scenario is exercised, including mismatches at creation
    # being accepted and detected only when the function is invoked.
    shapes=[[(3,3,False),(4,3,False)],[(5,2,False),(3,1,True)],[(3,2,False),(5,1,True)],
            [(5,3,False),(0,1,True)],[(0,3,False),(3,1,True)],
            [(5,3,False),(5,1,True),(3,1,True)],[(4,2,False),(5,3,False),(4,1,True)],
            [(5,2,False),(5,3,False),(4,1,True)],[(3,2,False),(3,3,False),(4,3,False)]]
    for shape in shapes:
        arrays=[make(r,c,vector=v) for r,c,v in shape]
        with Function(arrays,lambda gs:None) as f:
            assert f.call(*arrays)==-4
            good=[make(5,c,vector=v) for _,c,v in shape]
            f(*good)


@check('MismatchedColumnDimensions','WrongMatrixDimensions','WrongColumnCountAtInvocation')
def bad_columns(make):
    a=make(3,4,initial=1)
    with Function([a],lambda gs:None) as f:
        f(a);f(make(7,4))
        for cols in (2,3,5): assert f.call(make(3,cols))==-4
    with Function([a],lambda gs:group_view(C.byref(gs[0]),5,C.byref(View()))) as f:
        assert f.call(a)==-4


@check('MultipleVectorsSameSize')
def two_vectors(make):
    a=make(5,3);v=make(5,vector=True).initialize(range(1,6));w=make(5,vector=True).initialize(np.arange(1,6)*10)
    with Function([a,v,w],lambda g:assign(column(g[0],0),[column(g[1],0),column(g[2],0)],lambda x,y:x+y),[0,1,1]) as f:
        f(a,v,w)
    near(a.values(),[[11*i,0,0] for i in range(1,6)],rtol=0)


@check('MultipleMatricesOneVector')
def two_matrices_vector(make):
    a,b,v=make(4,2),make(4,3),make(4,vector=True)
    for i in range(4): a[i,0],b[i,0],v[i,0]=i+1,i*10,i*2
    def op(g):
        for m in g[:2]:assign(column(m,1),[column(m,0),column(g[2],0)],lambda x,y:x+y)
    with Function([a,b,v],op,[0,0,1]) as f:f(a,b,v)
    near(a.values(),[[1,1],[2,4],[3,7],[4,10]],rtol=0)
    near(b.values(),[[0,0,0],[10,12,0],[20,24,0],[30,36,0]],rtol=0)


@check('ConstVector','IntegerVector','FunctionWithConstSignature')
def const_vector(make):
    for integer in (False,True):
        a=make(3,2,integer=integer);v=make(3,vector=True,integer=integer).initialize([10,20,30])
        with Function([a,v],lambda g:assign(column(g[0],0),[column(g[1],0)],lambda x:2*x),[0,1]) as f:
            f(a,v)
        near(a.values(),[[20,0],[40,0],[60,0]],rtol=0)
        near(v.values().ravel(),[10,20,30],rtol=0)


@check('MutableVector')
def mutable_vector(make):
    a=make(3,2).initialize([[1,0],[2,0],[3,0]]);v=make(3,vector=True).initialize([5,10,15])
    with Function([a,v],lambda g:assign(column(g[1],0),[column(g[0],0)],lambda x:3*x),[1,0]) as f:f(a,v)
    near(v.values().ravel(),[3,6,9],rtol=0)
    near(a.values(),[[1,0],[2,0],[3,0]],rtol=0)


@check('FunctionReusabilityWithVectors')
def reuse_vector(make):
    a=make(3,2);v=make(3,vector=True).initialize([1,2,3]);w=make(3,vector=True).initialize([5,6,7])
    with Function([a,v],lambda g:assign(column(g[0],0),[column(g[1],0)],lambda x:10*x),[0,1]) as f:
        f(a,v);near(a.values(),[[10,0],[20,0],[30,0]],rtol=0)
        f(a,w);near(a.values(),[[50,0],[60,0],[70,0]],rtol=0)


@check('MixedVectorColumnViewRowVariable')
def mixed_views(make):
    a=make(4,3).initialize([[i,10*i,0] for i in range(1,5)])
    v=make(4,vector=True).initialize([100,200,300,400])
    def op(g):
        t=temporary(g[0]);assign(t,[column(g[0],0),column(g[0],1)],lambda x,y:x+y)
        assign(column(g[0],2),[t,column(g[1],0)],lambda x,y:x+y)
    with Function([a,v],op,[0,1]) as f:f(a,v)
    near(a.values(),[[i,10*i,111*i] for i in range(1,5)],rtol=0)


@check('TestFill','TestCopy')
def fill_copy(make):
    a=make(3,2);v=make(3,vector=True)
    def fill(g):
        require(view_fill(C.byref(column(g[0],1)),3.2,0))
        require(view_fill(C.byref(column(g[1],0)),3.2,1))
        t=temporary(g[0]);require(view_fill(C.byref(t),9.9,0))
        require(view_copy(C.byref(column(g[0],0)),C.byref(t),1))
    with Function([a,v],fill) as f:f(a,v)
    near(a.values(),[[9.9,3.2]]*3,rtol=0);near(v.values().ravel(),[3.2]*3,rtol=0)
    v.initialize([3.2,4.2,1.3]);w=make(3,vector=True)
    def copy(g):
        c0,c1=column(g[0],0),column(g[0],1)
        require(view_copy(C.byref(c1),C.byref(column(g[1],0)),1))
        require(view_copy(C.byref(column(g[2],0)),C.byref(c1),1))
        require(view_copy(C.byref(c0),C.byref(c1),1))
        require(view_fill(C.byref(c0),0,0))
        require(view_copy(C.byref(c0),C.byref(column(g[0],1,1)),1))
        t=temporary(g[0]);require(view_copy(C.byref(t),C.byref(c1),1))
        require(view_fill(C.byref(c1),0,0));assign(c1,[t],lambda x:x,1)
    with Function([a,v,w],copy,[0,1,0]) as f:f(a,v,w)
    near(a.values(),[[3.2,3.2],[4.2,4.2],[1.3,1.3]],rtol=0)
    near(w.values().ravel(),[3.2,4.2,1.3],rtol=0)


@check('ReduceSum','ReduceMax','ReduceLOr','ReduceLAnd','ReduceStrict')
def reductions(make):
    scenarios=[([[1,2],[3,4],[5,6]],0,lambda x,a:a+x*x,False,91),
               ([[1,8],[5,2],[3,4]],-1e300,lambda x,a:max(x,a),False,8),
               ([[1,1]]*3,0,lambda x,a:a or x>10,False,0),
               ([[1,1],[1,1],[1,42]],0,lambda x,a:a or x>10,False,1),
               ([[1,1]]*3,1,lambda x,a:a and math.isfinite(x),False,1),
               ([[1,1],[float('nan'),1],[1,1]],1,lambda x,a:a and math.isfinite(x),False,0),
               ([[0],[0],[0]],0,lambda x,a:a+1,True,3)]
    for values,start,expression,strict,expected in scenarios:
        a=make(3,len(values[0])).initialize(values);acc=D(start)
        def op(gs):
            for j in range(a.elements):reduce_view(column(gs[0],j),acc,expression,strict)
        with Function([a],op,[1]) as f:f(a)
        assert acc.value==expected,(acc.value,expected)


@check('VectorCapture')
def capture(make):
    a=make(3,3);v=make(3,vector=True,integer=True).initialize([2,0,1]);w=make(3,vector=True,integer=True).initialize([10,20,30])
    def op(gs):
        c0,c1=column(gs[0],0),column(gs[0],1)
        for i in range(v.rows):assign(c0,[c0],lambda x:x+w[int(v[i,0]),0])
        assign(c1,[c1],lambda x:x+w[int(v[1,0]),0])
    with Function([a],op) as f:f(a)
    near(a.values(),[[60,10,0]]*3,rtol=0)


def run_dense(id_):
    name=id_.split('.')[-1].replace('VectorMatrix','Matrix')
    widths=(1,2,3,4) if id_.startswith('vector_matrix:') else (0,)
    for width in widths:
        with ExitStack() as stack:
            def make(r,c=1,**kw):
                return stack.enter_context(Array(r,c,group=width,**kw))
            CHECKS[name](make)


for name in CHECKS:
    for prefix in ('matrix:Matrix.','vector_matrix:VectorMatrix.'):
        n=name
        if prefix.startswith('vector_'):
            if name=='IterateOverMatrix':continue
            if name in ('SmallMatrix','SmallConstMatrix','InitializeMatrix','InitializeConstMatrix','LoopOverMatrix','LoopOverConstMatrix'):
                n=name.replace('Matrix','VectorMatrix')
        case(prefix+n)(run_dense)
