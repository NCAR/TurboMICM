"""Raw-to-group view equivalence and C equivalents of the category traits."""
from contextlib import ExitStack
from arrays import *
from harness import case,group


@group('view_category_sfinae:ViewCategorySFINAE.',
       'VectorTypes MatrixViews SparseMatrixViews MatrixTypesExcluded RuntimeVectorUse')
@case('grouped_views:GroupedView.CategoryTagsAreDistinct')
def categories(_):
    with ExitStack() as stack:
        arrays=[stack.enter_context(Array(5,2,group=4)),
                stack.enter_context(Array(5,2,pattern=[(0,1)],group=4)),
                stack.enter_context(Array(5,vector=True))]
        views=[column(a.view_group(0),0) for a in arrays]
        t=temporary(arrays[0].view_group(0))
        assert [v.category for v in views]+[t.category]==[0,1,2,3]
        # A full array cannot be confused with a vector when dispatching a
        # function: kind checks replace C++ overload constraints.
        with Function([arrays[2]],lambda g:None) as f:assert f.call(arrays[0])<0
        arrays[2].initialize([1,2,3,4,5]);acc=D()
        with Function([arrays[2]],lambda g:reduce_view(column(g[0],0),acc,lambda x,a:a+x,True),[1]) as f:f(arrays[2])
        assert acc.value==15


@case('grouped_views:GroupedView.MatrixColumnViewRoundTrip')
def dense_roundtrip(_):
    with Array(3,4) as a:
        a.initialize([[10*r+c for c in range(4)] for r in range(3)])
        g=a.view_group(1);v,w=column(g,2),column(g,3,1)
        assert C.cast(v.data,DP)[0]==12 and C.cast(w.data,DP)[0]==13
        assign(v,[w],lambda x:2*x)
        near(a.values(),[[0,1,2,3],[10,11,26,13],[20,21,22,23]],rtol=0)


@case('grouped_views:GroupedView.VectorMatrixColumnViewContiguousBlock')
def dense_contiguous(_):
    with Array(8,3,group=4) as a:
        original=np.array([[100+10*r+c for c in range(3)] for r in range(8)])
        a.initialize(original);g=a.view_group(1)
        v,w=column(g,1),column(g,2,1)
        assert v.stride==w.stride==1
        near(np.ctypeslib.as_array(C.cast(v.data,DP),shape=(4,)),original[4:,1],rtol=0)
        near(np.ctypeslib.as_array(C.cast(w.data,DP),shape=(4,)),original[4:,2],rtol=0)
        assign(v,[w],lambda x:2*x)
        original[4:,1]=2*original[4:,2];near(a.values(),original,rtol=0)


@group('grouped_views:GroupedView.','DenseGroupedMatchesRawL1 DenseGroupedMatchesRawL2 DenseGroupedMatchesRawL4 DenseGroupedMatchesRawL8')
def dense_equal(id_):
    width=int(id_[-1]);rows=3*width
    with Array(rows,4,group=width) as a,Array(rows,4,group=width) as b:
        original=np.array([[r*10+c for c in range(4)] for r in range(rows)])
        a.initialize(original);b.initialize(original)
        with Function([a],lambda g:assign(column(g[0],3),[column(g[0],0)],lambda x:3*x+1)) as f:f(a)
        for r in range(rows):b[r,3]=b[r,0]*3+1
        near(a.values(),b.values(),rtol=0)


@case('grouped_views:GroupedView.PlainMatrixGroupedRoundTrip')
def plain(_):
    with Array(4,3) as a:
        a.initialize([[r+5*c for c in range(3)] for r in range(4)])
        with Function([a],lambda g:assign(column(g[0],2),[column(g[0],0),column(g[0],1)],lambda x,y:x-y)) as f:f(a)
        near(a.values(),[[r,r+5,-5] for r in range(4)],rtol=0)


@case('grouped_views:GroupedView.SparseStandardBlockViewRoundTrip')
def sparse_roundtrip(_):
    with Array(2,3,pattern=[(0,1),(1,2),(2,2)]) as a:
        a.initialize([[1,2,3],[10,20,30]])
        g=a.view_group(1);v,w=column(g,0),column(g,1,1)
        assert C.cast(v.data,DP)[0]==10 and C.cast(w.data,DP)[0]==20
        assign(v,[w],lambda x:x+5)
        near(a.values(),[[1,2,3],[25,20,30]],rtol=0)


@case('grouped_views:GroupedView.SparseVectorBlockViewContiguousBlock')
def sparse_contiguous(_):
    with Array(8,3,pattern=[(0,1),(1,2)],group=4) as a:
        original=np.array([[100+i,200+i] for i in range(8)])
        a.initialize(original);g=a.view_group(1);v,w=column(g,0),column(g,1)
        for view,j in ((v,0),(w,1)):
            assert view.stride==1
            near(np.ctypeslib.as_array(C.cast(view.data,DP),shape=(4,)),original[4:,j],rtol=0)
        assign(v,[v,w],lambda x,y:x+y)
        original[4:,0]+=original[4:,1];near(a.values(),original,rtol=0)


@group('grouped_views:GroupedView.',
       'MixedSparseAndDenseGroupedL1 MixedSparseAndDenseGroupedL2 MixedSparseAndDenseGroupedL4 MixedSparseAndDenseGroupedL8')
def mixed(id_):
    width=int(id_[-1]);rows=2*width
    with Array(rows,3,pattern=[(0,1)],group=width) as s,Array(rows,3,group=width) as d:
        s.initialize([[1+i] for i in range(rows)])
        original=np.array([[10+i+c for c in range(3)] for i in range(rows)]);d.initialize(original)
        with Function([s,d],lambda g:assign(column(g[1],0),[column(g[1],1),column(g[0],0)],lambda x,y:x+y),[1,0]) as f:f(s,d)
        original[:,0]=original[:,1]+1+np.arange(rows);near(d.values(),original,rtol=0)
        near(s.values().ravel(),np.arange(rows)+1,rtol=0)


@case('grouped_views:GroupedView.SparseVectorRowColOverloadStillFunctional')
def row_column(_):
    with Array(8,3,pattern=[(0,1)],group=4) as a:
        a.initialize(np.arange(1,9));g=a.view_group(1,1);v=View()
        require(group_sparse_view(C.byref(g),0,1,C.byref(v)))
        acc=D();reduce_view(v,acc,lambda x,a:a+x)
        assert acc.value==26


@case('grouped_views:GroupedView.ConstGroupViewProducesGroupedViews')
def const_group(_):
    with Array(4,3,pattern=[(0,1)],group=4) as s,Array(4,2,group=4) as d:
        s.initialize([1,2,3,4]);d.initialize([[0,100+i] for i in range(4)])
        for a,e,expected in ((s,0,[1,2,3,4]),(d,1,[100,101,102,103])):
            v=column(a.view_group(0,1),e)
            assert v.read_only==1
            assert view_set(C.byref(v),0,10)==-5
            near(np.ctypeslib.as_array(C.cast(v.data,DP),shape=(4,)),expected,rtol=0)


@group('grouped_views:GroupedView.',
       'FillAndCopySparseStandardCSR FillAndCopySparseStandardCSC FillAndCopySparseVectorCSR_L4 FillAndCopySparseVectorCSC_L4 FillAndCopySparseVectorCSR_L8')
def fillcopy(id_):
    width=int(id_[-1]) if '_L' in id_ else 1
    order=int('CSC' in id_)
    with Array(width+1,3,pattern=[(0,1),(2,0)],order=order,group=width if '_L' in id_ else 0) as a:
        e,h=a.element(0,1),a.element(2,0)
        for i in range(width+1):a[i,e],a[i,h]=i+1,100+i
        special=width==1 and order
        g=a.view_group(1 if special else 0);v,w=column(g,e),column(g,h)
        value=9 if special else 42
        require(view_fill(C.byref(v),value,0))
        for i in range(width+1):
            selected=(i==1 if special else i<width)
            assert a[i,e]==(value if selected else i+1)
            assert a[i,h]==100+i
        require(view_copy(C.byref(w),C.byref(v),0))
        for i in range(width+1):
            selected=(i==1 if special else i<width)
            assert a[i,h]==(value if selected else 100+i)
        if width>1:
            last=column(a.view_group(1),e)
            require(view_fill(C.byref(last),-7,0));assert a[width,e]==-7


@case('grouped_views:GroupedView.FillCopyCrossGroupViewSources')
def cross(_):
    with Array(4,2,pattern=[(0,1)],group=4) as s,Array(4,2,pattern=[(0,1)],group=4) as d:
        s.initialize([1,2,3,4]);src=column(s.view_group(0,1),0);dst=column(d.view_group(0),0)
        require(view_copy(C.byref(dst),C.byref(src),0))
        near(d.values().ravel(),[1,2,3,4],rtol=0)
        require(view_fill(C.byref(dst),55.5,0))
        near(d.values().ravel(),[55.5]*4,rtol=0)
