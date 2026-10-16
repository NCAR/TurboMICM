"""Sparse array policies: CSR/CSC and grouped CSR/CSC, with the MICM fixtures."""
from contextlib import ExitStack
from arrays import *
from harness import case
import cases_matrix as dense

@case('sparse_matrix_standard_compressed_column_ordering:SparseMatrixBuilder.BadConfiguration',
      'sparse_matrix_standard_compressed_row_ordering:SparseMatrixBuilder.BadConfiguration')
def bad_configuration(_):
    pattern=(C.c_ubyte*9)()
    for r,c in ((3,0),(2,4),(6,7)):
        assert pattern_set(3,pattern,r,c)==-4
    assert list(pattern)==[0]*9
    require(pattern_set(3,pattern,1,2))
    assert list(pattern)==[0,0,0,0,0,1,0,0,0]

CHECKS={}
def check(*names):
    def register(fun):
        for name in names:CHECKS[name]=fun
        return fun
    return register


@check('ZeroMatrix','ConstZeroMatrix','SetScalar')
def zero(make):
    a=make(1,3,pattern=[])
    assert a.elements==0 and len(a.raw)==0
    assert array_element(a.handle,0,0)==-3
    for r,c in ((6,0),(1,3),(6,3)):assert array_element(a.handle,r,c)==-4
    require(array_fill(a.handle,2))
    assert array_set(a.handle,0,0,2)==-4


SMALL=[(0,1),(3,2),(2,3),(1,1),(2,1)]
@check('SingleBlockMatrix','ConstSingleBlockMatrix','MultiBlockMatrix')
def blocks(make):
    for count in (1,53):
        a=make(count,4,pattern=SMALL,initial=24 if count>1 else 0,integer=count==1)
        assert a.elements==5
        for r,c in SMALL:assert array_element(a.handle,r,c)>=0
        for r,c in ((2,2),(0,0),(3,3)):assert array_element(a.handle,r,c)==-3
        e=a.element(2,1)
        assert a[0,e]==(24 if count>1 else 0)
        a[0,e]=45;assert a[0,e]==45
        for r,c in ((4,2),(1,5)):assert array_element(a.handle,r,c)==-4
        assert array_set(a.handle,count,e,2)==-4
        if count>1:
            for b,r,c,v in ((43,3,2,29),(33,2,3,79),(2,2,3,64)):
                a[b,a.element(r,c)]=v
                assert a[b,a.element(r,c)]==v
        else:
            a[0,a.element(3,2)]=42;a[0,a.element(2,3)]=21
            cloned=array_clone(a.handle);assert cloned
            try:
                for r,c,v in ((2,1,45),(3,2,42),(2,3,21)):
                    value=D();require(array_get(cloned,0,a.element(r,c),C.byref(value)));assert value.value==v
                a[0,e]=99
                value=D();require(array_get(cloned,0,e,C.byref(value)));assert value.value==45
            finally:array_free(cloned)
        for b in range(min(3,count)):
            index=Z();require(array_index(a.handle,b,a.element(1,1),C.byref(index)))
            assert a.raw[index.value]==(24 if count>1 else 0)


@check('AddToDiagonal','Print','PrintNonZero')
def diagonal_print(make):
    a=make(3,4,pattern=SMALL+[(0,1)],initial=24,integer=True)
    for b in range(3):a[b,a.element(1,1)]=1+2*b
    require(array_add_diagonal(a.handle,7))
    expected=''
    for b in range(3):
        expected+=f'Block {b}\n0,24,0,0\n0,{8+2*b},0,0\n0,24,0,24\n0,0,24,0\n'
        for r,c in SMALL:assert a[b,a.element(r,c)]==(8+2*b if (r,c)==(1,1) else 24)
    n=array_format_dense(a.handle,None,0);buf=C.create_string_buffer(n+1)
    array_format_dense(a.handle,buf,n+1);assert buf.value.decode()==expected
    require(array_fill(a.handle,32))
    for b in range(3):a[b,a.element(1,1)]=2+2*b
    require(array_add_diagonal(a.handle,5))
    expected=''.join(f'Block {b}\n0, 1, 32\n1, 1, {7+2*b}\n2, 1, 32\n2, 3, 32\n3, 2, 32\n' for b in range(3))
    assert a.text()==expected


@check('ArrayFunction')
def function(make):
    pattern=[(0,0),(1,1),(2,2),(2,3),(3,3)]
    a=make(3,4,pattern=pattern,initial=1)
    ids=[a.element(r,c) for r,c in pattern]
    a[1,ids[1]]=2;a[2,ids[2]]=3;a[2,ids[3]]=4;a[2,ids[4]]=5
    def op(gs):
        g=gs[0];t=temporary(g)
        assign(t,[column(g,e) for e in ids[:4]],lambda a,b,c,d:a+b+c+d)
        assign(column(g,ids[3]),[t],lambda x:2*x)
    original=a.values().copy()
    with Function([a],op) as f:
        f(a)
        original[:,ids[3]]=[8,10,18];near(a.values(),original,rtol=0)
        b=make(3,4,pattern=pattern,initial=-1);f(b)
        expected=np.full((3,5),-1.);expected[:,ids[3]]=-8
        near(b.values(),expected,rtol=0)


@check('MultiMatrixArrayFunction')
def multi_function(make):
    a=make(3,3,pattern=[(0,1),(1,2)])
    b=make(3,3,pattern=[(0,0),(1,1),(2,2)])
    a.initialize([[i+10,2*i+20] for i in range(3)])
    b.initialize([[4*i,3*i,5*i] for i in range(3)])
    original=b.values().copy()
    def op(gs):
        t=temporary(gs[0]);assign(t,[column(gs[0],0),column(gs[1],2)],lambda x,y:x+y)
        assign(column(gs[0],1),[t],lambda x:x)
    with Function([a,b],op,[0,1]) as f:f(a,b)
    near(a.values(),[[10,10],[11,16],[12,22]],rtol=0)
    near(b.values(),original,rtol=0)


@check('MismatchedElementDimensions')
def missing_element(make):
    a=make(3,4,pattern=[(0,1),(1,1),(2,2)])
    assert array_element(a.handle,3,3)==-3
    g=a.view_group(0);v=View()
    assert group_sparse_view(C.byref(g),3,3,C.byref(v))==-3
    assert group_view(C.byref(g),a.elements,C.byref(v))==-4


@check('TwoSparseMatricesDifferentStructure')
def structures(make):
    a=make(4,3,pattern=[(0,1),(1,1),(2,2)])
    b=make(4,5,pattern=[(0,0),(1,2),(2,3),(3,4)])
    a.initialize([[i*2+1,i*3+2,i*5+3] for i in range(4)])
    b.initialize([[i+10,i+20,i+30,i+40] for i in range(4)])
    original=b.values().copy()
    with Function([a,b],lambda g:assign(column(g[0],2),[column(g[0],0),column(g[0],1),column(g[1],0)],lambda x,y,z:x+y+z),[0,1]) as f:
        f(a,b)
    near(a.values(),[[i*2+1,i*3+2,13+6*i] for i in range(4)],rtol=0)
    near(b.values(),original,rtol=0)


@check('SparseAndDenseMatrixFunction','SparseAndVectorMatrixFunction')
def dense_sparse(make):
    a=make(3,4,pattern=[(0,1),(1,2),(2,3),(3,3)])
    b=make(3,4,dense=True)
    a.initialize([[i+1,i+10,i+20,i+30] for i in range(3)])
    b.initialize([[100*i+10*j for j in range(4)] for i in range(3)])
    original=b.values().copy()
    def op(g):assign(column(g[0],3),[column(g[0],0),column(g[0],1),column(g[1],0),column(g[1],1)],lambda a,b,c,d:a+b+c+d)
    with Function([a,b],op,[0,1]) as f:f(a,b)
    near(a.values(),[[1,10,20,21],[2,11,21,223],[3,12,22,425]],rtol=0)
    near(b.values(),original,rtol=0)


@check('IncompatibleOrdering','IncompatibleVectorOrdering','IncompatibleSparseOrdering')
def incompatible(make):
    a=make(3,4,pattern=[(0,1),(1,2)],group=4)
    for width in (0,2,3,8):
        for dense_array in (False,True):
            b=make(3,4,pattern=None if dense_array else [(0,1),(1,2)],dense=dense_array,group=width)
            result=P();callback=Callback(lambda n,g,c:0)
            assert array_function_create(2,(P*2)(a.handle,b.handle),None,callback,None,C.byref(result))==-1
            assert not result


@check('MultipleSparseMatricesDifferentBlocksFromCreation','SparseMatrixVectorDifferentBlocksFromCreation')
def flexible(make):
    a,b=make(3,3),make(3,3)
    def op(g):assign(column(g[0],2),[column(g[0],0),column(g[0],1),column(g[1],2)],lambda a,b,c:a+b+c)
    with Function([a,b],op,[0,1]) as f:
        c,d=make(4,3),make(4,3)
        c.initialize([[i+1,i+10,0] for i in range(4)]);d.initialize([[0,0,i+100] for i in range(4)])
        f(c,d);near(c.values()[:,2],[111,114,117,120],rtol=0)
    a,v=make(3,2),make(3,vector=True)
    with Function([a,v],lambda g:assign(column(g[0],1),[column(g[0],0),column(g[1],0)],lambda x,y:x+y),[0,1]) as f:
        c,d=make(5,2),make(5,vector=True).initialize([10,20,30,40,50])
        for i in range(5):c[i,0]=i+1
        f(c,d);near(c.values()[:,1],[11,22,33,44,55],rtol=0)


@check('VectorInSparseMatrixFunction')
def vector_function(make):
    a=make(3,3).initialize([[i+1,i+2,i+3] for i in range(3)])
    v=make(3,vector=True).initialize([5,10,15])
    with Function([a,v],lambda g:assign(column(g[0],0),[column(g[1],0)],lambda x:2*x),[0,1]) as f:f(a,v)
    near(a.values(),[[10,2,3],[20,3,4],[30,4,5]],rtol=0)


@check('EmptySparseMatrixFunction','EmptyVectorEmptySparseMatrix')
def empty(make):
    a=make(3,4,pattern=[])
    with Function([a],lambda gs:None) as f:f(a)
    b,v=make(0,4,pattern=[(0,1)]),make(0,vector=True)
    def never(gs):raise AssertionError('callback on empty array')
    with Function([b,v],never) as f:f(b,v)


@check('GetBlockViewByVectorIndex')
def flat_view(make):
    a=make(7,4,pattern=SMALL,initial=2)
    for i in range(array_groups(a.handle)):
        g=a.view_group(i)
        for r,c in SMALL:
            v,w=View(),View()
            require(group_view(C.byref(g),a.element(r,c),C.byref(v)))
            require(group_sparse_view(C.byref(g),r,c,C.byref(w)))
            assert (v.data,v.stride,v.size)==(w.data,w.stride,w.size)
            require(view_fill(C.byref(v),10*r+c,1))
    near(a.values(),[[10*r+c for r,c in sorted(SMALL,key=lambda rc:a.element(*rc))]]*7,rtol=0)


# These policies differ only in whether a record is a dense row or a sparse
# block. Use the same fixtures and operations through sparse block views.
SHARED={
 'MultipleTemporaries':'MultipleTemporaries','BlockViewReuse':'ColumnViewReuse',
 'FunctionReusability':'FunctionReusability','ConstSparseMatrixFunction':'ConstMatrixFunction',
 'MismatchedBlockDimensions':'MismatchedRowDimensions','MismatchedBlocksAtInvocation':'MismatchedRowsAtInvocation',
 'MultipleSparseMatricesMismatchedBlocksAtInvocation':'MultipleMatricesMismatchedRowsAtInvocation',
 'VectorTooSmall':'VectorTooSmall','VectorTooLarge':'VectorTooLarge',
 'EmptyVectorNonEmptySparseMatrix':'EmptyVectorNonEmptyMatrix','NonEmptyVectorEmptySparseMatrix':'NonEmptyVectorEmptyMatrix',
 'MultipleVectorsDifferentSizes':'MultipleVectorsDifferentSizes','MultipleVectorsSameSize':'MultipleVectorsSameSize',
 'MultipleSparseMatricesOneVector':'MultipleMatricesOneVector',
 'MultipleSparseMatricesDifferentBlocksVector':'MultipleMatricesDifferentRowsVector',
 'VectorSizeMatchesOneSparseMatrixOnly':'VectorSizeMatchesOneMatrixOnly',
 'ConstVectorSparse':'ConstVector','MutableVectorSparse':'MutableVector',
 'FunctionReusabilityWithVectorsSparse':'FunctionReusabilityWithVectors',
 'FunctionInvocationWithWrongSizedVectorSparse':'FunctionInvocationWithWrongSizedVector',
 'MixedVectorBlockViewBlockVariable':'MixedVectorColumnViewRowVariable','IntegerVectorSparse':'IntegerVector',
 'FunctionWithConstSignatureSparse':'FunctionWithConstSignature','TestFill':'TestFill','TestCopy':'TestCopy'}
for name,source in SHARED.items():CHECKS[name]=dense.CHECKS[source]


def run_sparse(id_):
    name=id_.split('.')[-1]
    widths=(1,2,3,4) if '_vector_' in id_ else (0,)
    order=int('_column_' in id_)
    for width in widths:
        with ExitStack() as stack:
            def make(r,c=1,pattern=None,vector=False,dense=False,group=None,**kw):
                width_used=width if group is None else group
                if not vector and not dense:
                    if pattern is None:
                        pattern=[(i,i) for i in range(c)]
                    return stack.enter_context(Array(r,c,group=width_used,order=order,pattern=pattern,**kw))
                return stack.enter_context(Array(r,c,group=width_used,vector=vector,**kw))
            CHECKS[name](make)


PREFIXES=('sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.',
          'sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.',
          'sparse_matrix_vector_compressed_column_ordering:SparseVectorCompressedColumnMatrix.',
          'sparse_matrix_vector_compressed_row_ordering:SparseVectorCompressedRowMatrix.')
for name in CHECKS:
    for prefix in PREFIXES:
        if 'standard' in prefix and name in ('SparseAndVectorMatrixFunction','IncompatibleOrdering','IncompatibleVectorOrdering','IncompatibleSparseOrdering'):
            continue
        case(prefix+name)(run_sparse)
