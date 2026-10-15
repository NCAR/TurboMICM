"""Ports of the sparse matrix policies, with CSR/CSC storage assertions."""
# Copyright (C) 2023-2026 University Corporation for Atmospheric Research
# SPDX-License-Identifier: Apache-2.0
from sparseffi import *
from contextlib import ExitStack

CASES={}
POLICIES=[('sparse_matrix_standard_compressed_row_ordering','SparseCompressedRowMatrix',0,False),
          ('sparse_matrix_standard_compressed_column_ordering','SparseCompressedColumnMatrix',1,False),
          ('sparse_matrix_vector_compressed_row_ordering','SparseVectorCompressedRowMatrix',0,True),
          ('sparse_matrix_vector_compressed_column_ordering','SparseVectorCompressedColumnMatrix',1,True)]


def port(name,groups=(1,2,3,4)):
    def decorate(fn):
        for executable,suite,compression,vector in POLICIES:
            def run(compression=compression,vector=vector):
                for g in groups if vector else (0,):fn(compression,g)
            CASES[f'{executable}:{suite}.{name}']=run
        return fn
    return decorate


PATTERN=[(0,1),(3,2),(2,3),(1,1),(2,1)]


def check_pattern(s):
    assert s.nnz==5
    assert s.diagonal(0)==[s.index(0,1,1)]
    for row,col in [(0,1),(3,2),(2,1)]:assert s.flat(row,col)!=MISSING
    for row,col in [(2,2),(0,0),(3,3)]:assert s.flat(row,col)==MISSING


@port('ZeroMatrix',groups=(2,))
def zero(compression,g):
    with Sparse(3,compression=compression,group=g) as s:
        assert s.nnz==0 and s.values.size==0
        for row,col in [(0,0),(6,0),(1,3),(6,3)]:assert s.index(0,row,col)==MISSING
        scalar=D(2)
        for key in [(0,0,4),(1,0,0),(0,3,0),(0,1,1)]:assert lib.chem_sparse_set(s.handle,*key,C.byref(scalar))==1


@port('ConstZeroMatrix',groups=(3,))
def const_zero(compression,g):
    with Sparse(3,compression=compression,group=g) as original,original.clone() as s:
        assert s.nnz==0
        for row,col in [(0,0),(6,0),(1,3),(6,3)]:assert s.index(0,row,col)==MISSING


@port('SetScalar',groups=(3,))
def scalar(compression,g):
    with Sparse(3,compression=compression,group=g) as s:
        s.fill(2);assert list(s.data[:s.values.size])==[2]*s.values.size


@port('SingleBlockMatrix',groups=(4,))
def single(compression,g):
    with Sparse(4,entries=PATTERN,compression=compression,group=g) as s:
        check_pattern(s)
        assert s[0,2,1]==0;s[0,2,1]=45;assert s[0,2,1]==45
        width=g or 1
        for row,col,index,value in [(3,2,(3 if compression else 4)*width,42),(2,3,(4 if compression else 3)*width,21)]:
            assert lib.chem_sparse_single_index(s.handle,row,col)==index
            s.data[index]=value;assert s.data[index]==value and s[0,row,col]==value
        if g:
            assert s.values.group==4 and s.values.columns*g==20 and lib.chem_array_groups(s.values.handle)==1
        for row,col in [(4,2),(1,5),(2,2),(0,2)]:assert s.index(0,row,col)==MISSING


@port('ConstSingleBlockMatrix',groups=(2,))
def const_single(compression,g):
    with Sparse(4,entries=PATTERN,compression=compression,group=g,type=1) as original:
        original[0,2,1]=45;original[0,3,2]=42;original[0,2,3]=21
        with original.clone() as s:
            check_pattern(s);assert s[0,2,1]==45
            for row,col,index,value in [(3,2,(3 if compression else 4)*(g or 1),42),(2,3,(4 if compression else 3)*(g or 1),21)]:
                assert lib.chem_sparse_single_index(s.handle,row,col)==index
                assert s.data[index]==value
            if g:assert s.group==2 and s.values.columns*g==10 and lib.chem_array_groups(s.values.handle)==1


@port('MultiBlockMatrix',groups=(2,))
def multiple(compression,g):
    with Sparse(4,blocks=53,entries=PATTERN,compression=compression,group=g,value=24) as s:
        assert s.nnz*53==265 and s.nnz==5
        assert s[0,2,1]==24;s[0,2,1]=45;assert s[0,2,1]==45
        for key,value in [((43,3,2),29),((33,2,3),79),((2,2,3),64)]:s[key]=value;assert s[key]==value
        for block in range(3):assert s.diagonal(block)==[s.index(block,1,1)]
        assert lib.chem_sparse_single_index(s.handle,0,1)==MISSING
        first_index=(4 if compression else 3)*(g or 1)
        assert s.index(0,2,3)==first_index;s.data[first_index]=21;assert s.data[first_index]==21
        expected=14 if g else 12
        assert s.index(2,2,1)==expected;s.data[expected]=31;assert s.data[expected]==31
        if g:assert s.group==2 and s.values.columns*g==10 and (4+g-1)//g==2


@port('AddToDiagonal')
def diagonal(compression,g):
    with Sparse(4,blocks=3,entries=PATTERN+[(0,1)],compression=compression,group=g,type=1,value=24) as s:
        assert s.nnz==5
        for block,value in enumerate([1,3,5]):s[block,1,1]=value
        s.add_diagonal(7)
        for block,value in enumerate([8,10,12]):
            assert s[block,1,1]==value
            for row,col in [(0,1),(2,1),(2,3),(3,2)]:assert s[block,row,col]==24


@port('Print')
def print_dense(compression,g):
    with Sparse(4,blocks=3,entries=PATTERN+[(0,1)],compression=compression,group=g,type=1,value=24) as s:
        for block,value in enumerate([1,3,5]):s[block,1,1]=value
        s.add_diagonal(7)
        expected=''.join(f'Block {b}\n0,24,0,0\n0,{v},0,0\n0,24,0,24\n0,0,24,0\n' for b,v in enumerate([8,10,12]))
        assert s.format()==expected,(s.format(),expected)


@port('PrintNonZero')
def print_nonzeros(compression,g):
    with Sparse(4,blocks=3,entries=PATTERN,compression=compression,group=g,type=1,value=32) as s:
        for block,value in enumerate([2,4,6]):s[block,1,1]=value
        s.add_diagonal(5)
        expected=''.join(f'Block {b}\n0, 1, 32\n1, 1, {v}\n2, 1, 32\n2, 3, 32\n3, 2, 32\n' for b,v in enumerate([7,9,11]))
        assert s.format(True)==expected,(s.format(True),expected)


def bad_configuration(compression):
    for row,col in [(3,0),(2,4),(6,7)]:
        with Sparse(3,entries=[(row,col)],compression=compression,status=1):pass
for executable,suite,compression,vector in POLICIES:
    if not vector:CASES[f'{executable}:SparseMatrixBuilder.BadConfiguration']=lambda c=compression:bad_configuration(c)


def make(compression,g,order=4,blocks=3,entries=((0,1),(1,2),(2,3)),value=0,type=0):
    return Sparse(order,blocks,entries,g,compression,type,value)


@port('ArrayFunction')
def array_function(c,g):
    pattern=[(0,0),(1,1),(2,2),(2,3),(3,3)]
    with make(c,g,entries=pattern,value=1) as s:
        for key,value in [((0,0,0),1),((1,1,1),2),((2,2,2),3),((2,2,3),4),((2,3,3),5)]:s[key]=value
        ids=[s.flat(*p) for p in pattern]
        def body(a):
            with Array(a.rows,1,g) as tmp:
                assign(tmp.col(0),[a.col(i,True) for i in ids[:4]],lambda x,y,z,w:x+y+z+w)
                assign(a.col(ids[3]),[tmp.col(0,True)],lambda x:2*x)
        with Function([s.values],body) as f:
            f(s.values)
            assert [s[b,2,3] for b in range(3)]==[8,10,18]
            for b in range(3):
                assert s[b,0,0]==1 and s[b,1,1]==(2 if b==1 else 1) and s[b,2,2]==(3 if b==2 else 1) and s[b,3,3]==(5 if b==2 else 1)
            with make(c,g,entries=pattern,value=-1) as t:
                f(t.values)
                for b in range(3):
                    assert t[b,2,3]==-8
                    for row,col in [(0,0),(1,1),(2,2),(3,3)]:assert t[b,row,col]==-1


@port('MultiMatrixArrayFunction')
def multi_function(c,g):
    with make(c,g,3,entries=[(0,1),(1,2)]) as s,make(c,g,3,entries=[(0,0),(1,1),(2,2)]) as t:
        for b in range(3):
            s[b,0,1]=b+10;s[b,1,2]=b*2+20
            t[b,0,0]=b*4;t[b,1,1]=b*3;t[b,2,2]=b*5
        i,j,k=s.flat(0,1),s.flat(1,2),t.flat(2,2)
        def body(a,b):assign(a.col(j),[a.col(i,True),b.col(k,True)],lambda x,y:x+y)
        with Function([s.values,t.values],body,readonly=[0,1]) as f:f(s.values,t.values)
        assert [s[b,1,2] for b in range(3)]==[10,16,22]
        assert [s[b,0,1] for b in range(3)]==[10,11,12]
        assert [t[b,0,0] for b in range(3)]==[0,4,8]


@port('MismatchedBlockDimensions')
def mismatch_blocks(c,g):
    with make(c,g,3,3,[(0,1),(1,1)]) as s,make(c,g,3,4,[(0,1),(1,1)]) as t:
        i,j=s.flat(0,1),s.flat(1,1)
        def body(a,b):assign(a.col(j),[a.col(i,True),b.col(i,True)],lambda x,y:x+y)
        with Function([s.values,s.values],body,readonly=[0,1]) as f:f(s.values,s.values);f(s.values,t.values,status=1)


@port('MismatchedElementDimensions')
def missing_element(c,g):
    with make(c,g,4,3,[(0,1),(1,1),(2,2)]) as s:
        assert s.index(0,3,3)==MISSING
        out=View();assert lib.chem_sparse_block_view(s.handle,3,3,0,C.byref(out))==1


@port('MultipleSparseMatricesDifferentBlocksFromCreation')
def different_blocks(c,g):
    with make(c,g) as s,make(c,g) as t:
        i,j,k=[s.flat(*p) for p in [(0,1),(1,2),(2,3)]]
        for b in range(3):s[b,0,1]=b+1;s[b,1,2]=b+10;t[b,2,3]=b+100
        def body(a,b):assign(a.col(k),[a.col(i,True),a.col(j,True),b.col(k,True)],lambda x,y,z:x+y+z)
        with Function([s.values,t.values],body,readonly=[0,1]) as f,make(c,g,blocks=4) as u,make(c,g,blocks=4) as v:
            for b in range(4):u[b,0,1]=b+1;u[b,1,2]=b+10;v[b,2,3]=b+100
            f(u.values,v.values);assert [u[b,2,3] for b in range(4)]==[111,114,117,120]


@port('SparseMatrixVectorDifferentBlocksFromCreation')
def different_vector_blocks(c,g):
    with make(c,g,3,3,[(0,1),(1,2)]) as s,Array.from_rows([[1],[2],[3]],g) as v:
        i,j=s.flat(0,1),s.flat(1,2)
        def body(a,v):assign(a.col(j),[a.col(i,True),v.col(0,True)],lambda x,y:x+y)
        with Function([s.values,v],body,readonly=[0,1]) as f,make(c,g,3,5,[(0,1),(1,2)]) as t,Array.from_rows([[10],[20],[30],[40],[50]],g) as w:
            for b in range(5):t[b,0,1]=b+1
            f(t.values,w);assert [t[b,1,2] for b in range(5)]==[11,22,33,44,55]


@port('MismatchedBlocksAtInvocation')
def bad_blocks_at_call(c,g):
    with make(c,g,3,3,[(0,1),(1,2)]) as s,make(c,g,3,4,[(0,1),(1,2)]) as t:
        i,j=s.flat(0,1),s.flat(1,2)
        with Function([s.values,s.values],lambda a,b:assign(a.col(j),[b.col(i,True)],lambda x:2*x),readonly=[0,1]) as f:
            f(s.values,s.values);f(s.values,t.values,status=1);f(t.values,s.values,status=1)


@port('MultipleSparseMatricesMismatchedBlocksAtInvocation')
def bad_multiple_blocks_at_call(c,g):
    with make(c,g) as a,make(c,g) as b,make(c,g,blocks=4) as d:
        i,j,k=[a.flat(*p) for p in [(0,1),(1,2),(2,3)]]
        def body(a,b,c):assign(a.col(k),[a.col(i,True),b.col(j,True),c.col(k,True)],lambda x,y,z:x+y+z)
        with Function([a.values,b.values,b.values],body,readonly=[0,1,1]) as f:
            f(a.values,b.values,b.values);f(a.values,b.values,d.values,status=1);f(d.values,b.values,b.values,status=1)


@port('MultipleTemporaries')
def temporaries(c,g):
    coords=[(0,1),(1,2),(2,3),(3,4),(4,4)]
    with make(c,g,5,4,coords) as s:
        for b in range(4):s[b,0,1]=b+1;s[b,1,2]=(b+1)*10
        ids=[s.flat(*p) for p in coords]
        def body(a):
            with Array(a.rows,1,g) as p,Array(a.rows,1,g) as t:
                assign(p.col(0),[a.col(ids[0],True),a.col(ids[1],True)],lambda x,y:x*y)
                assign(t.col(0),[a.col(ids[0],True),a.col(ids[1],True)],lambda x,y:x+y)
                for index,fn in zip(ids[2:],[lambda x,y:x+y,lambda x,y:x-y,lambda x,y:x*y]):assign(a.col(index),[p.col(0,True),t.col(0,True)],fn)
        with Function([s.values],body) as f:f(s.values)
        for b,want in [(0,[21,-1,110]),(1,[62,18,880]),(3,[204,116,7040])]:assert [s[b,*p] for p in coords[2:]]==want


@port('BlockViewReuse')
def view_reuse(c,g):
    coords=[(0,1),(1,2),(2,3),(3,3)]
    with make(c,g,entries=coords) as s:
        for b in range(3):s[b,0,1]=b+1
        ids=[s.flat(*p) for p in coords]
        def body(m):
            a,b,c,d=[m.col(i) for i in ids]
            assign(b,[a],lambda x:x*2);assign(c,[a,b],lambda x,y:x+y);assign(d,[c,b],lambda x,y:x*y)
        with Function([s.values],body) as f:f(s.values)
        for b,want in enumerate([[2,3,6],[4,6,24],[6,9,54]]):assert [s[b,*p] for p in coords[1:]]==want


@port('FunctionReusability')
def function_reuse(c,g):
    coords=[(0,0),(1,1),(2,2)]
    with make(c,g,3,2,coords) as s:
        ids=[s.flat(*p) for p in coords]
        for b in range(2):
            for j in range(3):s[b,j,j]=b+j
        def body(a):assign(a.col(ids[2]),[a.col(i,True) for i in ids],lambda x,y,z:2*(x+y+z))
        with Function([s.values],body) as f:
            f(s.values);assert [s[b,2,2] for b in range(2)]==[6,12]
            with make(c,g,3,2,coords,value=5) as t:f(t.values);assert [t[b,2,2] for b in range(2)]==[30,30]
            with make(c,g,3,2,coords) as u:u[1,0,0]=10;f(u.values);assert [u[b,2,2] for b in range(2)]==[0,20]


@port('TwoSparseMatricesDifferentStructure')
def different_structure(c,g):
    with make(c,g,3,4,[(0,1),(1,1),(2,2)]) as a,make(c,g,5,4,[(0,0),(1,2),(2,3),(3,4)]) as b:
        for k in range(4):
            for coord,v in [((0,1),2*k+1),((1,1),3*k+2),((2,2),5*k+3)]:a[k,*coord]=v
            for coord,v in [((0,0),k+10),((1,2),k+20),((2,3),k+30),((3,4),k+40)]:b[k,*coord]=v
        i,j,k,l=a.flat(0,1),a.flat(1,1),a.flat(2,2),b.flat(0,0)
        with Function([a.values,b.values],lambda x,y:assign(x.col(k),[x.col(i,True),x.col(j,True),y.col(l,True)],lambda a,b,c:a+b+c),readonly=[0,1]) as f:f(a.values,b.values)
        assert [a[n,2,2] for n in range(4)]==[13,19,25,31] and [a[n,0,1] for n in range(4)]==[1,3,5,7]
        assert [b[0,0,0],b[1,1,2],b[2,2,3],b[3,3,4]]==[10,21,32,43]


def mixed_dense(c,g,dense_group):
    with make(c,g,entries=[(0,1),(1,2),(2,3),(3,3)]) as s,Array(3,4,dense_group) as d:
        for b in range(3):
            for coord,v in [((0,1),b+1),((1,2),b+10),((2,3),b+20),((3,3),b+30)]:s[b,*coord]=v
            for j in range(4):d[b,j]=b*100+j*10
        i,j,k=[s.flat(*p) for p in [(0,1),(1,2),(3,3)]]
        with Function([s.values,d],lambda m,d:assign(m.col(k),[m.col(i,True),m.col(j,True),d.col(0,True),d.col(1,True)],lambda a,b,c,d:a+b+c+d),readonly=[0,1]) as f:f(s.values,d)
        assert [s[b,3,3] for b in range(3)]==[21,223,425]


@port('SparseAndDenseMatrixFunction',groups=(1,))
def sparse_dense(c,g):mixed_dense(c,g,0)


@port('SparseAndVectorMatrixFunction',groups=(1,2,3,4,10))
def sparse_vector_matrix(c,g):mixed_dense(c,g,g)


@port('ConstSparseMatrixFunction')
def const_sparse(c,g):
    with make(c,g,entries=[(0,1),(1,2),(2,3),(3,3)]) as s:
        for b in range(3):s[b,0,1]=b+1;s[b,1,2]=b+10
        i,j,k=[s.flat(*p) for p in [(0,1),(1,2),(2,3)]]
        def body(a):
            with Array(a.rows,1,g) as temp:assign(temp.col(0),[a.col(i,True),a.col(j,True)],lambda x,y:x+y)
            reduce([a.col(k,True)],lambda acc,x:acc+x)
        with Function([s.values],body,readonly=[1]) as f:f(s.values,readonly=[1])
        assert s[0,0,1]==1 and s[1,1,2]==11


@port('EmptySparseMatrixFunction')
def empty_sparse(c,g):
    with make(c,g,entries=[]) as s,Function([s.values],lambda a:None) as f:f(s.values)
    with make(c,g,blocks=0,entries=[(0,1)]) as s,Function([s.values],lambda a:None) as f:f(s.values)


def reject_vectors(c,g,blocks,lengths,other_blocks=None):
    with ExitStack() as stack:
        s=stack.enter_context(make(c,g,3,blocks,[(0,1)]));arrays=[s.values]
        if other_blocks is not None:arrays.append(stack.enter_context(make(c,g,3,other_blocks,[(0,1)])).values)
        arrays.extend(stack.enter_context(Array(n,1,g)) for n in lengths)
        def forbidden(*args):raise AssertionError('invalid dimensions reached callback')
        with Function(arrays,forbidden) as f:f(*arrays,status=1)


@port('VectorTooSmall')
def small_vector(c,g):reject_vectors(c,g,3,[2])
@port('VectorTooLarge')
def large_vector(c,g):reject_vectors(c,g,3,[4])
@port('EmptyVectorNonEmptySparseMatrix')
def empty_vector(c,g):reject_vectors(c,g,3,[0])
@port('NonEmptyVectorEmptySparseMatrix')
def empty_matrix(c,g):reject_vectors(c,g,0,[3])
@port('MultipleVectorsDifferentSizes')
def wrong_vectors(c,g):reject_vectors(c,g,3,[3,2])
@port('MultipleSparseMatricesDifferentBlocksVector')
def different_blocks_vector(c,g):reject_vectors(c,g,3,[3],other_blocks=4)
@port('VectorSizeMatchesOneSparseMatrixOnly')
def wrong_vector_pair(c,g):reject_vectors(c,g,3,[3,2],other_blocks=3)


@port('EmptyVectorEmptySparseMatrix')
def both_empty(c,g):
    with make(c,g,3,0,[(0,1)]) as s,Array(0,1,g) as v:
        with Function([s.values,v],lambda *args:None,readonly=[0,1]) as f:f(s.values,v)


@port('VectorInSparseMatrixFunction')
def vector_function(c,g):
    with make(c,g) as s,Array.from_rows([[5],[10],[15]],g) as v:
        for b in range(3):s[b,0,1]=b+1;s[b,1,2]=b+2;s[b,2,3]=b+3
        idx=s.flat(0,1)
        with Function([s.values,v],lambda m,v:assign(m.col(idx),[v.col(0,True)],lambda x:x*2),readonly=[0,1]) as f:f(s.values,v)
        assert [s[b,0,1] for b in range(3)]==[10,20,30]


@port('MultipleVectorsSameSize')
def multiple_vectors(c,g):
    with make(c,g,3,3,[(0,1)]) as s,Array.from_rows([[1],[2],[3]],g) as v,Array.from_rows([[4],[5],[6]],g) as w:
        idx=s.flat(0,1)
        with Function([s.values,v,w],lambda m,v,w:assign(m.col(idx),[v.col(0,True),w.col(0,True)],lambda x,y:x+y),readonly=[0,1,1]) as f:f(s.values,v,w)
        assert [s[b,0,1] for b in range(3)]==[5,7,9]


@port('MultipleSparseMatricesOneVector')
def multiple_matrices_vector(c,g):
    with make(c,g,3,3,[(0,1)]) as a,make(c,g,3,3,[(0,1)]) as b,Array.from_rows([[10],[20],[30]],g) as v:
        for block in range(3):a[block,0,1]=block+1
        idx=a.flat(0,1)
        with Function([a.values,b.values,v],lambda a,b,v:assign(b.col(idx),[a.col(idx,True),v.col(0,True)],lambda x,y:x+y),readonly=[1,0,1]) as f:f(a.values,b.values,v)
        assert [b[block,0,1] for block in range(3)]==[11,22,33]


@port('ConstVectorSparse')
def const_vector(c,g):
    with make(c,g,3,3,[(0,1)]) as s,Array.from_rows([[5],[10],[15]],g) as v:
        idx=s.flat(0,1)
        with Function([s.values,v],lambda m,v:assign(m.col(idx),[v.col(0,True)],lambda x:3*x),readonly=[0,1]) as f:f(s.values,v,readonly=[0,1])
        assert [s[b,0,1] for b in range(3)]==[15,30,45]


@port('MutableVectorSparse')
def mutable_vector(c,g):
    with make(c,g,3,3,[(0,1)]) as s,Array.from_rows([[5],[10],[15]],g) as v:
        for b in range(3):s[b,0,1]=b+1
        idx=s.flat(0,1)
        with Function([s.values,v],lambda m,v:assign(v.col(0),[m.col(idx,True)],lambda x:3*x),readonly=[1,0]) as f:f(s.values,v)
        assert [v[b] for b in range(3)]==[3,6,9]


@port('FunctionReusabilityWithVectorsSparse')
def vector_reuse(c,g):
    with make(c,g,3,3,[(0,1)]) as s,make(c,g,3,3,[(0,1)]) as t,Array.from_rows([[1],[2],[3]],g) as v,Array.from_rows([[10],[20],[30]],g) as w:
        idx=s.flat(0,1)
        with Function([s.values,v],lambda m,v:assign(m.col(idx),[v.col(0,True)],lambda x:2*x),readonly=[0,1]) as f:
            f(s.values,v);f(t.values,w)
        assert [s[b,0,1] for b in range(3)]==[2,4,6] and [t[b,0,1] for b in range(3)]==[20,40,60]


@port('FunctionInvocationWithWrongSizedVectorSparse')
def vector_wrong_invocation(c,g):
    with make(c,g,3,3,[(0,1)]) as s,Array.from_rows([[1],[2],[3]],g) as v,Array.from_rows([[1],[2]],g) as w:
        with Function([s.values,v],lambda *args:None,readonly=[1,1]) as f:f(s.values,w,status=1)


@port('MixedVectorBlockViewBlockVariable')
def mixed_vector(c,g):
    with make(c,g) as s,Array.from_rows([[100],[200],[300]],g) as v:
        for b in range(3):s[b,0,1]=b+1;s[b,1,2]=b+2;s[b,2,3]=b+3
        i,j,k=[s.flat(*p) for p in [(0,1),(1,2),(2,3)]]
        def body(a,v):
            with Array(a.rows,1,g) as temp:
                assign(temp.col(0),[a.col(i,True),a.col(j,True),v.col(0,True)],lambda x,y,z:x+y+z)
                assign(a.col(k),[temp.col(0,True)],lambda x:x*2)
        with Function([s.values,v],body,readonly=[0,1]) as f:f(s.values,v)
        assert [s[b,2,3] for b in range(3)]==[206,410,614]


@port('IntegerVectorSparse')
def integer_vector(c,g):
    with make(c,g,3,3,[(0,1)]) as s,Array.from_rows([[5],[10],[15]],g,type=1) as v:
        idx=s.flat(0,1)
        with Function([s.values,v],lambda m,v:assign(m.col(idx),[v.col(0,True)],lambda x:float(x)*1.5),readonly=[0,1]) as f:f(s.values,v)
        assert [s[b,0,1] for b in range(3)]==[7.5,15,22.5]


@port('FunctionWithConstSignatureSparse')
def const_signature(c,g):
    with make(c,g,2,3,[(0,0)]) as s,Array.from_rows([[1],[2],[3]],g) as v:
        idx=s.flat(0,0)
        with Function([s.values,v],lambda m,v:assign(m.col(idx),[v.col(0,True)],lambda x:x*2),readonly=[0,1]) as f:f(s.values,v,readonly=[0,1])


@port('GetBlockViewByVectorIndex')
def index_view(c,g):
    coords=[(0,0),(0,2),(1,0),(1,1),(2,1),(2,2)]
    with make(c,g,3,3,coords) as s,Array.from_rows([[0],[1],[2]],g,type=1) as values:
        indices=[s.index(0,*p) for p in coords]
        def body(m,v):
            for k,index in enumerate(indices):assign(s.index_view(index),[v.col(0,True)],lambda x:x+10*k)
        with Function([s.values,values],body,readonly=[0,1]) as f:f(s.values,values)
        for block in range(3):
            for k,coord in enumerate(coords):assert s[block,*coord]==block+k*10


@port('TestFill')
def fills(c,g):
    with make(c,g,3,3,[(0,1),(2,0)]) as s:
        i,j=s.flat(0,1),s.flat(2,0)
        with Function([s.values],lambda a:fill(a.col(i),3.2)) as f:f(s.values)
        assert [s[b,0,1] for b in range(3)]==[3.2]*3 and [s[b,2,0] for b in range(3)]==[0]*3
        with Array(3,1,g) as v,Function([s.values,v],lambda m,v:fill(v.col(0),3.2),readonly=[1,0]) as f:
            f(s.values,v);assert [v[b] for b in range(3)]==[3.2]*3
        def body(a):
            with Array(a.rows,1,g) as temp:fill(temp.col(0),9.9);assign(a.col(j),[temp.col(0,True)],lambda x:x)
        with Function([s.values],body) as f:f(s.values)
        assert [s[b,2,0] for b in range(3)]==[9.9]*3


@port('TestCopy')
def copies(c,g):
    def fresh():
        s=make(c,g,3,3,[(0,1),(2,0)])
        for b,v in enumerate([3.2,4.2,1.3]):s[b,0,1]=v
        return s
    with fresh() as s,Array.from_rows([[10],[20],[30]],g) as v:
        i=s.flat(0,1)
        for b in range(3):s[b,0,1]=0
        with Function([s.values,v],lambda a,v:copy(a.col(i),v.col(0,True)),readonly=[0,1]) as f:f(s.values,v)
        assert [s[b,0,1] for b in range(3)]==[10,20,30] and s[0,2,0]==0
    with fresh() as s,Array(3,1,g,-1) as v:
        i=s.flat(0,1)
        with Function([s.values,v],lambda a,v:copy(v.col(0),a.col(i,True)),readonly=[1,0]) as f:f(s.values,v)
        assert [v[b] for b in range(3)]==[3.2,4.2,1.3]
    for ro in [False,True]:
        with fresh() as s:
            i,j=s.flat(0,1),s.flat(2,0)
            with Function([s.values],lambda a:copy(a.col(j),a.col(i,ro))) as f:f(s.values)
            assert [s[b,2,0] for b in range(3)]==[3.2,4.2,1.3] and s[0,0,1]==3.2
    with fresh() as s:
        i=s.flat(0,1)
        def body(a):
            with Array(a.rows,1,g) as tmp:copy(tmp.col(0),a.col(i,True));fill(a.col(i),0);assign(a.col(i),[tmp.col(0,True)],lambda x:x)
        with Function([s.values],body) as f:f(s.values)
        assert [s[b,0,1] for b in range(3)]==[3.2,4.2,1.3]
    with fresh() as s,Array(3,2,g) as d:
        for b in range(3):d[b,1]=(b+1)*100
        j=s.flat(2,0)
        with Function([s.values,d],lambda a,d:assign(a.col(j),[d.col(1,True)],lambda x:x),readonly=[0,1]) as f:f(s.values,d)
        assert [s[b,2,0] for b in range(3)]==[100,200,300]
    with fresh() as s,Array(3,2,g) as d,Array(3,1,g) as tmp:
        i=s.flat(0,1)
        with Function([s.values,tmp],lambda a,t:copy(t.col(0),a.col(i,True)),readonly=[1,0]) as f:f(s.values,tmp)
        with Function([d,tmp],lambda d,t:assign(d.col(0),[t.col(0,True)],lambda x:x),readonly=[0,1]) as f:f(d,tmp)
        assert [d[b,0] for b in range(3)]==[3.2,4.2,1.3]
    with fresh() as s,make(c,g,3,3,[(1,2)]) as t:
        i,j=s.flat(0,1),t.flat(1,2)
        with Function([t.values,s.values],lambda t,s:copy(t.col(j),s.col(i,True)),readonly=[0,1]) as f:f(t.values,s.values)
        assert [t[b,1,2] for b in range(3)]==[3.2,4.2,1.3]


def incompatible_dense(c):
    for g in [2,3,4]:
        with make(c,g,entries=[(0,1),(1,2)]) as s,Array(3,4,0) as d,Function([s.values,d],lambda *a:None,status=1):pass


def incompatible_vector(c):
    for g,other in [(2,1),(2,3),(4,2),(4,10)]:
        with make(c,g,entries=[(0,1),(1,2)]) as s,Array(3,4,other) as d,Function([s.values,d],lambda *a:None,status=1):pass


def incompatible_sparse(c):
    for g,other in [(1,2),(2,3),(3,4),(4,10)]:
        with make(c,g,entries=[(0,1),(1,2)]) as s,make(c,other,entries=[(0,1),(1,2)]) as d,Function([s.values,d.values],lambda *a:None,status=1):pass

for executable,suite,compression,vector in POLICIES:
    if vector:
        for name,fn in [('IncompatibleOrdering',incompatible_dense),('IncompatibleVectorOrdering',incompatible_vector),('IncompatibleSparseOrdering',incompatible_sparse)]:
            CASES[f'{executable}:{suite}.{name}']=lambda fn=fn,c=compression:fn(c)
