"""Ports of the shared dense/vector-matrix policies and their storage checks.

The callbacks reproduce the user-supplied test lambdas. Storage, iteration,
views, ownership, validation, copying, and reductions execute in the C library.
"""
# Copyright (C) 2023-2026 University Corporation for Atmospheric Research
# SPDX-License-Identifier: Apache-2.0
from contextlib import ExitStack
import math
import sys
from arrayffi import *

CASES={}


def port(name,vector=None,groups=(1,2,3,4)):
    def decorate(fn):
        CASES['Matrix.'+name]=lambda:fn(0)
        def vector_case():
            for g in groups:fn(g)
        CASES['VectorMatrix.'+(vector or name)]=vector_case
        return fn
    return decorate


def near(a,b,tol=1e-5):assert abs(a-b)<=tol,(a,b,tol)

def check_rows(a,values):assert a.values()==values,(a.values(),values)


@port('SmallMatrix','SmallVectorMatrix',groups=(2,))
def small(g):
    with Array(3,5,g) as a:
        a[1,3]=64.7;a[0,0]=41.2;a[2,4]=102.3
        assert [a[1,3],a[0,0],a[2,4]]==[64.7,41.2,102.3]
        if not g:assert a.size==15 and [a.data[0],a.data[14],a.data[8]]==[41.2,102.3,64.7]
        else:
            assert a.size==20 and lib.chem_array_groups(a.handle)==2 and g*a.columns==10
            assert [a.data[0],a.data[18],a.data[7]]==[41.2,102.3,64.7]


@port('SmallConstMatrix','SmallConstVectorMatrix',groups=(4,))
def small_const(g):
    with Array(3,5,g) as a:
        a[1,3]=64.7;a[0,0]=41.2;a[2,4]=102.3
        with a.clone() as b:
            assert [b[1,3],b[0,0],b[2,4]]==[64.7,41.2,102.3]
            if not g:assert b.size==15 and [b.data[0],b.data[14],b.data[8]]==[41.2,102.3,64.7]
            else:
                assert b.size==20 and lib.chem_array_groups(b.handle)==1 and b.group==4
                assert [b.data[0],b.data[18],b.data[13]]==[41.2,102.3,64.7]


@port('InitializeMatrix','InitializeVectorMatrix')
def initialize(g):
    with Array(2,3,g,12.4) as a:assert [a[0,0],a[1,0],a[1,2]]==[12.4]*3


@port('InitializeConstMatrix','InitializeConstVectorMatrix')
def initialize_const(g):initialize(g)


@port('LoopOverMatrix','LoopOverVectorMatrix')
def loop(g):
    with Array(3,4,g,type=1) as a:
        for i in range(a.rows):
            for j in range(a.columns):a[i,j]=i*100+j
        assert [a[0,0],a[1,2],a[2,3],a[0,3]]==[0,102,203,3]


@port('LoopOverConstMatrix','LoopOverConstVectorMatrix')
def const_loop(g):
    with Array(3,4,g,type=1) as a:
        for i in range(3):
            for j in range(4):a[i,j]=i*100+j
        with a.clone() as b:assert [b[0,0],b[1,2],b[2,3],b[0,3]]==[0,102,203,3]


def iterate():
    with Array(3,4,type=1) as a:
        row=a.row(1)
        for j in range(row.length):element(row,j)[0]=42+j
        with a.clone() as b:
            assert b.read_row(1)==[42,43,44,45]
            assert [a[0,3],a[1,0],a[1,1],a[1,3],a[2,0]]==[0,42,43,45,0]
CASES['Matrix.IterateOverMatrix']=iterate


@port('Strides',groups=(3,4))
def strides(g):
    with Array(3,4,g,type=1) as a:
        rowstride=1 if g else 4;colstride=g or 1
        assert a.row(0).element_stride==colstride
        assert (a.col(0).element_stride if g else a.col(0).block_stride)==rowstride
        for i in range(3):
            for j in range(4):a.data[i*rowstride+j*colstride]=i*100+j
        check_rows(a,[[i*100+j for j in range(4)] for i in range(3)])


@port('ConversionToVector')
def to_vector(g):
    with Array(2,3,g) as a:
        for j,v in enumerate([13.2,31.2,314.2]):a[1,j]=v
        assert a.read_row(1)==[13.2,31.2,314.2]


@port('ConstConversionToVector')
def const_to_vector(g):
    with Array(2,3,g) as a:
        ok(a.assign_row(1,[13.2,31.2,314.2]))
        with a.clone() as b:assert b.read_row(1)==[13.2,31.2,314.2]


@port('ConversionFromVector')
def from_vector(g):
    with Array.from_rows([],g) as empty:assert empty.rows==0
    rows=[[412.3,32.4,41.3],[5.33,-.3,31.2]]
    with Array.from_rows(rows,g) as a:
        assert a.rows==2 and a.columns==3 and a.row(0).length==3 and a.row(1).length==3
        check_rows(a,rows)
    with Array.from_rows([[3],[4,5],[5]],g,type=1,status=1):pass


@port('AssignmentFromVector')
def row_assignment(g):
    with Array(4,3,g) as a:
        ok(a.assign_row(2,[12.3,15.1,24.3]));check_rows(a,[[0,0,0],[0,0,0],[12.3,15.1,24.3],[0,0,0]])
        ok(a.assign_row(2,[14.3,52.3,65.7,16.34]));check_rows(a,[[0,0,0],[0,0,0],[14.3,52.3,65.7],[0,0,0]])
        assert a.assign_row(2,[13.2,52.8])==1


@port('Axpy')
def axpy(g):
    with Array(4,3,g,100) as a,Array(4,3,g,200) as b:
        expected=0
        for i in range(4):
            for j in range(3):
                a[i,j]=i*10.3+j*100.5;b[i,j]=i*1.7+j*10.2;expected+=a[i,j]+1.39*b[i,j]
        ok(lib.chem_array_axpy(a.handle,1.39,b.handle))
        near(sum(sum(row) for row in a.values()),expected)


@port('ForEach')
def for_each(g):
    with Array(4,3,g) as a,Array(4,3,g) as b,Array(4,3,g) as c:
        for i in range(4):
            for j in range(3):a[i,j]=i*10.3+j*100.5;b[i,j]=i*1.7+j*10.2;c[i,j]=i*19.5+j*32.2
        for j in range(3):apply([a.col(j),b.col(j,True)],lambda p:p[0].__setitem__(0,p[0][0]+p[1][0]))
        for i in range(4):
            for j in range(3):near(a[i,j],(i*10.3+j*100.5)+(i*1.7+j*10.2));a[i,j]=i*10.3+j*100.5
        for j in range(3):apply([a.col(j),b.col(j,True),c.col(j,True)],lambda p:p[0].__setitem__(0,p[0][0]+p[1][0]-p[2][0]))
        for i in range(4):
            for j in range(3):near(a[i,j],(i*10.3+j*100.5)+(i*1.7+j*10.2)-(i*19.5+j*32.2))


@port('SetScaler')
def scalar(g):
    with Array(2,3,g) as a:a.fill(2);assert list(a.data[:a.size])==[2]*a.size


@port('Max')
def max_value(g):
    with Array(2,3,g) as a:
        ok(lib.chem_array_limit(a.handle,2,1));assert list(a.data[:a.size])==[2]*a.size
        a.fill(1);a[1,1]=3;ok(lib.chem_array_limit(a.handle,2,1));check_rows(a,[[2,2,2],[2,3,2]])


@port('Min')
def min_value(g):
    with Array(2,3,g) as a:
        ok(lib.chem_array_limit(a.handle,2,0));assert list(a.data[:a.size])==[0]*a.size
        a.fill(1);a[1,1]=3;ok(lib.chem_array_limit(a.handle,2,0));check_rows(a,[[1,1,1],[1,2,1]])


@port('Print')
def print_array(g):
    with Array(2,3,g) as a:a[1,1]=3;assert a.format()=='0,0,0\n0,3,0\n'


def sum_columns(m,scale=4):
    with Array(m.rows,1,m.group) as tmp:
        assign(tmp.col(0),[m.col(j,True) for j in range(3)],lambda a,b,c:a+b+c)
        assign(m.col(2),[tmp.col(0,True)],lambda x:scale*x)


@port('ArrayFunction')
def array_function(g):
    with Array(5,3,g,-1) as a:
        for i in range(5):
            for j in range(3):a[i,j]=i-2+10*j
        with Function([a],sum_columns) as f:
            f(a);check_rows(a,[[i-2,i+8,4*(3*i+24)] for i in range(5)])
            with Array(3,3,g,-1) as b:f(b);check_rows(b,[[-1,-1,-12]]*3)


@port('MultiMatrixArrayFunction')
def multi_matrix(g):
    with Array(3,2,g,1) as a,Array(3,3,g,2) as b:
        for i in range(3):
            for j in range(2):a[i,j]=i+10*j
            for j in range(3):b[i,j]=i*2+20*j
            b[i,2]=i*4
        def body(x,y):
            with Array(x.rows,1,x.group) as temp:
                assign(temp.col(0),[x.col(0,True),y.col(2,True)],lambda a,b:a+b)
                copy(x.col(1),temp.col(0,True))
        with Function([a,b],body,readonly=[0,1]) as f:f(a,b,readonly=[0,1])
        assert [a[i,1] for i in range(3)]==[0,5,10]
        assert [a[i,0] for i in range(3)]==[0,1,2] and [b[i,0] for i in range(3)]==[0,2,4]


@port('VectorInMatrixFunction')
def vector_in_matrix(g):
    with Array(5,3,g,-1) as a,Array(5,1,g) as v:
        for i in range(5):
            for j in range(3):a[i,j]=i-2+10*j
            v[i]=(i+1)*100
        def body(m,v):
            with Array(m.rows,1,g) as tmp:
                assign(tmp.col(0),[m.col(j,True) for j in range(3)],lambda x,y,z:x+y+z)
                assign(m.col(2),[v.col(0,True),tmp.col(0,True)],lambda x,y:x*y)
        with Function([a,v],body,readonly=[0,1]) as f:f(a,v)
        check_rows(a,[[i-2,i+8,(i+1)*100*(3*i+24)] for i in range(5)])


@port('MultiMatrixDifferentRowsFromCreation')
def different_rows(g):
    with ExitStack() as stack:
        a=stack.enter_context(Array(3,2,g));b=stack.enter_context(Array(3,3,g))
        def body(x,y):assign(x.col(1),[x.col(0,True),y.col(2,True)],lambda u,v:(u+v)*2)
        with Function([a,b],body,readonly=[0,1]) as f:
            x=stack.enter_context(Array(5,2,g));y=stack.enter_context(Array(5,3,g))
            for i in range(5):x[i,0]=i+1;y[i,2]=i*10
            f(x,y);assert [x[i,1] for i in range(5)]==[2,24,46,68,90]
            x2=stack.enter_context(Array(2,2,g));y2=stack.enter_context(Array(2,3,g))
            for i,(u,v) in enumerate([(10,5),(20,15)]):x2[i,0]=u;y2[i,2]=v
            f(x2,y2);assert [x2[i,1] for i in range(2)]==[30,70]


@port('MatrixVectorDifferentRowsFromCreation')
def different_vector_rows(g):
    with Array(3,3,g) as a,Array(3,1,g) as v:
        def body(m,x):assign(m.col(1),[m.col(0,True),x.col(0,True)],lambda a,b:(a+b)*3)
        with Function([a,v],body,readonly=[0,1]) as f,Array(5,3,g) as b,Array(5,1,g) as w:
            for i in range(5):b[i,0]=i+1;w[i]=i*10
            f(b,w);assert [b[i,1] for i in range(5)]==[3,36,69,102,135]


@port('MismatchedRowsAtInvocation')
def mismatch_at_call(g):
    with Array(3,2,g) as a,Array(3,1,g) as v,Function([a,v],lambda *a:None,readonly=[0,1]) as f:
        with Array(5,2,g) as b,Array(5,1,g) as w:
            f(b,v,status=1);f(a,w,readonly=[0,1],status=1)


@port('MultipleMatricesMismatchedRowsAtInvocation')
def mismatch_matrices(g):
    with Array(3,2,g) as a,Array(3,3,g) as b,Function([a,b],lambda *a:None) as f:
        with Array(5,2,g) as x:f(x,b,status=1)


@port('WrongColumnCountAtInvocation')
def wrong_columns(g):
    with Array(4,3,g) as a,Function([a],lambda m:sum_columns(m,1)) as f:
        for cols in [4,2]:
            with Array(4,cols,g) as x:f(x,status=1)
        with Array(7,3,g) as b:f(b)


@port('MismatchedRowDimensions')
def mismatched_creation(g):
    with Array(3,3,g,1) as a,Array(4,3,g,2) as b:
        def body(x,y):assign(x.col(1),[x.col(0,True),y.col(0,True)],lambda a,b:a+b)
        with Function([a,b],body,readonly=[0,1]) as f:
            with Array(5,3,g) as x,Array(5,3,g,1) as y:f(x,y)
            with Array(3,3,g) as x,Array(4,3,g,1) as y:f(x,y,status=1)


@port('MismatchedColumnDimensions')
def bad_column_access(g):
    with Array(3,4,g,1) as a:
        def body(m):
            view=View();return lib.chem_array_column(m.handle,5,0,C.byref(view))
        with Function([a],body) as f:f(a,status=1)


@port('WrongMatrixDimensions')
def wrong_dimensions(g):
    with Array(3,4,g,1) as a,Array(3,5,g,2) as b:
        def body(m):assign(m.col(3),[m.col(0,True)],lambda x:x*2)
        with Function([a],body) as f:
            f(a);f(a);f(b,status=1)
            with Array(7,4,g,1) as c:f(c)


@port('MultipleTemporaries')
def temporaries(g):
    with Array(4,5,g) as a:
        for i in range(4):a[i,0]=i+1;a[i,1]=(i+1)*10
        def body(m):
            with Array(m.rows,1,g) as p,Array(m.rows,1,g) as s:
                assign(p.col(0),[m.col(0,True),m.col(1,True)],lambda x,y:x*y)
                assign(s.col(0),[m.col(0,True),m.col(1,True)],lambda x,y:x+y)
                for j,fn in [(2,lambda x,y:x+y),(3,lambda x,y:x-y),(4,lambda x,y:x*y)]:assign(m.col(j),[p.col(0,True),s.col(0,True)],fn)
        with Function([a],body) as f:f(a)
        for i,want in [(0,[21,-1,110]),(1,[62,18,880]),(3,[204,116,7040])]:assert a.read_row(i)[2:]==want


@port('ColumnViewReuse')
def view_reuse(g):
    with Array(3,4,g) as a:
        for i in range(3):a[i,0]=i+1
        def body(m):
            a,b,c,d=[m.col(j) for j in range(4)]
            assign(b,[a],lambda x:x*2);assign(c,[a,b],lambda x,y:x+y);assign(d,[c,b],lambda x,y:x*y)
        with Function([a],body) as f:f(a)
        for i,want in enumerate([[2,3,6],[4,6,24],[6,9,54]]):assert a.read_row(i)[1:]==want


@port('FunctionReusability')
def reuse_function(g):
    with Array(2,3,g,1) as a:
        with Function([a],lambda m:sum_columns(m,2)) as f:
            for i in range(2):
                for j in range(3):a[i,j]=i+j
            f(a);assert [a[i,2] for i in range(2)]==[6,12]
            with Array(2,3,g,5) as b:f(b);assert [b[i,2] for i in range(2)]==[30,30]
            with Array(2,3,g) as c:
                c[1,0]=10;f(c);assert [c[i,2] for i in range(2)]==[0,20]


@port('ConstMatrixFunction')
def const_function(g):
    with Array(3,4,g) as a:
        for i in range(3):
            for j in range(4):a[i,j]=i*10+j
        def body(m):
            with Array(m.rows,1,g) as tmp:assign(tmp.col(0),[m.col(0,True),m.col(1,True)],lambda x,y:x+y)
            reduce([m.col(2,True)],lambda acc,x:acc+x)
        with Function([a],body,readonly=[1]) as f:f(a,readonly=[1])
        assert a[0,0]==0 and a[1,2]==12


@port('EmptyMatrixFunction')
def empty_function(g):
    with Array(0,3,g,1) as a:
        def never(p):raise AssertionError('empty row callback executed')
        with Function([a],lambda m:apply([m.col(0)],never)) as f:f(a)
    with Array(3,0,g,1) as a,Function([a],lambda m:None) as f:f(a)


def reject_mismatched(g,shapes):
    with ExitStack() as stack:
        arrays=[stack.enter_context(Array(r,c,g,1)) for r,c in shapes]
        # A mismatch must fail before this function can be invoked.
        def forbidden(*a):raise AssertionError('invalid argument shapes reached body')
        with Function(arrays,forbidden) as f:f(*arrays,status=1)


@port('VectorTooSmall')
def vector_small(g):reject_mismatched(g,[(5,3),(3,1)])
@port('VectorTooLarge')
def vector_large(g):reject_mismatched(g,[(5,3),(10,1)])
@port('EmptyVectorNonEmptyMatrix')
def vector_empty(g):reject_mismatched(g,[(5,3),(0,1)])
@port('NonEmptyVectorEmptyMatrix')
def matrix_empty(g):reject_mismatched(g,[(0,3),(5,1)])
@port('MultipleVectorsDifferentSizes')
def vectors_mismatch(g):reject_mismatched(g,[(5,3),(5,1),(3,1)])
@port('MultipleMatricesDifferentRowsVector')
def matrices_vector_mismatch(g):reject_mismatched(g,[(4,2),(5,3),(4,1)])
@port('VectorSizeMatchesOneMatrixOnly')
def vector_one_only(g):reject_mismatched(g,[(5,2),(5,3),(4,1)])


@port('EmptyVectorEmptyMatrix')
def all_empty(g):
    with Array(0,3,g,1) as a,Array(0,1,g) as v:
        def body(m,v):assign(m.col(0),[v.col(0,True)],lambda x:x)
        with Function([a,v],body,readonly=[0,1]) as f:f(a,v)


@port('MultipleVectorsSameSize')
def two_vectors(g):
    with Array(5,3,g) as a,Array(5,1,g) as u,Array(5,1,g) as v:
        for i in range(5):u[i]=i+1;v[i]=(i+1)*10
        def body(m,x,y):assign(m.col(0),[x.col(0,True),y.col(0,True)],lambda a,b:a+b)
        with Function([a,u,v],body,readonly=[0,1,1]) as f:f(a,u,v)
        assert [a[i,0] for i in range(5)]==[11,22,33,44,55]


@port('MultipleMatricesOneVector')
def matrices_one_vector(g):
    with Array(4,2,g) as a,Array(4,3,g) as b,Array(4,1,g) as v:
        for i in range(4):v[i]=i*2;a[i,0]=i+1;b[i,0]=i*10
        def body(x,y,v):
            assign(x.col(1),[x.col(0,True),v.col(0,True)],lambda a,b:a+b)
            assign(y.col(1),[y.col(0,True),v.col(0,True)],lambda a,b:a+b)
        with Function([a,b,v],body,readonly=[0,0,1]) as f:f(a,b,v)
        assert [a[i,1] for i in range(4)]==[1,4,7,10]
        assert [b[i,1] for i in range(4)]==[0,12,24,36]


@port('ConstVector')
def const_vector(g):
    with Array(3,2,g) as a,Array.from_rows([[10],[20],[30]],g) as v:
        with Function([a,v],lambda m,v:assign(m.col(0),[v.col(0,True)],lambda x:x*2),readonly=[0,1]) as f:f(a,v,readonly=[0,1])
        assert [a[i,0] for i in range(3)]==[20,40,60]


@port('MutableVector')
def mutable_vector(g):
    with Array(3,2,g) as a,Array.from_rows([[5],[10],[15]],g) as v:
        for i in range(3):a[i,0]=i+1
        with Function([a,v],lambda m,v:assign(v.col(0),[m.col(0,True)],lambda x:x*3),readonly=[1,0]) as f:f(a,v,readonly=[1,0])
        assert [v[i] for i in range(3)]==[3,6,9]


@port('FunctionReusabilityWithVectors')
def reusable_vectors(g):
    with Array(3,2,g) as a,Array.from_rows([[1],[2],[3]],g) as v:
        with Function([a,v],lambda m,v:assign(m.col(0),[v.col(0,True)],lambda x:x*10),readonly=[0,1]) as f:
            f(a,v);assert [a[i,0] for i in range(3)]==[10,20,30]
            with Array.from_rows([[5],[6],[7]],g) as w:f(a,w);assert [a[i,0] for i in range(3)]==[50,60,70]


@port('FunctionInvocationWithWrongSizedVector')
def wrong_vector_at_call(g):
    with Array(3,2,g) as a,Array(3,1,g) as v,Array(5,1,g) as w:
        with Function([a,v],lambda m,v:copy(m.col(0),v.col(0,True)),readonly=[0,1]) as f:
            f(a,v);f(a,w,status=1) # The C API also checks the commented-out reference assertion.


@port('MixedVectorColumnViewRowVariable')
def mixed_vector(g):
    with Array(4,3,g) as a,Array(4,1,g) as v:
        for i in range(4):a[i,0]=i+1;a[i,1]=(i+1)*10;v[i]=(i+1)*100
        def body(m,v):
            with Array(m.rows,1,g) as tmp:
                assign(tmp.col(0),[m.col(0,True),m.col(1,True)],lambda x,y:x+y)
                assign(m.col(2),[tmp.col(0,True),v.col(0,True)],lambda x,y:x+y)
        with Function([a,v],body,readonly=[0,1]) as f:f(a,v)
        assert [a[i,2] for i in range(4)]==[111,222,333,444]


@port('IntegerVector')
def integer_vector(g):
    with Array(3,2,g,type=1) as a,Array.from_rows([[10],[20],[30]],g,type=1) as v:
        with Function([a,v],lambda m,v:assign(m.col(0),[v.col(0,True)],lambda x:x*2),readonly=[0,1]) as f:f(a,v)
        assert [a[i,0] for i in range(3)]==[20,40,60]


@port('FunctionWithConstSignature')
def const_signature(g):
    with Array(3,2,g) as a,Array(3,1,g) as v:
        with Function([a,v],lambda m,v:assign(m.col(0),[v.col(0,True)],lambda x:x*2),readonly=[0,1]) as f:f(a,v,readonly=[0,1])


@port('TestFill')
def fills(g):
    with Array(3,2,g) as a:
        with Function([a],lambda m:fill(m.col(1),3.2)) as f:f(a)
        check_rows(a,[[0,3.2]]*3)
        with Array(3,1,g) as v:
            with Function([a,v],lambda m,v:fill(v.col(0),3.2),readonly=[1,0]) as f:f(a,v)
            assert [v[i] for i in range(3)]==[3.2]*3
        def body(m):
            with Array(m.rows,1,g) as tmp:
                fill(tmp.col(0),9.9);assign(m.col(0),[tmp.col(0,True)],lambda x:x)
        with Function([a],body) as f:f(a)
        assert [a[i,0] for i in range(3)]==[9.9]*3


@port('TestCopy')
def copies(g):
    with Array(3,2,g) as a,Array.from_rows([[3.2],[4.2],[1.3]],g) as v:
        with Function([a,v],lambda m,v:copy(m.col(1),v.col(0,True)),readonly=[0,1]) as f:f(a,v)
        check_rows(a,[[0,3.2],[0,4.2],[0,1.3]])
        with Array(3,1,g) as v2:
            with Function([a,v2],lambda m,v:copy(v.col(0),m.col(1,True)),readonly=[1,0]) as f:f(a,v2)
            assert [v2[i] for i in range(3)]==[3.2,4.2,1.3]
        with Function([a],lambda m:copy(m.col(0),m.col(1))) as f:f(a)
        check_rows(a,[[3.2,3.2],[4.2,4.2],[1.3,1.3]])
        for i in range(3):a[i,0]=0
        with Function([a],lambda m:copy(m.col(0),m.col(1,True))) as f:f(a)
        assert [a[i,0] for i in range(3)]==[3.2,4.2,1.3]
        def body(m):
            with Array(m.rows,1,g) as tmp:
                copy(tmp.col(0),m.col(1,True));fill(m.col(1),0);assign(m.col(1),[tmp.col(0,True)],lambda x:x)
        with Function([a],body) as f:f(a)
        assert [a[i,1] for i in range(3)]==[3.2,4.2,1.3]


@port('ReduceSum')
def sum_reduce(g):
    with Array.from_rows([[1,2],[3,4],[5,6]],g) as a:
        result=[0]
        def body(m):
            for j in range(2):result[0]=reduce([m.col(j,True)],lambda acc,x:acc+x*x,initial=result[0])
        with Function([a],body,readonly=[1]) as f:f(a)
        assert result[0]==91


@port('ReduceMax')
def max_reduce(g):
    with Array.from_rows([[1,8],[5,2],[3,4]],g) as a:
        result=[-sys.float_info.max]
        def body(m):
            for j in range(2):result[0]=reduce([m.col(j,True)],lambda acc,x:max(acc,x),initial=result[0])
        with Function([a],body,readonly=[1]) as f:f(a)
        assert result[0]==8


@port('ReduceLOr')
def or_reduce(g):
    for exceed in [False,True]:
        with Array(3,2,g,1) as a:
            if exceed:a[2,1]=42
            result=[False]
            def body(m):
                for j in range(2):result[0]=reduce([m.col(j,True)],lambda acc,x:acc or x>10,initial=result[0],type=2)
            with Function([a],body,readonly=[1]) as f:f(a)
            assert result[0]==exceed


@port('ReduceLAnd')
def and_reduce(g):
    for nonfinite in [False,True]:
        with Array(3,2,g,1) as a:
            if nonfinite:a[1,0]=math.nan
            result=[True]
            def body(m):
                for j in range(2):result[0]=reduce([m.col(j,True)],lambda acc,x:acc and math.isfinite(x),initial=result[0],type=2)
            with Function([a],body,readonly=[1]) as f:f(a)
            assert result[0]==(not nonfinite)


@port('ReduceStrict')
def strict_reduce(g):
    with Array(3,1,g) as a:
        result=[]
        with Function([a],lambda m:result.append(reduce([m.col(0,True)],lambda acc,x:acc+1,strict=True)),readonly=[1]) as f:f(a)
        assert result==[3]


@port('VectorCapture')
def vector_capture(g):
    with Array(3,3,g) as a,Array.from_rows([[2],[0],[1]],g,type=1) as v1,Array.from_rows([[10],[20],[30]],g,type=1) as v2:
        view1=v1.col(0);view2=v2.col(0,True)
        def body(m):
            for i in range(view1.length):
                value=element(view2,element(view1,i)[0])[0]
                apply([m.col(0)],lambda p:p[0].__setitem__(0,p[0][0]+value))
            value=element(view2,element(view1,1)[0])[0]
            apply([m.col(1)],lambda p:p[0].__setitem__(0,p[0][0]+value))
        with Function([a],body) as f:f(a)
        check_rows(a,[[60,10,0]]*3)
