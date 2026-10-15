"""Runtime grouped-view case ports; compile-time-only cases are documented N/A."""
# Copyright (C) 2023-2026 University Corporation for Atmospheric Research
# SPDX-License-Identifier: Apache-2.0
from sparseffi import *

CASES={}
def case(name):
    def decorate(fn):CASES['GroupedView.'+name]=fn;return fn
    return decorate


def group(s,number,row,column,ro=False):return s.values.group_col(number,s.flat(row,column),ro)


@case('MatrixColumnViewRoundTrip')
def dense_row():
    with Array(3,4) as a:
        for r in range(3):
            for c in range(4):a[r,c]=10*r+c
        x,y,z=a.group_col(1,2),a.group_col(1,2,True),a.group_col(1,3,True)
        assert element(x,0)[0]==12 and element(y,0)[0]==12 and element(z,0)[0]==13
        assign(x,[z],lambda v:v*2)
        assert [a[r,2] for r in range(3)]==[2,26,22]


@case('VectorMatrixColumnViewContiguousBlock')
def dense_contiguous():
    with Array(8,3,4) as a:
        for r in range(8):
            for c in range(3):a[r,c]=100+10*r+c
        x,y=a.group_col(1,1),a.group_col(1,2,True)
        assert x.element_stride==y.element_stride==1
        for i in range(4):assert element(x,i)[0]==100+10*(4+i)+1 and element(y,i)[0]==100+10*(4+i)+2
        assign(x,[y],lambda v:v*2)
        for r in range(8):assert a[r,1]==(2*(100+10*r+2) if r>=4 else 100+10*r+1)


def dense_equivalence(length):
    with Array(3*length,4,length) as raw,Array(3*length,4,length) as grouped:
        for r in range(3*length):
            for c in range(4):raw[r,c]=grouped[r,c]=r*10+c
        def body(m):
            for g in range(3):assign(m.group_col(g,3),[m.group_col(g,0,True)],lambda x:x*3+1)
        with Function([grouped],body) as f:f(grouped)
        for r in range(3*length):raw[r,3]=raw[r,0]*3+1
        assert raw.values()==grouped.values()
for length in [1,2,4,8]:CASES[f'GroupedView.DenseGroupedMatchesRawL{length}']=lambda n=length:dense_equivalence(n)


@case('PlainMatrixGroupedRoundTrip')
def plain_grouped():
    with Array(4,3) as raw,Array(4,3) as grouped:
        for r in range(4):
            for c in range(3):raw[r,c]=grouped[r,c]=r+c*5
        def body(m):
            for g in range(4):assign(m.group_col(g,2),[m.group_col(g,0,True),m.group_col(g,1,True)],lambda x,y:x-y)
        with Function([grouped],body) as f:f(grouped)
        for r in range(4):raw[r,2]=raw[r,0]-raw[r,1]
        assert raw.values()==grouped.values()


@case('SparseStandardBlockViewRoundTrip')
def sparse_standard():
    with Sparse(3,2,[(0,1),(1,2),(2,2)]) as s:
        for b,values in enumerate([[1,2,3],[10,20,30]]):
            for p,v in zip([(0,1),(1,2),(2,2)],values):s[b,*p]=v
        x,y=group(s,1,0,1),group(s,1,1,2,True)
        assert element(x,0)[0]==10 and element(y,0)[0]==20
        assign(x,[y],lambda v:v+5)
        assert s[1,0,1]==25 and s[0,0,1]==1


@case('SparseVectorBlockViewContiguousBlock')
def sparse_contiguous():
    with Sparse(3,8,[(0,1),(1,2)],group=4) as s:
        for b in range(8):s[b,0,1]=100+b;s[b,1,2]=200+b
        x,y=group(s,1,0,1),group(s,1,1,2,True)
        assert x.element_stride==y.element_stride==1
        for i in range(4):assert element(x,i)[0]==104+i and element(y,i)[0]==204+i
        apply([x,y],lambda p:p[0].__setitem__(0,p[0][0]+p[1][0]))
        for b in range(8):assert s[b,0,1]==((100+b)+(200+b) if b>=4 else 100+b)


def mixed_equivalence(length):
    with Sparse(3,2*length,[(0,1)],group=length) as s,Array(2*length,3,length) as d:
        for b in range(2*length):
            s[b,0,1]=1+b
            for c in range(3):d[b,c]=10+b+c
        reference=d.values()
        for b in range(2*length):reference[b][0]=reference[b][1]+s[b,0,1]
        def body(sm,dm):
            for g in range(2):assign(dm.group_col(g,0),[dm.group_col(g,1,True),sm.group_col(g,s.flat(0,1),True)],lambda a,b:a+b)
        with Function([s.values,d],body) as f:f(s.values,d)
        assert d.values()==reference
for length in [1,2,4,8]:CASES[f'GroupedView.MixedSparseAndDenseGroupedL{length}']=lambda n=length:mixed_equivalence(n)


@case('SparseVectorRowColOverloadStillFunctional')
def raw_sparse():
    with Sparse(3,8,[(0,1)],group=4) as s:
        for b in range(8):s[b,0,1]=b+1
        assert reduce([group(s,1,0,1,True)],lambda acc,x:acc+x)==26


@case('ConstGroupViewProducesGroupedViews')
def const_group():
    with Sparse(3,4,[(0,1)],group=4) as s,Array(4,2,4) as d:
        for b in range(4):s[b,0,1]=b+1;d[b,1]=100+b
        x,y=group(s,0,0,1,True),d.group_col(0,1,True)
        assert x.read_only==y.read_only==1
        for i in range(4):assert element(x,i)[0]==i+1 and element(y,i)[0]==100+i


def partial_fill_copy(compression,length):
    with Sparse(3,length+1,[(0,1),(2,0)],group=length,compression=compression) as s:
        for b in range(length+1):s[b,0,1]=b+1;s[b,2,0]=100+b
        fill(group(s,0,0,1),42)
        for b in range(length):assert s[b,0,1]==42 and s[b,2,0]==100+b
        assert s[length,0,1]==length+1
        copy(group(s,0,2,0),group(s,0,0,1,True))
        for b in range(length):assert s[b,2,0]==42 and s[b,0,1]==42
        fill(group(s,(length+1)//length,0,1),-7)
        assert s[length,0,1]==-7
for name,c,g in [('CSR_L4',0,4),('CSC_L4',1,4),('CSR_L8',0,8)]:
    CASES['GroupedView.FillAndCopySparseVector'+name]=lambda c=c,g=g:partial_fill_copy(c,g)


@case('FillAndCopySparseStandardCSR')
def csr_fill_copy():
    with Sparse(3,2,[(0,1),(2,0)]) as s:
        for b in range(2):s[b,0,1]=b+1;s[b,2,0]=100+b
        fill(group(s,0,0,1),42)
        assert s[0,0,1]==42 and s[0,2,0]==100 and s[1,0,1]==2
        copy(group(s,0,2,0),group(s,0,0,1,True));assert s[0,2,0]==42


@case('FillAndCopySparseStandardCSC')
def csc_fill_copy():
    with Sparse(3,2,[(0,1),(2,0)],compression=1) as s:
        for b in range(2):s[b,0,1]=b+1;s[b,2,0]=100+b
        fill(group(s,1,0,1),9)
        assert s[1,0,1]==9 and s[0,0,1]==1 and s[1,2,0]==101
        copy(group(s,1,2,0),group(s,1,0,1,True));assert s[1,2,0]==9


@case('FillCopyCrossGroupViewSources')
def cross_group_sources():
    with Sparse(2,4,[(0,1)],group=4) as src,Sparse(2,4,[(0,1)],group=4) as dst,Sparse(2,4,[(0,1)],group=4) as ref:
        for b in range(4):src[b,0,1]=b+1;ref[b,0,1]=999
        copy(group(dst,0,0,1),group(src,0,0,1,True))
        for b in range(4):assert dst[b,0,1]==b+1
        assign(group(ref,0,0,1),[],lambda:55.5);fill(group(dst,0,0,1),55.5)
        for b in range(4):assert dst[b,0,1]==ref[b,0,1]


def runtime_vector():
    values=(D*5)(1,2,3,4,5)
    view=lib.chem_vector_view(values,5,0,1)
    assert reduce([view],lambda acc,x:acc+x,strict=True)==15
CASES['ViewCategorySFINAE.RuntimeVectorUse']=runtime_vector
