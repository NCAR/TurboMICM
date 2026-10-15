"""Sparse matrix policy assertions for CSR, CSC, and grouped versions.

Original assertions: Copyright (C) 2023-2026 University Corporation for Atmospheric Research.
SPDX-License-Identifier: Apache-2.0
"""
from arrays import *
from harness import case, equal
from matrix_cases import Pool, map_to

POLICIES = {}


def policy(*names):
    def register(function):
        for name in names: POLICIES[name] = function
        return function
    return register


class SparsePool(Pool):
    def __init__(self,group,column_order): super().__init__(group); self.column_order = column_order
    def sparse(self,blocks,order,pattern,initial=0.0,group=None):
        a = Array(blocks,group=group or self.group,initial=initial,pattern=pattern,order=order,column_order=self.column_order)
        self.resources.append(a); return a


PATTERN = [(0,1),(3,2),(2,3),(1,1),(2,1)]


@policy('ZeroMatrix','ConstZeroMatrix','SetScalar','SingleBlockMatrix','ConstSingleBlockMatrix','MultiBlockMatrix','AddToDiagonal','BadConfiguration')
def construction(p):
    zero = p.sparse(1,3,[]); equal(zero.data.columns,0)
    equal(lib.chem_sparse_index(zero.p,0,1,1),MAX); equal(lib.chem_sparse_index(zero.p,0,6,0),MAX)
    zero.fill(2); equal(zero.data.allocated,0)
    a = p.sparse(1,4,PATTERN); equal(a.data.columns,5)
    for i,j in PATTERN: assert a.data.index[i*4+j] >= 0
    for i,j in [(2,2),(0,0),(3,3)]: equal(a.data.index[i*4+j],-1)
    equal(a[0,2,1],0); a[0,2,1] = 45; a[0,3,2] = 42; a[0,2,3] = 21
    equal(a[0,2,1],45)
    order = sorted(set(PATTERN),key=lambda rc:(rc[1],rc[0]) if p.column_order else rc)
    for i,(r,c) in enumerate(order): equal(lib.chem_sparse_index(a.p,0,r,c),i*p.group)
    equal(a.data.values[lib.chem_sparse_index(a.p,0,3,2)],42)
    equal(a.data.values[lib.chem_sparse_index(a.p,0,2,3)],21)
    a = p.sparse(53,4,PATTERN,24)
    for block,r,c,v in [(0,2,1,45),(43,3,2,29),(33,2,3,79),(2,2,3,64)]: a[block,r,c] = v; equal(a[block,r,c],v)
    for block in (0,1,2):
        j = order.index((1,1)); expected = (block//p.group*5+j)*p.group+block%p.group
        equal(lib.chem_sparse_index(a.p,block,1,1),expected)
    a = p.sparse(3,4,PATTERN+[(0,1)],24); equal(a.data.columns,5)
    for b,v in enumerate([1,3,5]): a[b,1,1] = v
    lib.chem_array_add_diagonal(a.p,7)
    for b,v in enumerate([8,10,12]):
        equal(a[b,1,1],v)
        for r,c in [(0,1),(2,1),(2,3),(3,2)]: equal(a[b,r,c],24)
    for r,c in [(3,0),(2,4),(6,7)]:
        pointer = lib.chem_array_sparse(1,3,p.group,p.column_order,inds([r]),inds([c]),1,0)
        assert not pointer


@policy('Print','PrintNonZero')
def printing(p):
    a = p.sparse(3,4,PATTERN,24)
    for b,v in enumerate([1,3,5]): a[b,1,1] = v
    lib.chem_array_add_diagonal(a.p,7)
    expected = ''.join(f'Block {b}\n0,24,0,0\n0,{v},0,0\n0,24,0,24\n0,0,24,0\n' for b,v in enumerate([8,10,12]))
    equal(a.text(),expected)
    a = p.sparse(3,4,PATTERN,32)
    for b,v in enumerate([2,4,6]): a[b,1,1] = v
    lib.chem_array_add_diagonal(a.p,5)
    expected = ''.join(f'Block {b}\n0, 1, 32\n1, 1, {v}\n2, 1, 32\n2, 3, 32\n3, 2, 32\n' for b,v in enumerate([7,9,11]))
    equal(a.text(True),expected)


@policy('ArrayFunction')
def array_function(p):
    pattern = [(0,0),(1,1),(2,2),(2,3),(3,3)]
    a = p.sparse(3,4,pattern,1)
    a[0,0,0] = 1; a[1,1,1] = 2; a[2,2,2] = 3; a[2,2,3] = 4; a[2,3,3] = 5
    def operation(a):
        tmp = Vector([0]*a.data.rows)
        map_to(tmp.view(),[a.block(*pair) for pair in pattern[:4]],lambda a,b,c,d:a+b+c+d)
        map_to(a.block(2,3),[tmp.view()],lambda x:2*x)
    f = p.function([a],operation); equal(f.invoke(a),0)
    for b,v in enumerate([8,10,18]): equal(a[b,2,3],v)
    for b in range(3):
        equal(a[b,0,0],1); equal(a[b,1,1],2 if b == 1 else 1); equal(a[b,2,2],3 if b == 2 else 1); equal(a[b,3,3],5 if b == 2 else 1)
    second = p.sparse(3,4,pattern,-1); equal(f.invoke(second),0)
    for b in range(3):
        equal(second[b,2,3],-8)
        for pair in pattern:
            if pair != (2,3): equal(second[(b,)+pair],-1)


@policy('MultiMatrixArrayFunction','TwoSparseMatricesDifferentStructure','SparseAndDenseMatrixFunction','SparseAndVectorMatrixFunction')
def multi(p):
    a,b = p.sparse(3,3,[(0,1),(1,2)]),p.sparse(3,3,[(0,0),(1,1),(2,2)])
    for i in range(3):
        a[i,0,1] = i+10; a[i,1,2] = i*2+20
        for j,factor in enumerate([4,3,5]): b[i,j,j] = i*factor
    f = p.function([a,b],lambda a,b:map_to(a.block(1,2),[a.block(0,1),b.block(2,2)],lambda x,y:x+y),[False,True]); equal(f.invoke(a,b),0)
    for i in range(3): equal(a[i,1,2],10+6*i); equal(a[i,0,1],10+i); equal(b[i,0,0],4*i)
    a,b = p.sparse(4,3,[(0,1),(1,1),(2,2)]),p.sparse(4,5,[(0,0),(1,2),(2,3),(3,4)])
    for i in range(4):
        for pair,value_ in [((0,1),i*2+1),((1,1),i*3+2),((2,2),i*5+3)]: a[(i,)+pair] = value_
        for j,pair in enumerate([(0,0),(1,2),(2,3),(3,4)]): b[(i,)+pair] = i+(j+1)*10
    f = p.function([a,b],lambda a,b:map_to(a.block(2,2),[a.block(0,1),a.block(1,1),b.block(0,0)],lambda x,y,z:x+y+z),[False,True]); equal(f.invoke(a,b),0)
    for i in range(4): equal(a[i,2,2],13+6*i); equal(a[i,0,1],2*i+1)
    for i,j in enumerate([(0,0),(1,2),(2,3),(3,4)]): equal(b[(i,)+j],i+(i+1)*10)
    a,b = p.sparse(3,4,[(0,1),(1,2),(2,3),(3,3)]),p.matrix(3,4)
    for i in range(3):
        for pair,v in zip([(0,1),(1,2),(2,3),(3,3)],[i+1,i+10,i+20,i+30]): a[(i,)+pair] = v
        for j in range(4): b[i,j] = i*100+j*10
    f = p.function([a,b],lambda a,b:map_to(a.block(3,3),[a.block(0,1),a.block(1,2),b.column(0),b.column(1)],lambda a,b,c,d:a+b+c+d),[False,True])
    equal(f.invoke(a,b),0)
    for i,v in enumerate([21,223,425]): equal(a[i,3,3],v)
    equal(a[0,0,1],1); equal(a[1,1,2],11); equal(a[2,2,3],22); equal(b[0,0],0); equal(b[1,2],120); equal(b[2,3],230)


@policy('MultipleTemporaries','BlockViewReuse','FunctionReusability')
def temporaries_reuse(p):
    pattern = [(0,1),(1,2),(2,3),(3,4),(4,4)]
    a = p.sparse(4,5,pattern)
    for i in range(4): a[i,0,1] = i+1; a[i,1,2] = 10*(i+1)
    def operation(a):
        one,two = Vector([0]*a.data.rows),Vector([0]*a.data.rows)
        map_to(one.view(),[a.block(0,1),a.block(1,2)],lambda x,y:x*y)
        map_to(two.view(),[a.block(0,1),a.block(1,2)],lambda x,y:x+y)
        for pair,fn in [((2,3),lambda x,y:x+y),((3,4),lambda x,y:x-y),((4,4),lambda x,y:x*y)]: map_to(a.block(*pair),[one.view(),two.view()],fn)
    f = p.function([a],operation); equal(f.invoke(a),0)
    for i in range(4):
        product,total = 10*(i+1)**2,11*(i+1)
        for pair,v in [((2,3),product+total),((3,4),product-total),((4,4),product*total)]: equal(a[(i,)+pair],v)
    a = p.sparse(3,4,[(0,1),(1,2),(2,3),(3,3)])
    for i in range(3): a[i,0,1] = i+1
    def operation(a):
        c0,c1,c2,c3 = [a.block(*pair) for pair in [(0,1),(1,2),(2,3),(3,3)]]
        map_to(c1,[c0],lambda x:2*x); map_to(c2,[c0,c1],lambda x,y:x+y); map_to(c3,[c2,c1],lambda x,y:x*y)
    f = p.function([a],operation); equal(f.invoke(a),0)
    for i in range(3): equal(a[i,1,2],2*(i+1)); equal(a[i,2,3],3*(i+1)); equal(a[i,3,3],6*(i+1)**2)
    pattern = [(0,0),(1,1),(2,2)]; a,b,c = [p.sparse(2,3,pattern,v) for v in [1,5,0]]
    for i in range(2):
        for j in range(3): a[i,j,j] = i+j
        c[i,0,0] = i*10
    def operation(a): map_to(a.block(2,2),[a.block(0,0),a.block(1,1),a.block(2,2)],lambda x,y,z:2*(x+y+z))
    f = p.function([a],operation)
    for obj,expected in [(a,[6,12]),(b,[30,30]),(c,[0,20])]:
        equal(f.invoke(obj),0)
        for i,v in enumerate(expected): equal(obj[i,2,2],v)


@policy('MismatchedBlockDimensions','MismatchedElementDimensions','MismatchedBlocksAtInvocation',
        'MultipleSparseMatricesMismatchedBlocksAtInvocation','VectorTooSmall','VectorTooLarge',
        'EmptyVectorNonEmptySparseMatrix','NonEmptyVectorEmptySparseMatrix','EmptyVectorEmptySparseMatrix',
        'MultipleVectorsDifferentSizes','MultipleSparseMatricesDifferentBlocksVector','VectorSizeMatchesOneSparseMatrixOnly',
        'FunctionInvocationWithWrongSizedVectorSparse')
def dimensions(p):
    pattern = [(0,1),(1,1)]
    a,b = p.sparse(3,3,pattern),p.sparse(4,3,pattern)
    f = p.function([a,a],lambda *objects:None); equal(f.invoke(a,a),0); equal(f.invoke(a,b),4)
    a = p.sparse(3,4,[(0,1),(1,1),(2,2)])
    equal(lib.chem_sparse_index(a.p,0,3,3),MAX)
    for rows,nvec,expected in [(3,2,4),(3,5,4),(3,0,4),(0,5,4),(0,0,0),(5,3,4),(5,10,4)]:
        a,v = p.sparse(rows,3,[(0,1)]),Vector([0]*nvec)
        f = p.function([a,v],lambda *objects:None); equal(f.invoke(a,v),expected)
    a,b,v,w = p.sparse(3,3,pattern),p.sparse(4,3,pattern),Vector([0]*3),Vector([0]*2)
    f = p.function([a,v,w],lambda *objects:None); equal(f.invoke(a,v,w),4)
    f = p.function([a,b,v],lambda *objects:None); equal(f.invoke(a,b,v),4)
    f = p.function([a,v],lambda *objects:None); equal(f.invoke(a,v),0); equal(f.invoke(a,Vector([0]*5)),4)


@policy('MultipleSparseMatricesDifferentBlocksFromCreation','SparseMatrixVectorDifferentBlocksFromCreation')
def flexible(p):
    pattern = [(0,1),(1,2),(2,3)]
    a,b = p.sparse(3,4,pattern),p.sparse(3,4,pattern)
    def operation(a,b): map_to(a.block(2,3),[a.block(0,1),a.block(1,2),b.block(2,3)],lambda x,y,z:x+y+z)
    f = p.function([a,b],operation,[False,True]); x,y = p.sparse(4,4,pattern),p.sparse(4,4,pattern)
    for i in range(4): x[i,0,1] = i+1; x[i,1,2] = i+10; y[i,2,3] = i+100
    equal(f.invoke(x,y),0)
    for i in range(4): equal(x[i,2,3],111+3*i)
    a,v = p.sparse(3,3,[(0,1),(1,2)]),Vector([1,2,3])
    f = p.function([a,v],lambda a,v:map_to(a.block(1,2),[a.block(0,1),v.view()],lambda x,y:x+y),[False,True])
    x,w = p.sparse(5,3,[(0,1),(1,2)]),Vector([10,20,30,40,50])
    for i in range(5): x[i,0,1] = i+1
    equal(f.invoke(x,w),0)
    for i in range(5): equal(x[i,1,2],11*(i+1))


@policy('VectorInSparseMatrixFunction','MultipleVectorsSameSize','MultipleSparseMatricesOneVector','ConstVectorSparse',
        'MutableVectorSparse','FunctionWithConstSignatureSparse','FunctionReusabilityWithVectorsSparse',
        'MixedVectorBlockViewBlockVariable','IntegerVectorSparse')
def vectors(p):
    a = p.sparse(3,4,[(0,1),(1,2),(2,3)])
    for i in range(3):
        for j,pair in enumerate([(0,1),(1,2),(2,3)]): a[(i,)+pair] = i+j+1
    for factor,values in [(2,[5,10,15]),(3,[5,10,15]),(1.5,[5,10,15]),(2,[1,2,3]),(2,[10,20,30])]:
        v = Vector(values); f = p.function([a,v],lambda a,v:map_to(a.block(0,1),[v.view()],lambda x:factor*x),[False,True])
        equal(f.invoke(a,v),0)
        for i,vv in enumerate(values): equal(a[i,0,1],vv*factor)
    v,w = Vector([1,2,3]),Vector([4,5,6])
    f = p.function([a,v,w],lambda a,v,w:map_to(a.block(0,1),[v.view(),w.view()],lambda x,y:x+y),[False,True,True])
    equal(f.invoke(a,v,w),0)
    for i,val in enumerate([5,7,9]): equal(a[i,0,1],val)
    b,v = p.sparse(3,3,[(0,1)]),Vector([10,20,30])
    for i in range(3): a[i,0,1] = i+1
    f = p.function([a,b,v],lambda a,b,v:map_to(b.block(0,1),[a.block(0,1),v.view()],lambda x,y:x+y),[True,False,True])
    equal(f.invoke(a,b,v),0)
    for i in range(3): equal(b[i,0,1],11*(i+1))
    f = p.function([a,v],lambda a,v:map_to(v.view(),[a.block(0,1)],lambda x:3*x),[True,False]); equal(f.invoke(a,v),0)
    equal(list(v.values[:3]),[3,6,9])
    v = Vector([100,200,300])
    for i in range(3): a[i,0,1] = i+1; a[i,1,2] = i+2
    def mixed(a,v):
        tmp = Vector([0]*a.data.rows)
        map_to(tmp.view(),[a.block(0,1),a.block(1,2),v.view()],lambda x,y,z:x+y+z)
        map_to(a.block(2,3),[tmp.view()],lambda x:2*x)
    f = p.function([a,v],mixed,[False,True]); equal(f.invoke(a,v),0)
    for i,val in enumerate([206,410,614]): equal(a[i,2,3],val)


@policy('GetBlockViewByVectorIndex')
def vector_indices(p):
    pattern = [(0,0),(0,2),(1,0),(1,1),(2,1),(2,2)]
    a,v = p.sparse(3,3,pattern),Vector([0,1,2])
    raw = [lib.chem_sparse_index(a.p,0,*pair) for pair in pattern]
    def operation(a,v):
        for j,offset in enumerate(raw):
            map_to(a.column(offset//p.group),[v.view()],lambda x:x+j*10)
    f = p.function([a,v],operation,[False,True]); equal(f.invoke(a,v),0)
    for i in range(3):
        for j,pair in enumerate(pattern): equal(a[(i,)+pair],i+j*10)


@policy('ConstSparseMatrixFunction','EmptySparseMatrixFunction')
def const_empty(p):
    a = p.sparse(3,4,[(0,1),(1,2),(2,3),(3,3)])
    for i in range(3): a[i,0,1] = i+1; a[i,1,2] = i+10
    def read(a):
        t = Vector([0]*a.data.rows)
        map_to(t.view(),[a.block(0,1),a.block(1,2)],lambda x,y:x+y)
        foreach([a.block(2,3,True)],lambda x:None)
    f = p.function([a],read,[True]); equal(f.invoke(a),0); equal(a[0,0,1],1); equal(a[1,1,2],11)
    for a in [p.sparse(3,4,[]),p.sparse(0,4,[(0,1)])]:
        f = p.function([a],lambda a:None,[True]); equal(f.invoke(a),0)


@policy('TestFill','TestCopy')
def fill_copy(p):
    a = p.sparse(3,3,[(0,1),(2,0)]); col = a.block(0,1)
    equal(lib.chem_view_fill(C.byref(col),3.2),0)
    for i in range(3): equal(a[i,0,1],3.2); equal(a[i,2,0],0)
    v = Vector([0]*3); vv = v.view(); equal(lib.chem_view_fill(C.byref(vv),3.2),0); equal(list(v.values[:3]),[3.2]*3)
    tmp = Vector([0]*3); tv = tmp.view(); equal(lib.chem_view_fill(C.byref(tv),9.9),0)
    map_to(a.block(2,0),[tv],lambda x:x)
    for i in range(3): equal(a[i,2,0],9.9)
    v = Vector([10,20,30]); vv = v.view(True); col = a.block(0,1)
    equal(lib.chem_view_copy(C.byref(col),C.byref(vv)),0)
    for i,vv in enumerate([10,20,30]): equal(a[i,0,1],vv)
    for i,val in enumerate([3.2,4.2,1.3]): a[i,0,1] = val
    v = Vector([-1]*3); vv = v.view(); col = a.block(0,1,True)
    equal(lib.chem_view_copy(C.byref(vv),C.byref(col)),0); equal(list(v.values[:3]),[3.2,4.2,1.3])
    for ro in (False,True):
        to,source = a.block(2,0),a.block(0,1,ro); equal(lib.chem_view_copy(C.byref(to),C.byref(source)),0)
        for i,val in enumerate([3.2,4.2,1.3]): equal(a[i,2,0],val); equal(a[i,0,1],val)
    temp = Vector([0]*3); tv,source = temp.view(),a.block(0,1,True)
    equal(lib.chem_view_copy(C.byref(tv),C.byref(source)),0)
    target = a.block(0,1); equal(lib.chem_view_fill(C.byref(target),0),0); map_to(target,[tv],lambda x:x)
    for i,val in enumerate([3.2,4.2,1.3]): equal(a[i,0,1],val)
    dense = p.matrix(3,2)
    for i in range(3): dense[i,1] = 100*(i+1)
    map_to(a.block(2,0),[dense.column(1)],lambda x:x)
    for i in range(3): equal(a[i,2,0],100*(i+1))
    vv,source = temp.view(),a.block(0,1,True); equal(lib.chem_view_copy(C.byref(vv),C.byref(source)),0)
    map_to(dense.column(0),[temp.view()],lambda x:x)
    b = p.sparse(3,3,[(1,2)])
    to,source = b.block(1,2),a.block(0,1,True); equal(lib.chem_view_copy(C.byref(to),C.byref(source)),0)
    for i,val in enumerate([3.2,4.2,1.3]): equal(dense[i,0],val); equal(b[i,1,2],val)


@policy('IncompatibleOrdering','IncompatibleVectorOrdering','IncompatibleSparseOrdering')
def incompatible(p):
    a = p.sparse(3,4,[(0,1),(1,2)])
    for group in (1,2,4,8):
        if group == p.group: continue
        dense = Array(3,4,group)
        sparse = p.sparse(3,4,[(0,1),(1,2)],group=group)
        p.resources.append(dense)
        for other in (dense,sparse):
            f = Function([a,other],lambda *args:None); p.resources.append(f); equal(f.status,1)


def run_sparse(identifier):
    family,fullname = identifier.split(':'); name = fullname.split('.')[-1]
    column_order = 'column_ordering' in family
    groups = (1,) if '_standard_' in family else (1,2,3,4,8)
    for group in groups:
        with SparsePool(group,column_order) as pool: POLICIES[name](pool)


for _line in (ROOT/'case-list.txt').read_text().splitlines():
    if _line.startswith('sparse_matrix_'):
        _id = _line.split()[0]; _name = _id.split('.')[-1]
        if _name in POLICIES: case(_id)(run_sparse)
