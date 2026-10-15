"""Dense matrix policy tests, shared by ordinary and grouped layouts.

Original assertions: Copyright (C) 2023-2026 University Corporation for Atmospheric Research.
SPDX-License-Identifier: Apache-2.0
"""
import math
from pathlib import Path
from arrays import *
from harness import case, equal, near

POLICIES = {}


def policy(*names):
    def register(function):
        for name in names: POLICIES[name] = function
        return function
    return register


class Pool:
    def __init__(self,group): self.group,self.resources = group,[]
    def matrix(self,rows,columns,initial=0.0):
        a = Array(rows,columns,self.group,initial); self.resources.append(a); return a
    def function(self,objects,function,readonly=()):
        f = Function(objects,function,readonly); self.resources.append(f); equal(f.status,0); return f
    def __enter__(self): return self
    def __exit__(self,*args):
        for resource in reversed(self.resources): resource.close()


def map_to(output,inputs,function):
    inputs = [View.from_buffer_copy(v) for v in inputs]
    for v in inputs: v.readonly = 1
    foreach([output]+inputs,lambda old,*values: (function(*values),)+values)


def check_rows(matrix,rows):
    equal(matrix.data.rows,len(rows))
    for i,row in enumerate(rows):
        equal(matrix.data.columns,len(row))
        for j,value_ in enumerate(row): equal(matrix[i,j],value_)


@policy('SmallMatrix','SmallConstMatrix')
def small(p):
    a = p.matrix(3,5); a[1,3] = 64.7; a[0,0] = 41.2; a[2,4] = 102.3
    equal(a[1,3],64.7); equal(a[0,0],41.2); equal(a[2,4],102.3)
    equal(a.data.allocated,((3+p.group-1)//p.group)*p.group*5)
    for i,j,v in [(1,3,64.7),(0,0,41.2),(2,4,102.3)]:
        at = (i//p.group*5+j)*p.group+i%p.group
        equal(lib.chem_array_index(a.p,i,j),at); equal(a.data.values[at],v)


@policy('InitializeMatrix','InitializeConstMatrix')
def initialize(p):
    a = p.matrix(2,3,12.4)
    for i,j in [(0,0),(1,0),(1,2)]: equal(a[i,j],12.4)


@policy('LoopOverMatrix','LoopOverConstMatrix','Strides','IterateOverMatrix')
def loops(p):
    a = p.matrix(3,4)
    for i in range(3):
        for j in range(4): a[i,j] = i*100+j
    check_rows(a,[[i*100+j for j in range(4)] for i in range(3)])
    values = arr([42,43,44,45]); equal(lib.chem_array_row(a.p,1,values,4,1),0)
    equal(a[0,3],3); equal(a[1,0],42); equal(a[1,1],43); equal(a[1,3],45); equal(a[2,0],200)


@policy('ConversionToVector','ConstConversionToVector')
def to_vector(p):
    a = p.matrix(2,3); values = [13.2,31.2,314.2]
    for j,v in enumerate(values): a[1,j] = v
    output = arr([0]*3); equal(lib.chem_array_row(a.p,1,output,3,0),0); equal(list(output),values)


@policy('ConversionFromVector')
def from_vector(p):
    data = [[412.3,32.4,41.3],[5.33,-.3,31.2]]
    rows = [arr(row) for row in data]; result = AP()
    equal(lib.chem_array_from_rows((DP*2)(*rows),inds([3,3]),2,p.group,C.byref(result)),0)
    a = Array(pointer=result); a.owned = True; p.resources.append(a); check_rows(a,data)
    bad = [arr([3]),arr([4,5]),arr([5])]; result = AP()
    equal(lib.chem_array_from_rows((DP*3)(*bad),inds([1,2,1]),3,p.group,C.byref(result)),4)
    empty = AP(); equal(lib.chem_array_from_rows(None,None,0,p.group,C.byref(empty)),0)
    equal(empty.contents.rows,0); lib.chem_array_free(empty)


@policy('AssignmentFromVector')
def assignment(p):
    a = p.matrix(4,3)
    equal(lib.chem_array_row(a.p,2,arr([12.3,15.1,24.3]),3,1),0)
    equal(a[0,0],0); equal(a[3,0],0)
    for j,v in enumerate([12.3,15.1,24.3]): equal(a[2,j],v)
    equal(lib.chem_array_row(a.p,1,arr([14.3,52.3,65.7,16.34]),4,1),4)
    equal(lib.chem_array_row(a.p,1,arr([13.2,52.8]),2,1),4)


@policy('Axpy','ForEach')
def axpy(p):
    a,b,c = p.matrix(4,3),p.matrix(4,3),p.matrix(4,3)
    for i in range(4):
        for j in range(3): a[i,j] = i*10.3+j*100.5; b[i,j] = i*1.7+j*10.2; c[i,j] = i*19.5+j*32.2
    expected = sum(i*10.3+j*100.5+1.39*(i*1.7+j*10.2) for i in range(4) for j in range(3))
    equal(lib.chem_array_axpy(a.p,1.39,b.p),0)
    near(sum(a[i,j] for i in range(4) for j in range(3)),expected,1e-5)
    for i in range(4):
        for j in range(3): a[i,j] = i*10.3+j*100.5
    for j in range(3): map_to(a.column(j),[a.column(j),b.column(j)],lambda x,y:x+y)
    for i in range(4):
        for j in range(3): near(a[i,j],i*10.3+j*100.5+i*1.7+j*10.2,1e-5)
    for i in range(4):
        for j in range(3): a[i,j] = i*10.3+j*100.5
    for j in range(3): map_to(a.column(j),[a.column(j),b.column(j),c.column(j)],lambda x,y,z:x+y-z)
    for i in range(4):
        for j in range(3): near(a[i,j],(i*10.3+j*100.5)+(i*1.7+j*10.2)-(i*19.5+j*32.2),1e-5)


@policy('SetScaler')
def scalar(p):
    a = p.matrix(2,3); a.fill(2); equal(list(a.data.values[:a.data.allocated]),[2]*a.data.allocated)


@policy('Max','Min')
def bounds(p):
    a = p.matrix(2,3)
    lib.chem_array_bound(a.p,2,1); check_rows(a,[[2]*3]*2)
    a.fill(1); a[1,1] = 3; lib.chem_array_bound(a.p,2,1)
    check_rows(a,[[2,2,2],[2,3,2]])
    a.fill(0); lib.chem_array_bound(a.p,2,0); check_rows(a,[[0]*3]*2)
    a.fill(3); a[1,1] = 1; lib.chem_array_bound(a.p,2,0)
    check_rows(a,[[2,2,2],[2,1,2]])


@policy('Print')
def printing(p):
    a = p.matrix(2,3); a[1,1] = 3; equal(a.text(),'0,0,0\n0,3,0\n')


def three_sum(a,factor=4):
    tmp = Vector([0]*a.data.rows)
    map_to(tmp.view(),[a.column(j) for j in range(3)],lambda x,y,z:x+y+z)
    map_to(a.column(2),[tmp.view()],lambda x:factor*x)


@policy('ArrayFunction')
def array_function(p):
    a = p.matrix(5,3,-1)
    for i in range(5):
        for j in range(3): a[i,j] = i-2+10*j
    f = p.function([a],three_sum); equal(f.invoke(a),0)
    for i in range(5): equal(a[i,2],96+12*i); equal(a[i,0],i-2); equal(a[i,1],i+8)
    b = p.matrix(3,3,-1); equal(f.invoke(b),0); check_rows(b,[[-1,-1,-12]]*3)


@policy('MultiMatrixArrayFunction','MultiMatrixDifferentRowsFromCreation')
def two_matrices(p):
    a,b = p.matrix(3,2),p.matrix(3,3)
    def operation(a,b,factor=1):
        tmp = Vector([0]*a.data.rows)
        map_to(tmp.view(),[a.column(0),b.column(2)],lambda x,y:x+y)
        map_to(a.column(1),[tmp.view()],lambda x:x*factor)
    for i in range(3):
        for j in range(2): a[i,j] = i+10*j
        for j in range(3): b[i,j] = i*2+20*j
        b[i,2] = i*4
    f = p.function([a,b],operation,[False,True]); equal(f.invoke(a,b),0)
    for i in range(3): equal(a[i,1],5*i); equal(a[i,0],i); equal(b[i,0],2*i)
    f2 = p.function([a,b],lambda x,y:operation(x,y,2),[False,True])
    for rows in (5,2):
        x,y = p.matrix(rows,2),p.matrix(rows,3)
        for i in range(rows): x[i,0] = i+1 if rows == 5 else 10*(i+1); y[i,2] = i*10 if rows == 5 else 5+10*i
        equal(f2.invoke(x,y),0)
        for i in range(rows): equal(x[i,1],((i+1)+10*i)*2 if rows == 5 else (10*(i+1)+5+10*i)*2)


@policy('VectorInMatrixFunction','MatrixVectorDifferentRowsFromCreation')
def matrix_vector(p):
    a = p.matrix(5,3,-1); v = Vector([100,200,300,400,500])
    for i in range(5):
        for j in range(3): a[i,j] = i-2+10*j
    def operation(a,v):
        tmp = Vector([0]*a.data.rows)
        map_to(tmp.view(),[a.column(j) for j in range(3)],lambda x,y,z:x+y+z)
        map_to(a.column(2),[v.view(),tmp.view()],lambda x,y:x*y)
    f = p.function([a,v],operation,[False,True]); equal(f.invoke(a,v),0)
    for i in range(5): equal(a[i,2],100*(i+1)*(24+3*i)); equal(a[i,0],i-2); equal(a[i,1],i+8)
    x,v0 = p.matrix(3,3),Vector([0]*3)
    def flexible(x,v):
        tmp = Vector([0]*x.data.rows)
        map_to(tmp.view(),[x.column(0),v.view()],lambda a,b:a+b)
        map_to(x.column(1),[tmp.view()],lambda t:3*t)
    ff = p.function([x,v0],flexible,[False,True]); y,w = p.matrix(5,3),Vector([i*10 for i in range(5)])
    for i in range(5): y[i,0] = i+1
    equal(ff.invoke(y,w),0)
    for i in range(5): equal(y[i,1],3*(i+1+i*10))


@policy('MultipleTemporaries')
def temporaries(p):
    a = p.matrix(4,5)
    for i in range(4): a[i,0] = i+1; a[i,1] = (i+1)*10
    def operation(a):
        t1,t2 = Vector([0]*a.data.rows),Vector([0]*a.data.rows)
        map_to(t1.view(),[a.column(0),a.column(1)],lambda x,y:x*y)
        map_to(t2.view(),[a.column(0),a.column(1)],lambda x,y:x+y)
        for j,fn in [(2,lambda x,y:x+y),(3,lambda x,y:x-y),(4,lambda x,y:x*y)]: map_to(a.column(j),[t1.view(),t2.view()],fn)
    f = p.function([a],operation); equal(f.invoke(a),0)
    for i in range(4):
        product,total = 10*(i+1)**2,11*(i+1)
        equal(a[i,2],product+total); equal(a[i,3],product-total); equal(a[i,4],product*total)


@policy('ColumnViewReuse')
def column_reuse(p):
    a = p.matrix(3,4)
    for i in range(3): a[i,0] = i+1
    def operation(a):
        c0,c1,c2,c3 = [a.column(j) for j in range(4)]
        map_to(c1,[c0],lambda x:2*x); map_to(c2,[c0,c1],lambda x,y:x+y); map_to(c3,[c2,c1],lambda x,y:x*y)
    f = p.function([a],operation); equal(f.invoke(a),0)
    for i in range(3): equal(a[i,1],2*(i+1)); equal(a[i,2],3*(i+1)); equal(a[i,3],6*(i+1)**2)


@policy('FunctionReusability')
def reuse(p):
    a,b,c = p.matrix(2,3,1),p.matrix(2,3,5),p.matrix(2,3)
    for i in range(2):
        for j in range(3): a[i,j] = i+j
        c[i,0] = i*10
    f = p.function([a],lambda a:three_sum(a,2))
    for obj,expected in [(a,[6,12]),(b,[30,30]),(c,[0,20])]:
        equal(f.invoke(obj),0)
        for i,v in enumerate(expected): equal(obj[i,2],v)


@policy('MismatchedRowsAtInvocation','MultipleMatricesMismatchedRowsAtInvocation','WrongColumnCountAtInvocation',
        'MismatchedRowDimensions','MismatchedColumnDimensions','WrongMatrixDimensions')
def dimensions(p):
    a,b = p.matrix(3,3,1),p.matrix(4,3,2)
    f = p.function([a,b],lambda a,b:map_to(a.column(1),[a.column(0),b.column(0)],lambda x,y:x+y),[False,True])
    equal(f.invoke(a,b),4)
    c,d = p.matrix(5,3),p.matrix(5,3,1); equal(f.invoke(c,d),0)
    v = Vector([0]*3)
    f2 = p.function([a,v],lambda a,v:None,[False,True])
    equal(f2.invoke(c,v),4); equal(f2.invoke(a,Vector([0]*5)),4)
    f3 = p.function([a],lambda a:None)
    for cols in (4,2): equal(f3.invoke(p.matrix(4,cols)),4)
    equal(f3.invoke(p.matrix(7,3)),0)
    invalid = a.column(5)
    equal(lib.chem_view_fill(C.byref(invalid),1),1)
    x,y,z = p.matrix(3,4,1),p.matrix(3,5,2),p.matrix(7,4,1)
    f4 = p.function([x],lambda a:map_to(a.column(3),[a.column(0)],lambda x:2*x))
    equal(f4.invoke(x),0); equal(f4.invoke(y),4); equal(f4.invoke(z),0)


@policy('ConstMatrixFunction','EmptyMatrixFunction')
def const_empty(p):
    a = p.matrix(3,4)
    for i in range(3):
        for j in range(4): a[i,j] = i*10+j
    def operation(a):
        t = Vector([0]*a.data.rows)
        map_to(t.view(),[a.column(0),a.column(1)],lambda x,y:x+y)
        reduce([a.column(2,True)],lambda x,y:x+y)
    f = p.function([a],operation,[True]); equal(f.invoke(a),0); equal(a[0,0],0); equal(a[1,2],12)
    x,y = p.matrix(0,3,1),p.matrix(3,0,1)
    fx = p.function([x],lambda a:map_to(a.column(0),[],lambda:99)); equal(fx.invoke(x),0)
    fy = p.function([y],lambda a:None); equal(fy.invoke(y),0)


@policy('VectorTooSmall','VectorTooLarge','EmptyVectorNonEmptyMatrix','NonEmptyVectorEmptyMatrix','EmptyVectorEmptyMatrix',
        'MultipleVectorsDifferentSizes','MultipleMatricesDifferentRowsVector','VectorSizeMatchesOneMatrixOnly',
        'FunctionInvocationWithWrongSizedVector')
def vector_dimensions(p):
    for rows,vecsize,expected in [(5,3,4),(5,10,4),(5,0,4),(0,5,4),(0,0,0)]:
        a,v = p.matrix(rows,3,1),Vector([0]*vecsize)
        f = p.function([a,v],lambda a,v:map_to(a.column(0),[v.view()],lambda x:x),[False,True])
        equal(f.invoke(a,v),expected)
    a,b,v,w = p.matrix(5,3,1),p.matrix(5,2),Vector([0]*5),Vector([0]*3)
    f = p.function([a,v,w],lambda *args:None); equal(f.invoke(a,v,w),4)
    c,d = p.matrix(4,2),p.matrix(5,3)
    f = p.function([c,d,Vector([0]*4)],lambda *args:None); equal(f.invoke(c,d,Vector([0]*4)),4)
    f = p.function([a,b,Vector([0]*4)],lambda *args:None); equal(f.invoke(a,b,Vector([0]*4)),4)
    e,v = p.matrix(3,2),Vector([0]*3)
    f = p.function([e,v],lambda *args:None); equal(f.invoke(e,v),0); equal(f.invoke(e,Vector([0]*5)),4)


@policy('MultipleVectorsSameSize','MultipleMatricesOneVector','ConstVector','MutableVector','FunctionReusabilityWithVectors',
        'MixedVectorColumnViewRowVariable','IntegerVector','FunctionWithConstSignature')
def vectors(p):
    a,v,w = p.matrix(5,3),Vector([1,2,3,4,5]),Vector([10,20,30,40,50])
    f = p.function([a,v,w],lambda a,v,w:map_to(a.column(0),[v.view(),w.view()],lambda x,y:x+y),[False,True,True])
    equal(f.invoke(a,v,w),0)
    for i in range(5): equal(a[i,0],11*(i+1))
    a,b,v = p.matrix(4,2),p.matrix(4,3),Vector([0,2,4,6])
    for i in range(4): a[i,0] = i+1; b[i,0] = i*10
    def two(a,b,v):
        for obj in (a,b): map_to(obj.column(1),[obj.column(0),v.view()],lambda x,y:x+y)
    f = p.function([a,b,v],two,[False,False,True]); equal(f.invoke(a,b,v),0)
    for i in range(4): equal(a[i,1],i+1+2*i); equal(b[i,1],12*i)
    a,v = p.matrix(3,2),Vector([10,20,30])
    f = p.function([a,v],lambda a,v:map_to(a.column(0),[v.view()],lambda x:2*x),[False,True]); equal(f.invoke(a,v),0)
    for i in range(3): equal(a[i,0],20*(i+1))
    for i in range(3): a[i,0] = i+1
    v = Vector([5,10,15]); f = p.function([a,v],lambda a,v:map_to(v.view(),[a.column(0)],lambda x:3*x),[True,False]); equal(f.invoke(a,v),0)
    equal(list(v.values[:3]),[3,6,9])
    f = p.function([a,v],lambda a,v:map_to(a.column(0),[v.view()],lambda x:10*x),[False,True])
    for values in ([1,2,3],[5,6,7]):
        v = Vector(values); equal(f.invoke(a,v),0)
        for i,val in enumerate(values): equal(a[i,0],val*10)
    a,v = p.matrix(4,3),Vector([100,200,300,400])
    for i in range(4): a[i,0] = i+1; a[i,1] = 10*(i+1)
    def mixed(a,v):
        tmp = Vector([0]*a.data.rows)
        map_to(tmp.view(),[a.column(0),a.column(1)],lambda x,y:x+y)
        map_to(a.column(2),[tmp.view(),v.view()],lambda x,y:x+y)
    f = p.function([a,v],mixed,[False,True]); equal(f.invoke(a,v),0)
    for i in range(4): equal(a[i,2],111*(i+1))


@policy('Fill','Copy')
def fill_copy(p):
    a,v = p.matrix(3,2),Vector([0]*3)
    col = a.column(1); equal(lib.chem_view_fill(C.byref(col),3.2),0)
    for i in range(3): equal(a[i,0],0); equal(a[i,1],3.2)
    vv = v.view(); equal(lib.chem_view_fill(C.byref(vv),3.2),0); equal(list(v.values[:3]),[3.2]*3)
    tmp = Vector([0]*3); tv = tmp.view(); equal(lib.chem_view_fill(C.byref(tv),9.9),0)
    map_to(a.column(0),[tv],lambda x:x)
    for i in range(3): equal(a[i,0],9.9)
    v = Vector([3.2,4.2,1.3]); vv = v.view(True); dest = a.column(1)
    equal(lib.chem_view_copy(C.byref(dest),C.byref(vv)),0)
    out = Vector([0]*3); ov = out.view(); source = a.column(1,True)
    equal(lib.chem_view_copy(C.byref(ov),C.byref(source)),0); equal(list(out.values[:3]),[3.2,4.2,1.3])
    for readonly in (False,True):
        dest,source = a.column(0),a.column(1,readonly)
        equal(lib.chem_view_copy(C.byref(dest),C.byref(source)),0)
        for i,val in enumerate([3.2,4.2,1.3]): equal(a[i,0],val); equal(a[i,1],val)
    tmp = Vector([0]*3); tv = tmp.view(); source = a.column(1,True)
    equal(lib.chem_view_copy(C.byref(tv),C.byref(source)),0)
    dest = a.column(1); equal(lib.chem_view_fill(C.byref(dest),0),0); map_to(dest,[tv],lambda x:x)
    for i,val in enumerate([3.2,4.2,1.3]): equal(a[i,1],val)


@policy('ReduceSum','ReduceMax','ReduceLOr','ReduceLAnd','ReduceStrict')
def reductions(p):
    a = p.matrix(3,2)
    for i in range(3):
        for j in range(2): a[i,j] = 1+i*2+j
    equal(reduce([a.column(0,True),a.column(1,True)],lambda acc,v:acc+v*v),91)
    for i,row in enumerate([[1,8],[5,2],[3,4]]):
        for j,v in enumerate(row): a[i,j] = v
    equal(reduce([a.column(0,True),a.column(1,True)],max,-sys.float_info.max),8)
    for large in (False,True):
        a.fill(1)
        if large: a[2,1] = 42
        equal(reduce([a.column(0,True),a.column(1,True)],lambda acc,v:bool(acc) or v>10),float(large))
    for nan in (False,True):
        a.fill(1)
        if nan: a[1,0] = math.nan
        equal(reduce([a.column(0,True),a.column(1,True)],lambda acc,v:bool(acc) and math.isfinite(v),1),float(not nan))
    b = p.matrix(3,1)
    equal(reduce([b.column(0,True)],lambda acc,v:acc+1,strict=True),3)
    b.fill(99)
    equal(reduce([b.column(0,True)],lambda acc,v:acc+1,strict=True),3)


@policy('VectorCapture')
def capture(p):
    a = p.matrix(3,3); order,values = [2,0,1],[10,20,30]
    def operation(a):
        for i in order: map_to(a.column(0),[a.column(0)],lambda x:x+values[i])
        map_to(a.column(1),[a.column(1)],lambda x:x+values[order[1]])
    f = p.function([a],operation); equal(f.invoke(a),0); check_rows(a,[[60,10,0]]*3)


def run_dense(identifier):
    family,fullname = identifier.split(':')
    name = fullname.split('.')[-1].replace('VectorMatrix','Matrix')
    if name.startswith('Test'): name = name[4:]
    for group in ((1,) if family == 'matrix' else (1,2,3,4)):
        with Pool(group) as pool: POLICIES[name](pool)


for _line in (ROOT/'case-list.txt').read_text().splitlines():
    if _line.startswith(('matrix:','vector_matrix:')):
        _id = _line.split()[0]
        _name = _id.split('.')[-1].replace('VectorMatrix','Matrix')
        if _name.startswith('Test'): _name = _name[4:]
        if _name in POLICIES: case(_id)(run_dense)
