"""Bindings for dense/sparse arrays and reusable C callback functions."""
from api import *


class ArrayData(C.Structure):
    _fields_ = [(name,Z) for name in ('rows','columns','group','allocated')] + [('values',DP), ('order',Z), ('index',C.POINTER(I)), ('column_order',I)]


AP = C.POINTER(ArrayData)


class View(C.Structure):
    _fields_ = [('array',AP), ('vector',DP), ('column',Z), ('start',Z), ('count',Z), ('readonly',I)]


class Argument(C.Structure):
    _fields_ = [('array',AP), ('vector',DP), ('count',Z), ('readonly',I)]


ElementFn = C.CFUNCTYPE(None,DP,Z,P)
ArrayFn = C.CFUNCTYPE(None,C.POINTER(Argument),Z,P)
ReduceFn = C.CFUNCTYPE(None,D,DP,P)
bind('chem_array_dense',AP,Z,Z,Z,D)
bind('chem_array_from_rows',I,C.POINTER(DP),ZP,Z,Z,C.POINTER(AP))
bind('chem_array_sparse',AP,Z,Z,Z,I,ZP,ZP,Z,D)
bind('chem_array_free',None,AP)
bind('chem_array_index',Z,AP,Z,Z)
bind('chem_sparse_index',Z,AP,Z,Z,Z)
bind('chem_array_get',I,AP,Z,Z,DP)
bind('chem_array_set',I,AP,Z,Z,D)
bind('chem_array_row',I,AP,Z,DP,Z,I)
bind('chem_array_fill',None,AP,D)
bind('chem_array_bound',None,AP,D,I)
bind('chem_array_copy',I,AP,AP)
bind('chem_array_axpy',I,AP,D,AP)
bind('chem_array_add_diagonal',None,AP,D)
bind('chem_array_format',Z,AP,I,C.c_char_p,Z)
bind('chem_array_view',View,AP,Z,I)
bind('chem_vector_view',View,DP,Z,I)
bind('chem_view_fill',I,C.POINTER(View),D)
bind('chem_view_copy',I,C.POINTER(View),C.POINTER(View))
bind('chem_view_for_each',I,C.POINTER(View),Z,ElementFn,P)
bind('chem_view_reduce',I,C.POINTER(View),I,I,DP)
bind('chem_view_reduce_function',I,C.POINTER(View),I,DP,ReduceFn,P)
bind('chem_function_create',I,C.POINTER(Argument),Z,ArrayFn,P,C.POINTER(P))
bind('chem_function_invoke',I,P,C.POINTER(Argument),Z)
bind('chem_function_free',None,P)


class Array:
    def __init__(self,rows=0,columns=0,group=1,initial=0.0,pattern=None,order=0,column_order=False,pointer=None):
        self.owned = pointer is None
        if pointer is not None:
            self.p = pointer
        elif pattern is None:
            self.p = lib.chem_array_dense(rows,columns,group,initial)
        else:
            self.p = lib.chem_array_sparse(rows,order,group,column_order,inds([r for r,c in pattern]),inds([c for r,c in pattern]),len(pattern),initial)
        assert self.p

    @property
    def data(self): return self.p.contents

    def close(self):
        if self.p and self.owned:
            lib.chem_array_free(self.p)
        self.p = None

    def __getitem__(self,index):
        if len(index) == 3:
            at = lib.chem_sparse_index(self.p,*index); assert at != MAX
            return self.data.values[at]
        value = D(); assert lib.chem_array_get(self.p,*index,C.byref(value)) == 0
        return value.value

    def __setitem__(self,index,value):
        if len(index) == 3:
            at = lib.chem_sparse_index(self.p,*index); assert at != MAX
            self.data.values[at] = value
        else:
            assert lib.chem_array_set(self.p,*index,value) == 0

    def column(self,index,readonly=False):
        return lib.chem_array_view(self.p,index,readonly)

    def block(self,row,column,readonly=False):
        index = self.data.index[row*self.data.order+column]; assert index >= 0
        return self.column(index,readonly)

    def fill(self,value): lib.chem_array_fill(self.p,value)

    def copy(self,other): assert lib.chem_array_copy(self.p,other.p) == 0

    def text(self,nonzero=False):
        n = lib.chem_array_format(self.p,nonzero,None,0)
        buffer = C.create_string_buffer(n+1)
        assert lib.chem_array_format(self.p,nonzero,buffer,n+1) == n
        return buffer.value.decode()


class Vector:
    def __init__(self,values=(),pointer=None,count=None):
        self.values = pointer if pointer is not None else arr(values)
        self.count = len(values) if count is None else count

    def view(self,readonly=False): return lib.chem_vector_view(self.values,self.count,readonly)


def foreach(views,function):
    errors = []
    def run(values,n,context):
        try:
            result = function(*values[:n])
            if result is not None:
                for i,v in enumerate(result): values[i] = v
        except Exception as error: errors.append(error)
    callback = ElementFn(run)
    assert lib.chem_view_for_each((View*len(views))(*views),len(views),callback,None) == 0
    if errors: raise errors[0]


def argument(obj,readonly=False):
    if isinstance(obj,Array): return Argument(obj.p,None,0,readonly)
    return Argument(None,obj.values,obj.count,readonly)


def reduce(views,function,initial=0.0,strict=False):
    acc = D(initial)
    callback = ReduceFn(lambda val,a,ctx: a.__setitem__(0,function(a[0],val)))
    for view in views:
        assert lib.chem_view_reduce_function(C.byref(view),strict,C.byref(acc),callback,None) == 0
    return acc.value


class Function:
    def __init__(self,objects,function,readonly=()):
        self.function = function; self.errors = []; self.p = P()
        self.readonly = readonly or [False]*len(objects)
        def run(arguments,n,context):
            try:
                objs = [Array(pointer=a.array) if a.array else Vector(pointer=a.vector,count=a.count) for a in arguments[:n]]
                function(*objs)
            except Exception as error: self.errors.append(error)
        self.callback = ArrayFn(run)
        args = (Argument*len(objects))(*(argument(obj,ro) for obj,ro in zip(objects,self.readonly)))
        self.status = lib.chem_function_create(args,len(objects),self.callback,None,C.byref(self.p))

    def invoke(self,*objects,readonly=None):
        ro = readonly or self.readonly
        args = (Argument*len(objects))(*(argument(obj,r) for obj,r in zip(objects,ro)))
        status = lib.chem_function_invoke(self.p,args,len(objects))
        if self.errors: raise self.errors.pop(0)
        return status

    def close(self):
        lib.chem_function_free(self.p); self.p = P()
