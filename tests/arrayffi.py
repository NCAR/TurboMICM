"""C array bindings used by matrix and grouped-view reference case ports."""
from chemffi import C, lib, bind, I, D, P, Z

TYPES=[D,I,C.c_bool]

class View(C.Structure):
    _fields_=[('data',P)]+[(name,Z) for name in ('length','capacity','block_length','block_stride','element_stride')]+[('type',I),('read_only',I)]

class Arg(C.Structure):
    _fields_=[('array',P),('read_only',I)]

ELEMENT=C.CFUNCTYPE(None,Z,C.POINTER(P),P)
REDUCE=C.CFUNCTYPE(None,Z,C.POINTER(P),P,P)
BODY=C.CFUNCTYPE(I,Z,C.POINTER(Arg),P)
bind('chem_array_create',I,Z,Z,Z,I,P,C.POINTER(P))
bind('chem_array_free',None,P)
bind('chem_array_clone',I,P,C.POINTER(P))
bind('chem_array_move',I,C.POINTER(P),C.POINTER(P))
for suffix in ['rows','columns','group_length','groups','storage_size']:
    bind('chem_array_'+suffix,Z,P)
bind('chem_array_index',Z,P,Z,Z)
bind('chem_array_data',P,P)
bind('chem_array_const_data',P,P)
for name in ['get','set']:
    bind('chem_array_'+name,I,P,Z,Z,P)
bind('chem_array_fill',I,P,P)
bind('chem_array_copy',I,P,P)
bind('chem_array_axpy',I,P,D,P)
bind('chem_array_limit',I,P,D,I)
for name in ['assign_row','read_row']:
    bind('chem_array_'+name,I,P,Z,P,Z)
bind('chem_array_from_rows',I,Z,C.POINTER(Z),C.POINTER(P),Z,I,C.POINTER(P))
bind('chem_array_format',Z,P,P,Z)
for name in ['row','column']:
    bind('chem_array_'+name,I,P,Z,I,C.POINTER(View))
bind('chem_array_group_column',I,P,Z,Z,I,C.POINTER(View))
bind('chem_vector_view',View,P,Z,I,I)
bind('chem_view_element',P,C.POINTER(View),Z)
bind('chem_view_fill',I,View,P)
bind('chem_view_copy',I,View,View)
bind('chem_views_apply',I,Z,C.POINTER(View),Z,ELEMENT,P)
bind('chem_views_reduce',I,Z,C.POINTER(View),I,REDUCE,P,P)
bind('chem_array_function_create',I,Z,C.POINTER(Arg),BODY,P,C.POINTER(P))
bind('chem_array_function_call',I,P,Z,C.POINTER(Arg))
bind('chem_array_function_free',None,P)


def ok(status):
    assert status==0, status


class Array:
    def __init__(self,rows,columns,group=0,value=0,type=0):
        self.type=type;self.handle=P();self.owned=True;self.read_only=0
        initial=TYPES[type](value)
        ok(lib.chem_array_create(rows,columns,group,type,C.byref(initial),C.byref(self.handle)))

    @classmethod
    def borrowed(cls,ptr,readonly=0):
        obj=object.__new__(cls);obj.handle=P(ptr);obj.type=0;obj.owned=False;obj.read_only=readonly
        return obj

    @classmethod
    def from_rows(cls,rows,group=0,type=0,status=0):
        obj=object.__new__(cls);obj.type=type;obj.handle=P();obj.owned=True;obj.read_only=0
        lengths=(Z*len(rows))(*map(len,rows))
        data=[(TYPES[type]*len(row))(*row) for row in rows]
        pointers=(P*len(rows))(*[C.cast(a,P) for a in data])
        assert lib.chem_array_from_rows(len(rows),lengths,pointers,group,type,C.byref(obj.handle))==status
        return obj

    def __enter__(self):return self
    def __exit__(self,*args):self.close()
    def close(self):
        if self.owned:lib.chem_array_free(self.handle)
        self.handle=P()
    @property
    def rows(self):return lib.chem_array_rows(self.handle)
    @property
    def columns(self):return lib.chem_array_columns(self.handle)
    @property
    def group(self):return lib.chem_array_group_length(self.handle)
    @property
    def size(self):return lib.chem_array_storage_size(self.handle)
    @property
    def data(self):return C.cast(lib.chem_array_data(self.handle),C.POINTER(TYPES[self.type]))
    def __getitem__(self,key):
        if not isinstance(key,tuple):key=(key,0)
        value=TYPES[self.type]();ok(lib.chem_array_get(self.handle,*key,C.byref(value)));return value.value
    def __setitem__(self,key,value):
        if not isinstance(key,tuple):key=(key,0)
        scalar=TYPES[self.type](value);ok(lib.chem_array_set(self.handle,*key,C.byref(scalar)))
    def values(self):return [[self[i,j] for j in range(self.columns)] for i in range(self.rows)]
    def fill(self,value):
        scalar=TYPES[self.type](value);ok(lib.chem_array_fill(self.handle,C.byref(scalar)))
    def clone(self):
        obj=object.__new__(Array);obj.type=self.type;obj.handle=P();obj.owned=True;obj.read_only=0
        ok(lib.chem_array_clone(self.handle,C.byref(obj.handle)));return obj
    def assign_row(self,row,values):
        data=(TYPES[self.type]*len(values))(*values)
        return lib.chem_array_assign_row(self.handle,row,data,len(values))
    def read_row(self,row):
        data=(TYPES[self.type]*self.columns)();ok(lib.chem_array_read_row(self.handle,row,data,self.columns));return list(data)
    def col(self,column,ro=False):
        out=View();ok(lib.chem_array_column(self.handle,column,ro or self.read_only,C.byref(out)));return out
    def row(self,row,ro=False):
        out=View();ok(lib.chem_array_row(self.handle,row,ro or self.read_only,C.byref(out)));return out
    def group_col(self,group,column,ro=False):
        out=View();ok(lib.chem_array_group_column(self.handle,group,column,ro or self.read_only,C.byref(out)));return out
    def format(self):
        length=lib.chem_array_format(self.handle,None,0)
        data=C.create_string_buffer(length+1)
        assert lib.chem_array_format(self.handle,data,len(data))==length
        return data.value.decode()


def element(view,index):
    ptr=lib.chem_view_element(C.byref(view),index)
    assert ptr
    return C.cast(ptr,C.POINTER(TYPES[view.type]))


def apply(views,fn,writable=1,status=0):
    errors=[]
    def callback(n,pointers,context):
        try:fn([C.cast(pointers[j],C.POINTER(TYPES[views[j].type])) for j in range(n)])
        except Exception as error:errors.append(error)
    cb=ELEMENT(callback)
    assert lib.chem_views_apply(len(views),(View*len(views))(*views),writable,cb,None)==status
    if errors:raise errors[0]


def assign(destination,sources,fn):
    apply([destination]+list(sources),lambda p:p[0].__setitem__(0,fn(*[q[0] for q in p[1:]])))


def fill(view,value):
    scalar=TYPES[view.type](value);ok(lib.chem_view_fill(view,C.byref(scalar)))


def copy(destination,source):ok(lib.chem_view_copy(destination,source))


def reduce(views,fn,initial=0,strict=False,type=0):
    acc=TYPES[type](initial);errors=[]
    def callback(n,pointers,target,context):
        try:
            ptr=C.cast(target,C.POINTER(TYPES[type]))
            ptr[0]=fn(ptr[0],*[C.cast(pointers[j],C.POINTER(TYPES[views[j].type]))[0] for j in range(n)])
        except Exception as error:errors.append(error)
    cb=REDUCE(callback)
    ok(lib.chem_views_reduce(len(views),(View*len(views))(*views),strict,cb,C.byref(acc),None))
    if errors:raise errors[0]
    return acc.value


class Function:
    def __init__(self,arrays,body,readonly=None,status=0):
        self.handle=P();self.errors=[];self.types=[a.type for a in arrays]
        readonly=readonly or [0]*len(arrays)
        def callback(n,args,context):
            try:
                views=[Array.borrowed(args[i].array,args[i].read_only) for i in range(n)]
                for i,a in enumerate(views):a.type=self.types[i]
                result=body(*views)
                return result or 0
            except Exception as error:self.errors.append(error);return 1
        self.body=BODY(callback)
        args=(Arg*len(arrays))(*[Arg(a.handle,ro) for a,ro in zip(arrays,readonly)])
        assert lib.chem_array_function_create(len(arrays),args,self.body,None,C.byref(self.handle))==status
    def __enter__(self):return self
    def __exit__(self,*args):self.close()
    def close(self):lib.chem_array_function_free(self.handle);self.handle=P()
    def __call__(self,*arrays,readonly=None,status=0):
        readonly=readonly or [0]*len(arrays)
        args=(Arg*len(arrays))(*[Arg(a.handle,ro) for a,ro in zip(arrays,readonly)])
        result=lib.chem_array_function_call(self.handle,len(arrays),args)
        if self.errors:raise self.errors.pop(0)
        assert result==status,result
