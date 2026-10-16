"""Bindings for the independent C array/view and reusable dispatch APIs."""
from ffi import *
from ffi import bind as ffi_bind


def bind(name,result,*args):
    globals()[name] = ffi_bind(name,result,*args)


class View(C.Structure):
    _fields_ = [('data',P),('size',Z),('capacity',Z),('stride',Z),('type',I),('read_only',I),('category',I)]


class Group(C.Structure):
    _fields_ = [('array',P),('start',Z),('size',Z),('capacity',Z),('read_only',I)]


Element = C.CFUNCTYPE(None,DP,Z,P)
Reduction = C.CFUNCTYPE(None,D,DP,P)
Callback = C.CFUNCTYPE(I,Z,C.POINTER(Group),P)
AP = C.POINTER(P)
VP = C.POINTER(View)
GP = C.POINTER(Group)
bind('pattern_set',I,Z,BP,Z,Z)
bind('array_dense',P,Z,Z,Z,I,D)
bind('array_sparse',P,Z,Z,BP,I,Z,I,D)
bind('array_vector',P,Z,I,D)
bind('array_clone',P,P)
bind('array_free',None,P)
for name in ('records','elements','storage_size','group_size','groups','row_stride','column_stride'):
    bind('array_'+name,Z,P)
bind('array_data',P,P)
bind('array_element',I,P,Z,Z)
bind('array_index',I,P,Z,Z,C.POINTER(Z))
bind('array_get',I,P,Z,Z,DP)
bind('array_set',I,P,Z,Z,D)
bind('array_set_row',I,P,Z,DP,Z)
for name in ('fill','max','min','add_diagonal'):
    bind('array_'+name,I,P,D)
bind('array_copy',I,P,P)
bind('array_axpy',I,P,D,P)
bind('array_format',Z,P,P,Z)
bind('array_format_dense',Z,P,P,Z)
bind('array_group',I,P,Z,I,GP)
bind('group_view',I,GP,Z,VP)
bind('group_sparse_view',I,GP,Z,Z,VP)
bind('temporary_view',View,DP,Z,Z)
bind('view_get',D,VP,Z)
bind('view_set',I,VP,Z,D)
bind('view_fill',I,VP,D,I)
bind('view_copy',I,VP,VP,I)
bind('views_map',I,Z,VP,Element,P,I)
bind('view_reduce',I,VP,Reduction,DP,P,I)
bind('array_map',I,Z,AP,BP,Element,P)
bind('array_function_create',I,Z,AP,BP,Callback,P,C.POINTER(P))
bind('array_function_call',I,P,Z,AP)
bind('array_function_free',None,P)


class Array:
    def __init__(self,rows,columns=1,group=0,initial=0,integer=False,pattern=None,order=0,vector=False):
        self.rows,self.columns,self.group,self.integer = rows,columns,group,integer
        if vector:
            self.handle = array_vector(rows,integer,initial)
        elif pattern is None:
            self.handle = array_dense(rows,columns,group,integer,initial)
        else:
            data = np.zeros((columns,columns),dtype=np.uint8)
            for r,c in pattern:
                data[r,c] = 1
            self.handle = array_sparse(rows,columns,data.ctypes.data_as(BP),order,group,integer,initial)
        assert self.handle
        self.elements = array_elements(self.handle)

    def close(self):
        if self.handle:
            array_free(self.handle)
            self.handle = None

    def __enter__(self): return self
    def __exit__(self,*_): self.close()
    def element(self,r,c): return require(array_element(self.handle,r,c))
    def __getitem__(self,key):
        r,e = key
        value = D()
        require(array_get(self.handle,r,e,C.byref(value)))
        return value.value
    def __setitem__(self,key,value):
        r,e = key
        require(array_set(self.handle,r,e,value))
    def values(self):
        return np.array([[self[r,e] for e in range(self.elements)] for r in range(self.rows)])
    def initialize(self,values):
        values = np.asarray(values).reshape((self.rows,self.elements))
        for r in range(self.rows):
            require(array_set_row(self.handle,r,doubles(values[r]),self.elements))
        return self
    @property
    def raw(self):
        p = C.cast(array_data(self.handle), IP if self.integer else DP)
        return np.ctypeslib.as_array(p,shape=(array_storage_size(self.handle),))
    def view_group(self,index,read_only=0):
        result = Group()
        require(array_group(self.handle,index,read_only,C.byref(result)))
        return result
    def text(self):
        n = array_format(self.handle,None,0)
        buf = C.create_string_buffer(n+1)
        assert array_format(self.handle,buf,n+1) == n
        return buf.value.decode()


def column(group,index,read_only=None):
    v = View()
    require(group_view(C.byref(group),index,C.byref(v)))
    if read_only is not None:
        assert not group.read_only or read_only
        v.read_only = read_only
    return v


def temporary(group):
    data = doubles([0]*group.capacity)
    v = temporary_view(data,group.size,group.capacity)
    v.backing = data
    return v


def assign(destination, inputs, expression, strict=0):
    views = []
    for source in inputs:
        view = View.from_buffer_copy(source)
        view.read_only = 1
        views.append(view)
    views.append(destination)
    failures = []
    def evaluate(values,count,_):
        try: values[count-1] = expression(*values[:count-1])
        except Exception as error: failures.append(error)
    require(views_map(len(views),(View*len(views))(*views),Element(evaluate),None,strict))
    if failures: raise failures[0]


class Function:
    def __init__(self,arrays,operation,readonly=None):
        self.handle,self.failures = P(),[]
        def invoke(n,groups,_):
            try:
                result = operation(list(groups[:n]))
                return result if isinstance(result,int) else 0
            except Exception as error:
                self.failures.append(error)
                return -1
        self.callback = Callback(invoke)
        read_only = (C.c_ubyte*len(arrays))(*(readonly or [0]*len(arrays)))
        require(array_function_create(len(arrays),(P*len(arrays))(*(a.handle for a in arrays)),
                    read_only,self.callback,None,C.byref(self.handle)))

    def call(self,*arrays):
        status = array_function_call(self.handle,len(arrays),(P*len(arrays))(*(a.handle for a in arrays)))
        if self.failures: raise self.failures.pop(0)
        return status

    def __call__(self,*arrays): require(self.call(*arrays))
    def close(self):
        if self.handle:
            array_function_free(self.handle)
            self.handle = P()
    def __enter__(self): return self
    def __exit__(self,*_): self.close()


def reduce_view(view,accumulator,expression,strict=False):
    def evaluate(value,acc,_): acc[0] = expression(value,acc[0])
    require(view_reduce(C.byref(view),Reduction(evaluate),C.byref(accumulator),None,strict))
