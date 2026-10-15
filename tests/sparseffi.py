"""Sparse CSR/CSC bindings; data storage uses the shared C array/view API."""
from arrayffi import *
from chemffi import MatrixIndex

bind('chem_sparse_create',I,Z,Z,Z,I,I,Z,C.POINTER(MatrixIndex),P,C.POINTER(P))
bind('chem_sparse_free',None,P)
bind('chem_sparse_clone',I,P,C.POINTER(P))
for name in ['order','nonzeros']:
    bind('chem_sparse_'+name,Z,P)
bind('chem_sparse_compression',I,P)
bind('chem_sparse_values',P,P)
for name in ['outer_starts','inner_indices']:
    bind('chem_sparse_'+name,C.POINTER(Z),P)
bind('chem_sparse_flat_index',Z,P,Z,Z)
bind('chem_sparse_index',Z,P,Z,Z,Z)
bind('chem_sparse_single_index',Z,P,Z,Z)
for name in ['get','set']:
    bind('chem_sparse_'+name,I,P,Z,Z,Z,P)
bind('chem_sparse_diagonal_indices',Z,P,Z,C.POINTER(Z),Z)
bind('chem_sparse_add_diagonal',I,P,P)
bind('chem_sparse_block_view',I,P,Z,Z,I,C.POINTER(View))
bind('chem_sparse_index_view',I,P,Z,I,C.POINTER(View))
bind('chem_sparse_format',Z,P,P,Z,I)
MISSING=Z(-1).value


class Sparse:
    def __init__(self,order,blocks=1,entries=(),group=0,compression=0,type=0,value=0,status=0):
        self.handle=P();self.type=type;self.owned=True
        coordinates=(MatrixIndex*len(entries))(*[MatrixIndex(r,c) for r,c in entries])
        initial=TYPES[type](value)
        result=lib.chem_sparse_create(order,blocks,group,compression,type,len(entries),coordinates,C.byref(initial),C.byref(self.handle))
        assert result==status,result
        self.values=Array.borrowed(lib.chem_sparse_values(self.handle)) if result==0 else None
        if self.values:self.values.type=type
    def __enter__(self):return self
    def __exit__(self,*args):self.close()
    def close(self):
        if self.owned:lib.chem_sparse_free(self.handle)
        self.handle=P()
    @property
    def order(self):return lib.chem_sparse_order(self.handle)
    @property
    def nnz(self):return lib.chem_sparse_nonzeros(self.handle)
    @property
    def blocks(self):return self.values.rows
    @property
    def group(self):return self.values.group
    @property
    def data(self):return self.values.data
    def flat(self,row,col):return lib.chem_sparse_flat_index(self.handle,row,col)
    def index(self,block,row,col):return lib.chem_sparse_index(self.handle,block,row,col)
    def __getitem__(self,key):
        scalar=TYPES[self.type]();ok(lib.chem_sparse_get(self.handle,*key,C.byref(scalar)));return scalar.value
    def __setitem__(self,key,value):
        scalar=TYPES[self.type](value);ok(lib.chem_sparse_set(self.handle,*key,C.byref(scalar)))
    def fill(self,value):self.values.fill(value)
    def diagonal(self,block):
        count=lib.chem_sparse_diagonal_indices(self.handle,block,None,0);assert count!=MISSING
        data=(Z*count)();assert lib.chem_sparse_diagonal_indices(self.handle,block,data,count)==count
        return list(data)
    def add_diagonal(self,value):
        scalar=TYPES[self.type](value);ok(lib.chem_sparse_add_diagonal(self.handle,C.byref(scalar)))
    def block_view(self,row,col,ro=False):
        out=View();ok(lib.chem_sparse_block_view(self.handle,row,col,ro,C.byref(out)));return out
    def index_view(self,index,ro=False):
        out=View();ok(lib.chem_sparse_index_view(self.handle,index,ro,C.byref(out)));return out
    def format(self,nonzeros=False):
        size=lib.chem_sparse_format(self.handle,None,0,nonzeros);assert size!=MISSING
        data=C.create_string_buffer(size+1)
        assert lib.chem_sparse_format(self.handle,data,len(data),nonzeros)==size
        return data.value.decode()
    def clone(self):
        obj=object.__new__(Sparse);obj.handle=P();obj.type=self.type;obj.owned=True
        ok(lib.chem_sparse_clone(self.handle,C.byref(obj.handle)))
        obj.values=Array.borrowed(lib.chem_sparse_values(obj.handle));obj.values.type=self.type
        return obj
