"""C slices, mutable/const views, and grouped storage equivalence checks."""
import re
from arrays import *
from harness import case, equal, NotApplicable
from matrix_cases import Pool, map_to
from sparse_cases import SparsePool


def slice_(view,start,count):
    view.start,view.count = start,count
    return view


@case('grouped_views:GroupedView.CategoryTagsAreDistinct',
      'view_category_sfinae:ViewCategorySFINAE.VectorTypes','view_category_sfinae:ViewCategorySFINAE.MatrixViews',
      'view_category_sfinae:ViewCategorySFINAE.SparseMatrixViews','view_category_sfinae:ViewCategorySFINAE.MatrixTypesExcluded')
def type_categories(_):
    raise NotApplicable('C has no C++ concepts, template category tags, or overload resolution')


@case('view_category_sfinae:ViewCategorySFINAE.RuntimeVectorUse')
def runtime_vector(_):
    vector = Vector([1,2,3,4,5])
    equal(reduce([vector.view(True)],lambda acc,v:acc+v),15)


@case('grouped_views:GroupedView.MatrixColumnViewRoundTrip','grouped_views:GroupedView.VectorMatrixColumnViewContiguousBlock')
def dense_roundtrip(identifier):
    vector = 'VectorMatrix' in identifier
    with Pool(4 if vector else 1) as p:
        a = p.matrix(8 if vector else 3,3 if vector else 4)
        for i in range(a.data.rows):
            for j in range(a.data.columns): a[i,j] = (100 if vector else 0)+10*i+j
        first,count = (4,4) if vector else (1,1)
        dest,source = slice_(a.column(1 if vector else 2),first,count),slice_(a.column(2 if vector else 3,True),first,count)
        offset = lib.chem_array_index(a.p,first,dest.column)
        for i in range(count): equal(a.data.values[offset+i],(100 if vector else 0)+10*(first+i)+dest.column)
        map_to(dest,[source],lambda v:2*v)
        for i in range(a.data.rows):
            expected = 2*((100 if vector else 0)+10*i+source.column) if first <= i < first+count else (100 if vector else 0)+10*i+dest.column
            equal(a[i,dest.column],expected)


@case(*['grouped_views:GroupedView.DenseGroupedMatchesRawL'+str(g) for g in (1,2,4,8)],'grouped_views:GroupedView.PlainMatrixGroupedRoundTrip')
def dense_equivalence(identifier):
    match = re.search(r'L(\d+)$',identifier); group = int(match[1]) if match else 1
    plain = match is None
    with Pool(group) as p:
        rows,cols = (4,3) if plain else (3*group,4)
        a,b = p.matrix(rows,cols),p.matrix(rows,cols)
        for i in range(rows):
            for j in range(cols): a[i,j] = b[i,j] = i+j*5 if plain else i*10+j
        if plain:
            fn = lambda a:map_to(a.column(2),[a.column(0),a.column(1)],lambda x,y:x-y)
            for i in range(rows): b[i,2] = b[i,0]-b[i,1]
        else:
            fn = lambda a:map_to(a.column(3),[a.column(0)],lambda x:3*x+1)
            for i in range(rows): b[i,3] = b[i,0]*3+1
        f = p.function([a],fn); equal(f.invoke(a),0)
        for i in range(rows):
            for j in range(cols): equal(a[i,j],b[i,j])


@case('grouped_views:GroupedView.SparseStandardBlockViewRoundTrip','grouped_views:GroupedView.SparseVectorBlockViewContiguousBlock',
      'grouped_views:GroupedView.SparseVectorRowColOverloadStillFunctional','grouped_views:GroupedView.ConstGroupViewProducesGroupedViews')
def sparse_roundtrip(identifier):
    standard = 'Standard' in identifier
    with SparsePool(1 if standard else 4,False) as p:
        if standard:
            a = p.sparse(2,3,[(0,1),(1,2),(2,2)])
            for b,factor in enumerate([1,10]):
                for j,pair in enumerate([(0,1),(1,2),(2,2)]): a[(b,)+pair] = factor*(j+1)
            dest,source = slice_(a.block(0,1),1,1),slice_(a.block(1,2,True),1,1)
            map_to(dest,[source],lambda v:v+5); equal(a[1,0,1],25); equal(a[0,0,1],1)
        elif 'RowColOverload' in identifier:
            a = p.sparse(8,3,[(0,1)])
            for b in range(8): a[b,0,1] = b+1
            equal(reduce([slice_(a.block(0,1,True),4,4)],lambda a,b:a+b),26)
        elif 'ConstGroup' in identifier:
            a,dense = p.sparse(4,3,[(0,1)]),p.matrix(4,2)
            for b in range(4): a[b,0,1] = b+1; dense[b,1] = 100+b
            for b in range(4): equal(a.data.values[b],b+1); equal(dense.data.values[4+b],100+b)
            v = a.block(0,1,True); equal(lib.chem_view_fill(C.byref(v),0),1)
        else:
            a = p.sparse(8,3,[(0,1),(1,2)])
            for b in range(8): a[b,0,1] = 100+b; a[b,1,2] = 200+b
            dest,source = slice_(a.block(0,1),4,4),slice_(a.block(1,2,True),4,4)
            for b in range(4):
                equal(a.data.values[lib.chem_sparse_index(a.p,4,0,1)+b],104+b)
                equal(a.data.values[lib.chem_sparse_index(a.p,4,1,2)+b],204+b)
            map_to(dest,[dest,source],lambda x,y:x+y)
            for b in range(8): equal(a[b,0,1],300+2*b if b >= 4 else 100+b)


@case(*['grouped_views:GroupedView.MixedSparseAndDenseGroupedL'+str(g) for g in (1,2,4,8)])
def mixed(identifier):
    group = int(identifier[-1])
    with SparsePool(group,False) as p:
        a,d = p.sparse(2*group,3,[(0,1)]),p.matrix(2*group,3)
        for b in range(2*group):
            a[b,0,1] = b+1
            for c in range(3): d[b,c] = 10+b+c
        f = p.function([a,d],lambda a,d:map_to(d.column(0),[d.column(1),a.block(0,1)],lambda x,y:x+y),[True,False])
        equal(f.invoke(a,d),0)
        for b in range(2*group):
            for c in range(3): equal(d[b,c],12+2*b if c == 0 else 10+b+c)


@case('grouped_views:GroupedView.FillAndCopySparseStandardCSR','grouped_views:GroupedView.FillAndCopySparseStandardCSC',
      'grouped_views:GroupedView.FillAndCopySparseVectorCSR_L4','grouped_views:GroupedView.FillAndCopySparseVectorCSC_L4',
      'grouped_views:GroupedView.FillAndCopySparseVectorCSR_L8','grouped_views:GroupedView.FillCopyCrossGroupViewSources')
def fill_copy(identifier):
    group = 8 if identifier.endswith('L8') else (1 if 'Standard' in identifier else 4)
    column = 'CSC' in identifier
    with SparsePool(group,column) as p:
        a = p.sparse(group+1,3,[(0,1),(2,0)])
        for b in range(group+1): a[b,0,1] = b+1; a[b,2,0] = 100+b
        start,value_ = (1,9) if group == 1 and column else (0,42)
        dest = slice_(a.block(0,1),start,group); source = slice_(a.block(0,1,True),start,group)
        equal(lib.chem_view_fill(C.byref(dest),value_),0)
        for b in range(start,start+group): equal(a[b,0,1],value_); equal(a[b,2,0],100+b)
        if not start: equal(a[group,0,1],group+1)
        to = slice_(a.block(2,0),start,group); equal(lib.chem_view_copy(C.byref(to),C.byref(source)),0)
        for b in range(start,start+group): equal(a[b,2,0],value_)
        if group > 1:
            last = slice_(a.block(0,1),group,group); equal(lib.chem_view_fill(C.byref(last),-7),0); equal(a[group,0,1],-7)
        if 'CrossGroup' in identifier:
            src,dst,ref = [p.sparse(4,2,[(0,1)]) for _ in range(3)]
            for b in range(4): src[b,0,1] = b+1; ref[b,0,1] = 999
            to,source = dst.block(0,1),src.block(0,1,True); equal(lib.chem_view_copy(C.byref(to),C.byref(source)),0)
            for b in range(4): equal(dst[b,0,1],b+1)
            map_to(ref.block(0,1),[],lambda:55.5); equal(lib.chem_view_fill(C.byref(to),55.5),0)
            for b in range(4): equal(dst[b,0,1],ref[b,0,1])
