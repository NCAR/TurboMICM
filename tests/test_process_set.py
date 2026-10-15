"""Process-set forcing, declared pattern, Jacobian, and masking ports."""
from functools import partial
import math
from chemffi import *
from reference_random import uniform_integers, lognormals
from test_rate_store import multiplier
from test_jacobian import finite_difference, Sparse, Comparison

bind('chem_add_forcing',None,P,DP,DP)
bind('chem_subtract_jacobian',None,P,DP,DP)
bind('chem_kinetic_nonzero',I,P,Z,Z)
bind('chem_reaction_species_used',I,P,C.c_char_p)


def set_rates(m,rows):
    values = C.cast(lib.chem_rate_constants(m.state),DP)
    for c,row in enumerate(rows):
        for j,v in enumerate(row):values[j*m.nc+c] = v


def forcing(m,initial):
    result = (D*(m.n*m.nc))(*([initial]*(m.n*m.nc)))
    lib.chem_add_forcing(m.state,m.y,result)
    return [[result[i*m.nc+c] for i in range(m.n)] for c in range(m.nc)]


def jacobian(m,initial):
    result = (D*(lib.chem_jacobian_size(m.handle)*m.nc))(*([initial]*(lib.chem_jacobian_size(m.handle)*m.nc)))
    lib.chem_subtract_jacobian(m.state,m.y,result)
    matrices = []
    for c in range(m.nc):
        rows = []
        for i in range(m.n):
            row = []
            for j in range(m.n):
                index = lib.chem_jacobian_index(m.handle,i,j)
                row.append(0 if index == Z(-1).value else result[index*m.nc+c])
            rows.append(row)
        matrices.append(rows)
    return matrices


def real_eq(a,b):assert abs(a-b) <= 4*max(math.ulp(a),math.ulp(b)),(a,b)


def process_set(vector=False):
    p = Parameterization(air_density=.72,has_value=1)
    r = rate(0,[12.2,0,300.,300.,0.])
    reactions = [([0,2],[(1,1.),(4,2.4)],r),([1],[(0,1.),(3,1.4)],multiplier(rate(0),p,'qux')),
                 ([3],[],r),([2],[(1,1.),(3,2.5)],multiplier(rate(0),p,'qux'))]
    all_y = [[.1,.2,.3,.4,.5,0.],[1.1,1.2,1.3,1.4,1.5,0.]]
    all_k = [[10.,20.*70.*.72,30.,40.*70.*.72],[110.,120.*80.*.72,130.,140.*80.*.72]]
    pattern = [(0,0),(0,1),(0,2),(1,0),(1,1),(1,2),(2,0),(2,2),(3,1),(3,2),(3,3),(4,0),(4,2)]
    for group in ([1,2,3,4] if vector else [2]):
        for start in range(0,2,group):
            yy,kk = all_y[start:start+group],all_k[start:start+group]
            with Mechanism(['foo','bar','baz','quz','quuz','corge'],reactions,cells=len(yy)) as m:
                m.set_y(yy);set_rates(m,kk)
                assert [n for n in ['foo','bar','baz','quz','quuz','qux','corge'] if lib.chem_reaction_species_used(m.handle,n.encode())] == ['foo','bar','baz','quz','quuz','qux']
                assert [(i,j) for i in range(6) for j in range(6) if lib.chem_kinetic_nonzero(m.handle,i,j)] == pattern
                f,jac = forcing(m,1000.),jacobian(m,100.)
                for c,(y,k) in enumerate(zip(yy,kk)):
                    a,b,d,e = k
                    x,v,z,w,q,_ = y
                    expected = [1000-a*x*z+b*v,1000+a*x*z-b*v+e*z,1000-a*x*z-e*z,
                                1000+b*v*1.4-d*w+e*2.5*z,1000+a*x*z*2.4,1000.]
                    for i,value in enumerate(expected):
                        if start+c==1 and i==1:assert abs(f[c][i]-value) <= value*1e-12
                        else:real_eq(f[c][i],value)
                    values = [100+a*z,100-b,100+a*x,100-a*z,100+b,100-a*x-e,
                              100+a*z,100+a*x+e,100-1.4*b,100-2.5*e,100+d,100-2.4*a*z,100-2.4*a*x]
                    for (i,j),value in zip(pattern,values):real_eq(jac[c][i][j],value)


def random_processes():
    for nc,nr,n in [(200,50,40),(300,30,20),(400,100,80)]:
        counts,products,ids = uniform_integers(3),uniform_integers(10),uniform_integers(n-1)
        numbers = lognormals(4.)
        reactions = []
        for _ in range(nr):
            rr = [next(ids) for _ in range(next(counts))]
            pp = [(next(ids),1.2) for _ in range(next(products))]
            reactions.append((rr,pp,rate(0,[12.2,0,300,300,0])))
        with Mechanism([str(i) for i in range(n)],reactions,cells=nc) as m:
            m.set_y([[next(numbers) for _ in range(n)] for _ in range(nc)])
            set_rates(m,[[next(numbers) for _ in range(nr)] for _ in range(nc)])
            result = forcing(m,1000.)
            # The source random cases only require successful execution.
            assert all(math.isfinite(value) for row in result for value in row)


def algebraic():
    for group in [1,2,3,4]:
        all_y,all_k = [[1.,2.,3.,4.],[5.,6.,7.,8.]],[[10.],[20.]]
        for start in range(0,2,group):
            yy,kk = all_y[start:start+group],all_k[start:start+group]
            constraints = [linear(1,[(1,1.)],0),linear(3,[(3,1.)],0)]
            with Mechanism(['A','B','C','D'],[([0,1],[(2,1.),(3,1.)],rate(0))],constraints,cells=len(yy)) as m:
                m.set_y(yy);set_rates(m,kk)
                f,jac = forcing(m,1000.),jacobian(m,500.)
                for c,(y,k) in enumerate(zip(yy,kk)):
                    amount = k[0]*y[0]*y[1]
                    assert f[c] == [1000-amount,1000.,1000+amount,1000.]
                    for row in [0,1,2,3]:
                        for col in [0,1]:
                            derivative = k[0]*y[1-col]
                            expected = 500+derivative if row==0 else (500-derivative if row==2 else 500.)
                            real_eq(jac[c][row][col],expected)


def finite_difference_jacobian():
    with Mechanism(['A','B','C'],[([0,1],[(2,2.)],rate(0)),([2],[(0,1.)],rate(0))],cells=2) as m:
        yy = [[1.,2.,3.],[.5,4.,1.5]]
        m.set_y(yy);set_rates(m,[[5.,10.],[8.,12.]])
        analytical = jacobian(m,0.)
        fd = finite_difference(yy,lambda nc,n,y,f:lib.chem_add_forcing(m.state,y,f))
        pattern = [(i,j) for i in range(3) for j in range(3) if lib.chem_kinetic_nonzero(m.handle,i,j)]
        with Sparse(3,blocks=2,entries=pattern) as a:
            for c in range(2):
                for i,j in pattern:a[c,i,j] = analytical[c][i][j]
            comparison = Comparison()
            assert lib.chem_compare_jacobian(a.handle,fd,1e-7,1e-7,C.byref(comparison)) == 0
            assert comparison.passed
            assert lib.chem_check_jacobian_sparsity(a.handle,fd,1e-6,C.byref(comparison)) == 0
            assert comparison.passed


CASES = {
    'ProcessSet.Matrix':process_set,
    'ProcessSet.VectorMatrix':partial(process_set,True),
    'RandomProcessSet.Matrix':random_processes,
    'ProcessSetAlgebraicVariables.CudaMatrix':algebraic,
    'ProcessSetFiniteDifferenceJacobian.Matrix':finite_difference_jacobian,
}
