"""Independent Chapman test oracle: explicit chemistry and dense stage solves.

The three-stage method and constants follow the supplied legacy reference;
this test implementation calls no library kinetics, Jacobian, or solver code.
"""
import math
import numpy as np

NAMES = ['M','Ar','CO2','H2O','N2','O1D','O','O2','O3']
POSITIONS = [(0,0),(6,0),(7,0),(8,0),(1,1),(2,2),(3,3),(4,4),(5,4),(6,4),(5,5),(6,5),
             (6,6),(7,6),(8,6),(5,7),(6,7),(7,7),(8,7),(5,8),(6,8),(7,8),(8,8)]


def rates(t,photo):
    return np.array([2.15e-11*math.exp(110/t),3.3e-11*math.exp(55/t),8e-12*math.exp(-2060/t),
                     6e-34*math.pow(t/300,2.4),*photo])


def forcing(k,y):
    f = np.zeros(9)
    f[5] = k[5]*y[8]-k[0]*y[4]*y[5]-k[1]*y[5]*y[7]
    f[6] = 2*k[4]*y[7]+k[6]*y[8]+k[0]*y[4]*y[5]+k[1]*y[5]*y[7]-k[2]*y[6]*y[8]-k[3]*y[0]*y[6]*y[7]
    f[7] = -k[4]*y[7]+k[5]*y[8]+k[6]*y[8]+2*k[2]*y[6]*y[8]-k[3]*y[0]*y[6]*y[7]
    f[8] = -k[5]*y[8]-k[6]*y[8]-k[2]*y[6]*y[8]+k[3]*y[0]*y[6]*y[7]
    return f


def jacobian(k,y):
    a = np.array([0,-k[3]*y[6]*y[7],-k[3]*y[6]*y[7],k[3]*y[6]*y[7],0,0,0,0,
        -k[0]*y[5],k[0]*y[5],-k[0]*y[4]-k[1]*y[7],k[0]*y[4]+k[1]*y[7],
        -k[2]*y[8]-k[3]*y[0]*y[7],2*k[2]*y[8]-k[3]*y[0]*y[7],-k[2]*y[8]+k[3]*y[0]*y[7],
        -k[1]*y[5],2*k[4]+k[1]*y[5]-k[3]*y[0]*y[6],-k[4]-k[3]*y[0]*y[6],k[3]*y[0]*y[6],
        k[5],k[6]-k[2]*y[6],k[5]+k[6]+2*k[2]*y[6],-k[5]-k[6]-k[2]*y[6]])
    result = np.zeros((9,9))
    for (i,j),v in zip(POSITIONS,a): result[i,j] = v
    return result


def solve(y,k,dt):
    y = np.array(y,copy=True); time,h = 0.,1e-6; reject_last = reject_more = False; steps = 0
    gamma = .43586652150845899941601945119356
    c = [-1.0156171083877702091975600115545,4.0759956452537699824805835358067,9.2076794298330791242156818474003]
    m = [1.,6.1697947043828245592553615689730,-.42772256543218573326238373806514]
    e = [.5,-2.9079558716805469821718236208017,.22354069897811569627360909276199]
    while time < dt:
        assert steps <= 100 and h > np.finfo(float).eps
        h = min(h,dt-time); f0 = forcing(k,y)
        while True:
            matrix = np.eye(9)/(h*gamma)-jacobian(k,y)
            k0 = np.linalg.solve(matrix,f0)
            f1 = forcing(k,y+k0)
            k1 = np.linalg.solve(matrix,f1+c[0]/h*k0)
            k2 = np.linalg.solve(matrix,f1+c[1]/h*k0+c[2]/h*k1)
            new = y.copy(); error = np.zeros(9)
            for a,b,stage in zip(m,e,(k0,k1,k2)): new += a*stage; error += b*stage
            norm = max(np.sqrt(np.mean((error/(.001+.0001*np.maximum(abs(y),abs(new))))**2)),1e-10)
            proposed = h*min(6,max(.2,.9/math.pow(norm,1/3))); steps += 1
            if norm < 1:
                time += h; y = new
                h = min(proposed,dt,h if reject_last else dt); reject_last = reject_more = False
                break
            if reject_more: proposed = h*.1
            reject_more,reject_last,h = reject_last,True,proposed
    return y
