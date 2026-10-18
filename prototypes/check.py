#!/usr/bin/env python3
"""Independent dense Python integrator vs C on a reversible nonlinear dimer.

Pure Python, no C solver pieces used. C prints its retained MICM coefficient
Tables, which are shared data, not shared stage/controller/linear algebra code.
The Tier 1 table tests independently check these data. No Python timing is used.
"""
import itertools
import math
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[1]
INITIAL=[1.,.5,.1,.2,.01,.03]

def forcing(y):
    out=[0.]*6
    for c in range(2):
        a,b,d=y[c],y[2+c],y[4+c]
        # Preserve individual reaction accumulation, including duplicate terms.
        r=20*a*a
        out[c]-=r;out[c]-=r;out[2+c]+=r
        r=2*b
        out[2+c]-=r;out[c]+=2*r
        r=.4*b
        out[2+c]-=r;out[4+c]+=r
        r=.1*d
        out[4+c]-=r;out[c]+=2*r
    return out

def jacobian(y):
    return [[[80*y[c],-4.,-.2],[-40*y[c],2.4,0.],[0.,-.4,.1]] for c in range(2)]

def factor(jac,alpha):
    lu=[[row[:] for row in cell] for cell in jac]
    for a in lu:
        for i in range(3): a[i][i]+=alpha
        for k in range(3):
            inv=1/a[k][k]
            for i in range(k+1,3):
                a[i][k]*=inv
                for j in range(k+1,3): a[i][j]-=a[i][k]*a[k][j]
    return lu

def solve(lu,rhs):
    x=rhs[:]
    for c,a in enumerate(lu):
        for i in range(3):
            for j in range(i): x[2*i+c]-=a[i][j]*x[2*j+c]
        for i in range(2,-1,-1):
            for j in range(i+1,3): x[2*i+c]-=a[i][j]*x[2*j+c]
            x[2*i+c]/=a[i][i]
    return x

def add(y,a,x):
    for i in range(6): y[i]+=a*x[i]

def norm(y,z,err,kind):
    total=worst=negative=0.
    for c in range(2):
        ss=0.
        for j in range(3):
            i=2*j+c
            scale=1e-11+1e-5*max(abs(y[i]),abs(z[i]))
            e=err[i]/scale
            ss+=e*e
            if kind==2: worst=max(worst,abs(e))
            negative=max(negative,-z[i]/scale)
        total+=ss
        if kind==1: worst=max(worst,math.sqrt(ss/3))
    if kind==0: worst=math.sqrt(total/6)
    return max(negative,worst,1e-12)

def integrate(method,ctl,jage,deadband,kind,warm,p):
    y=INITIAL[:];proposed=0.;counts=[0]*7;outputs=[]
    for call in range(2):
        time=0.;h=min(.2,proposed if warm and proposed>0 else .01)
        previous=1.;previous_h=0.;factored_h=0.;age=jage;rejected=False;fresh=False;accepted=0;rhs_valid=False
        while time<.2:
            h=min(h,.2-time)
            if not rhs_valid:
                f=forcing(y);counts[3]+=1;rhs_valid=True
            if age>=jage:
                jac=jacobian(y);counts[4]+=1;age=0;fresh=True;factored_h=0.
            if method in (12,13):
                order=method-10;table=[]
                for k in range(order):
                    steps=k+1;lu=factor(jac,steps/h);counts[5]+=1
                    z=y[:]
                    for j in range(steps):
                        rhs=f if j==0 else forcing(z)
                        if j: counts[3]+=1
                        add(z,1.,solve(lu,rhs));counts[6]+=1
                    for j in range(k):
                        old=z[:];ratio=steps/(k-j)-1
                        z=[z[i]+(old[i]-table[j][i])/ratio for i in range(6)]
                        table[j]=old
                    if k==order-1: err=[z[i]-table[k-1][i] for i in range(6)]
                    table.append(z[:])
            else:
                stages,order,gamma=p['parameters']
                if h!=factored_h: lu=factor(jac,1/(h*gamma));counts[5]+=1
                stage=[None]*int(stages)
                for k in range(int(stages)):
                    off=k*(k-1)//2
                    if k==0: stage[k]=f[:]
                    elif p['new'][k]:
                        z=y[:]
                        for j in range(k): add(z,p['a'][off+j],stage[j])
                        stage[k]=forcing(z);counts[3]+=1
                    if k+1<stages and not p['new'][k+1]: stage[k+1]=stage[k][:]
                    for j in range(k): add(stage[k],p['c'][off+j]/h,stage[j])
                    stage[k]=solve(lu,stage[k]);counts[6]+=1
                z=y[:];err=[0.]*6
                for k in range(int(stages)):
                    add(z,p['m'][k],stage[k]);add(err,p['e'][k],stage[k])
            factored_h=h
            e=norm(y,z,err,kind);counts[0]+=1
            fac=.9*e**(-1/order)
            if ctl==1 and not rejected and accepted: fac=.9*e**(-.7/order)*previous**(.4/order)
            if ctl==2 and not rejected and accepted:
                fac=.9*e**(-1/(4*order))*previous**(-1/(4*order))*(h/previous_h)**(-1/4)
            fac=min(6.,max(.1,fac))
            if e<=1:
                counts[1]+=1;accepted+=1
                time=.2 if h>=.2-time else time+h
                y=[max(0.,v) for v in z]
                rhs_valid=False
                age+=1;fresh=False
                if rejected: fac=min(1.,fac)
                previous=max(e,1e-12);previous_h=h;rejected=False
                if abs(fac-1)<=deadband: fac=1.
                proposed=h*fac;h=proposed
            else:
                counts[2]+=1
                if not fresh: age=jage
                else: h*=min(.5,fac)
                rejected=True
        outputs.append(y[:])
    return outputs,counts

def main():
    total=0;worst_error=0.
    for method,(ctl,jage,deadband,kind,warm) in itertools.product([2,3,4,43,6,12,13,31],
        [(0,1,0.,1,1),(1,1,0.,1,1),(2,1,0.,1,1),(0,4,.2,1,1),(0,1,0.,0,0),(0,2,.5,2,1)]):
        text=subprocess.check_output([str(ROOT/'bin/integrator_probe'),*map(str,[method,ctl,jage,deadband,kind,warm])],text=True)
        p={};ys=[]
        for line in text.splitlines():
            key,*values=line.split();v=list(map(float,values))
            if key=='y': ys.append(v)
            else: p[key]=v
        expected,counts=integrate(method,ctl,jage,deadband,kind,warm,p)
        assert counts==p['counts'],(method,ctl,jage,kind,counts,p['counts'])
        error=max(abs(a-b) for x,z in zip(expected,ys) for a,b in zip(x,z))
        worst_error=max(worst_error,error)
        assert error<3e-13,(method,ctl,jage,kind,error)
        uncached=subprocess.check_output([str(ROOT/'bin/integrator_probe'),
            *map(str,[method,ctl,jage,deadband,kind,warm,0])],text=True)
        lines=uncached.splitlines()
        uncached_y=[list(map(float,line.split()[1:])) for line in lines if line.startswith('y ')]
        uncached_counts=[int(x) for line in lines if line.startswith('counts ') for x in line.split()[1:]]
        assert uncached_y==ys, 'Forcing reuse changed the C trajectory'
        counts[3]+=counts[2]
        assert counts==uncached_counts, 'Forcing reuse must save exactly one call per rejection'
        total+=1
    print(f'{total} C/prototype comparisons PASS: eight methods, exact operation counts, solutions within 3e-13')
    print(f'Maximum C/prototype absolute difference: {worst_error:.3g}')
    print(f'{total} forcing-cache checks PASS: identical trajectories, one evaluation saved per rejection')
if __name__=='__main__': main()
