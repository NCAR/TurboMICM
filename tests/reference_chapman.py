"""Independent test-only oracle for the supplied legacy Chapman regression.

Closed-form forcing and a small dense Gaussian elimination reproduce the
legacy three-stage test oracle; no production code imports this module.
Numerical parameters are retained from the source regression helper.
"""
import math


def rates(temperature,photos):
    return [2.15e-11*math.exp(110/temperature),3.3e-11*math.exp(55/temperature),
            8e-12*math.exp(-2060/temperature),6e-34*(temperature/300)**2.4]+list(photos)


def force(k,y):
    a=k[0]*y[4]*y[5];b=k[1]*y[5]*y[7];c=k[2]*y[6]*y[8];d=k[3]*y[0]*y[6]*y[7]
    p=k[4]*y[7];q=k[5]*y[8];r=k[6]*y[8]
    return [0.]*5+[q-a-b,2*p+r+a+b-c-d,-p+q+r+2*c-d,-q-r-c+d]


def jacobian(k,y):
    # Positive dF/dy, in the legacy species order.
    j=[[0.]*9 for _ in range(9)]
    a=k[0]*y[4];b=k[1]*y[7];c=k[2]*y[8];d=k[3]*y[0]*y[7];e=k[3]*y[0]*y[6]
    j[5][4]=-k[0]*y[5];j[6][4]=k[0]*y[5]
    j[5][5]=-a-b;j[6][5]=a+b
    j[6][0]=-k[3]*y[6]*y[7];j[7][0]=j[6][0];j[8][0]=-j[6][0]
    j[6][6]=-c-d;j[7][6]=2*c-d;j[8][6]=-c+d
    j[5][7]=-k[1]*y[5];j[6][7]=2*k[4]+k[1]*y[5]-e;j[7][7]=-k[4]-e;j[8][7]=e
    j[5][8]=k[5];j[6][8]=k[6]-k[2]*y[6];j[7][8]=k[5]+k[6]+2*k[2]*y[6];j[8][8]=-k[5]-k[6]-k[2]*y[6]
    return j


def factor(matrix):
    a=[list(row) for row in matrix]
    for pivot in range(9):
        for row in range(pivot+1,9):
            a[row][pivot]/=a[pivot][pivot]
        for col in range(pivot+1,9):
            for row in range(pivot+1,9):a[row][col]-=a[row][pivot]*a[pivot][col]
    return a


def linear_solve(a,b):
    x=list(b)
    for i in range(9):
        for j in range(i):x[i]-=a[i][j]*x[j]
    for i in range(8,-1,-1):
        for j in range(i+1,9):x[i]-=a[i][j]*x[j]
        x[i]/=a[i][i]
    return x


def solve(y,k,dt=500.):
    y=list(y);time=0.;h=1e-6;reject_last=reject_more=False
    gamma=.43586652150845899941601945119356
    c=[-1.0156171083877702091975600115545,4.0759956452537699824805835358067,9.2076794298330791242156818474003]
    m=[1.,6.169794704382824559255361568973,-.42772256543218573326238373806514]
    e=[.5,-2.9079558716805469821718236208017,.22354069897811569627360909276199]
    iterations=0
    while time-dt+2.220446049250313e-16<=0:
        h=min(h,abs(dt-time));f0=force(k,y)
        while True:
            iterations+=1;assert iterations<101
            j=jacobian(k,y);a=factor([[(1/(h*gamma) if i==l else 0)-j[i][l] for l in range(9)] for i in range(9)])
            k0=linear_solve(a,f0)
            f1=force(k,[v+w for v,w in zip(y,k0)])
            k1=linear_solve(a,[v+c[0]/h*w for v,w in zip(f1,k0)])
            k2=linear_solve(a,[v+c[1]/h*w+c[2]/h*z for v,w,z in zip(f1,k0,k1)])
            new=list(y);err=[0.]*9
            for stage,ks in enumerate([k0,k1,k2]):
                for i in range(9):new[i]+=m[stage]*ks[i];err[i]+=e[stage]*ks[i]
            error=max(math.sqrt(sum((err[i]/(1e-3+1e-4*max(abs(y[i]),abs(new[i]))))**2 for i in range(9))/9),1e-10)
            hnew=h*min(6.,max(.2,.9/error**(1/3)))
            if error<1:
                time+=h;y=new;hnew=min(hnew,dt)
                if reject_last:hnew=min(hnew,h)
                reject_last=reject_more=False;h=hnew;break
            if reject_more:hnew=.1*h
            reject_more=reject_last;reject_last=True;h=hnew
    return y
