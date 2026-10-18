"""Dense Python prototype of the adaptive ODE algorithms (not a speed test).

The caller supplies f, its true Jacobian, and a mathematical Rosenbrock table.
No C kinetics, Jacobian, factorization, or integration routine is used here.
The retained Tier 1 tables can be supplied by the test harness; RODAS3.1's new
table is independently transcribed below from the cited mathematical table.
"""
from dataclasses import dataclass
import math
import numpy as np


@dataclass
class Options:
    method: int = 43
    controller: int = 0
    norm: int = 1
    jacobian_steps: int = 1
    warm_start: int = 1
    max_steps: int = 100000
    newton_iterations: int = 8
    h_start: float = .01
    h_max: float = 0
    safety: float = .9
    factor_min: float = .2
    factor_max: float = 6
    pi_beta: float = .4
    lu_band: float = 1
    positivity: float = .1
    newton_tolerance: float = .03
    filter_b: float = 4
    filter_k: float = 0


def rodas31():
    return dict(stages=4, gamma=.515, order=3, new=[True,False,True,True],
                a=[0,.646601929740551,.409567801987914,
                   .646601929740551,.409567801987914,1],
                c=[4.198495621784201,3.711590161613010,-1.787771994729384,
                   4.458898153216104,-2.024095448516552,-2.626700600119396],
                m=[.646601929740551,.409567801987914,1,1],e=[0,0,0,1])


COUNTERS = ('internal_steps','accepted','rejected','function_calls','jacobian_updates',
            'decompositions','solves','jacobian_reuses','factor_reuses','positivity_rejections',
            'nonfinite_rejections','clipped_values','nonlinear_iterations','nonlinear_failures')


class Integrator:
    def __init__(self, f, jac, y, options, table=None, rtol=1e-4, atol=1e-18):
        self.f, self.jac, self.y, self.o, self.table = f, jac, y.copy(), options, table
        self.rtol, self.atol, self.next_h = rtol, np.broadcast_to(atol,(y.shape[0],)), 0.

    def norm(self, v, candidate):
        scale = self.atol[:,None]+self.rtol*np.maximum(abs(self.y),abs(candidate))
        with np.errstate(all='ignore'):
            e = np.where(v == 0,0,abs(v)/scale)
        if not np.all(np.isfinite(e)) or not np.all(np.isfinite(candidate)):
            return math.inf
        if self.o.norm == 2:
            return float(np.max(e))
        # Match the specified cell-then-species summation order.
        sums = [sum(float(e[i,c])**2 for i in range(e.shape[0])) for c in range(e.shape[1])]
        return (math.sqrt(sum(sums)/e.size) if self.o.norm == 0 else
                max(math.sqrt(s/e.shape[0]) for s in sums))

    def force(self, y):
        self.counts['function_calls'] += 1
        return self.f(y)

    def solve(self, rhs):
        self.counts['solves'] += 1
        # Factorization is mathematically reused even though NumPy's small-case
        # dense solve refactors internally. Counts model the numerical algorithm,
        # not BLAS calls; no Python timing is used.
        return np.stack([np.linalg.solve(self.matrix[c],rhs[:,c])
                         for c in range(rhs.shape[1])],axis=1)

    def ros_trial(self,h):
        p = self.table
        stages, forcing = [], self.force(self.y)
        for i in range(p['stages']):
            offset = i*(i-1)//2
            if i and p['new'][i]:
                stage_y = self.y.copy()
                for j in range(i):
                    stage_y += p['a'][offset+j]*stages[j]
                forcing = self.force(stage_y)
            rhs = forcing.copy()
            for j in range(i):
                rhs += (p['c'][offset+j]/h)*stages[j]
            stages.append(self.solve(rhs))
        candidate, error = self.y.copy(), np.zeros_like(self.y)
        for i in range(p['stages']):
            candidate += p['m'][i]*stages[i]
        for i in range(p['stages']):
            error += p['e'][i]*stages[i]
        return candidate,error,True

    def sdirk_trial(self,h):
        gamma = 1-1/math.sqrt(2)
        alpha = 1/(gamma*h)
        candidate,base = self.y.copy(),self.y.copy()
        for stage in range(2):
            converged = False
            for _ in range(self.o.newton_iterations):
                delta = self.solve(self.force(candidate)-alpha*(candidate-base))
                self.counts['nonlinear_iterations'] += 1
                candidate += delta
                e = self.norm(delta,candidate)
                if not math.isfinite(e):
                    break
                if e <= self.o.newton_tolerance:
                    converged = True
                    break
            if not converged:
                self.counts['nonlinear_failures'] += 1
                return candidate,np.zeros_like(candidate),False
            if not stage:
                increment = candidate-self.y
                base = self.y+((1-gamma)/gamma)*increment
        return candidate,increment/gamma-(candidate-self.y),True

    def step(self,dt):
        self.counts = {name:0 for name in COUNTERS}
        if dt == 0:
            return self.y,self.counts
        o = self.o
        hmax = min(dt,o.h_max) if o.h_max > 0 else dt
        h = min(hmax,self.next_h if o.warm_start and self.next_h > 0 else o.h_start)
        time, factor_h, previous_error, previous_h = 0.,0.,1.,h
        age, have_factor, rejected = o.jacobian_steps,False,False
        gamma = 1-1/math.sqrt(2) if o.method == 102 else self.table['gamma']
        order = 2 if o.method == 102 else self.table['order']
        while time < dt:
            assert self.counts['internal_steps'] < o.max_steps
            proposed = h
            h = min(h,dt-time)
            assert time+h != time
            fresh = age >= o.jacobian_steps
            if fresh:
                negative_jacobian = -self.jac(self.y)
                self.counts['jacobian_updates'] += 1
                age,have_factor = 0,False
            else:
                self.counts['jacobian_reuses'] += 1
            if not have_factor or factor_h != h:
                self.matrix = negative_jacobian.copy()
                for c in range(self.y.shape[1]):
                    for i in range(self.y.shape[0]):
                        self.matrix[c,i,i] += 1/(gamma*h)
                self.counts['decompositions'] += 1
                have_factor,factor_h = True,h
            else:
                self.counts['factor_reuses'] += 1
            self.counts['internal_steps'] += 1
            candidate,error,ok = self.sdirk_trial(h) if o.method == 102 else self.ros_trial(h)
            e = self.norm(error,candidate) if ok else math.inf
            negative = ok and np.any(candidate < -o.positivity*self.atol[:,None])
            if not math.isfinite(e):
                self.counts['nonfinite_rejections'] += 1
            if negative:
                self.counts['positivity_rejections'] += 1
            safe_error = max(e,1e-12)
            factor = o.safety*math.pow(safe_error,-1/order)
            if o.controller == 1 and self.counts['accepted'] and not rejected and e <= 1:
                factor = o.safety*math.pow(safe_error,-(1-o.pi_beta)/order)*math.pow(previous_error,o.pi_beta/order)
            if o.controller == 2 and self.counts['accepted'] and not rejected and e <= 1:
                k = o.filter_k if o.filter_k > 0 else order
                factor = (math.pow(o.safety,2/(o.filter_b*k)) *
                          math.pow(safe_error*previous_error,-1/(o.filter_b*k)) *
                          math.pow(h/previous_h,-1/o.filter_b))
            factor = min(o.factor_max,max(o.factor_min,factor))
            if math.isfinite(e) and e <= 1 and not negative:
                self.counts['clipped_values'] += int(np.sum(candidate < 0))
                self.y = np.maximum(candidate,0)
                self.counts['accepted'] += 1
                age += 1
                time = dt if h >= dt-time else time+h
                if rejected:
                    factor = min(1.,factor)
                next_h = min(hmax,h*factor)
                if o.lu_band > 1 and age < o.jacobian_steps and 1/o.lu_band <= factor <= o.lu_band:
                    next_h = h
                self.next_h = max(next_h,proposed) if time == dt and h < proposed and e < .5 else next_h
                previous_h,h,rejected,previous_error = h,next_h,False,max(e,1e-12)
            else:
                self.counts['rejected'] += 1
                rejected = True
                if fresh:
                    h *= o.factor_min if negative or not math.isfinite(e) else min(.9,factor)
                age,have_factor = o.jacobian_steps,False
        return self.y,self.counts
