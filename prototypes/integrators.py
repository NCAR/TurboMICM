"""Independent small-system reference integrators; no Python speed claims.

The chemistry is supplied by ordinary Python functions and Jacobian functions.
Dense Gaussian elimination is implemented here, independently of the sparse C
library. Coefficient tables for the retained Rosenbrock formulas are data passed
by the tests; all stage arithmetic and adaptive decisions run in Python.
"""
import math

COUNTERS = ('internal_steps accepted rejected function_calls jacobian_updates '
            'decompositions solves nonlinear_iterations matrix_reuses '
            'positivity_rejections failed_iterations clipped_values').split()


def lu(matrix):
    a = [row[:] for row in matrix]
    for k in range(len(a)):
        inverse = 1 / a[k][k]
        for i in range(k + 1, len(a)):
            a[i][k] *= inverse
        for i in range(k + 1, len(a)):
            for j in range(k + 1, len(a)):
                a[i][j] -= a[i][k] * a[k][j]
    return a


def solve(a, b):
    x = b[:]
    for i in range(len(x)):
        for j in range(i):
            x[i] -= a[i][j] * x[j]
    for i in reversed(range(len(x))):
        for j in range(i + 1, len(x)):
            x[i] -= a[i][j] * x[j]
        x[i] /= a[i][i]
    return x


def norm(y, yn, e, rtol, atol, mode=1):
    z = [abs(v) / (atol + rtol * max(abs(a), abs(b))) for a, b, v in zip(y, yn, e)]
    if not all(math.isfinite(v) for v in z):
        return math.inf
    return max(z) if mode == 2 else math.sqrt(sum(v*v for v in z) / len(z))


def positive(y, yn, rtol, atol, allowance):
    return all(math.isfinite(b) and b >= -allowance * (atol + rtol * abs(a)) for a, b in zip(y, yn))


def axpy(y, a, x):
    return [v + a*w for v, w in zip(y, x)]


def rosenbrock(f, jac, y0, dt, table, p, rtol, atol):
    counts = dict.fromkeys(COUNTERS, 0)
    y = y0[:]; t = 0.; h = p['h_start']; previous = 1.
    fresh = True; rejected = False; age = 0; matrix_h = None; previous_h = 0.
    while t < dt:
        h = min(h, p.get('h_max', 0) or dt, dt-t)
        if fresh:
            force = f(y); counts['function_calls'] += 1
            if not age or age >= p['jacobian_age']:
                j = jac(y); counts['jacobian_updates'] += 1; age = 0; matrix_h = None
            fresh = False
        if not p['reuse'] or matrix_h != h:
            a = [[-v for v in row] for row in j]
            for i in range(len(y)):
                a[i][i] += 1 / (h * table.gamma[0])
            a = lu(a); counts['decompositions'] += 1; matrix_h = h
        else:
            counts['matrix_reuses'] += 1
        stages = []
        reused = None
        for k in range(table.stages):
            offset = k*(k-1)//2
            if k == 0:
                rhs = force[:]
            elif table.new_function[k]:
                v = y[:]
                for l in range(k):
                    v = axpy(v, table.a[offset+l], stages[l])
                rhs = f(v); counts['function_calls'] += 1
            else:
                rhs = reused[:]
            reused = rhs[:]
            for l in range(k):
                rhs = axpy(rhs, table.c[offset+l]/h, stages[l])
            stages.append(solve(a, rhs)); counts['solves'] += 1
        yn = y[:]; e = [0.] * len(y)
        for k in range(table.stages):
            yn = axpy(yn, table.m[k], stages[k]); e = axpy(e, table.e[k], stages[k])
        error = norm(y, yn, e, rtol, atol, p['norm'])
        good = positive(y, yn, rtol, atol, p['positivity'])
        counts['internal_steps'] += 1
        if p['controller'] == 2 and not rejected and previous_h > 0:
            exponent = 1/(p['controller_b']*p['controller_k'])
            factor = p['safety'] * max(error, 1e-12)**(-exponent) * previous**(-exponent) * (h/previous_h)**(-1/p['controller_b'])
        elif p['controller'] == 1 and not rejected:
            factor = p['safety'] * max(error, 1e-12)**(-.7/table.order) * previous**(.4/table.order)
        else:
            factor = p['safety'] * max(error, 1e-12)**(-1/table.order)
        factor = min(p['growth'], max(p['shrink'], factor))
        if error <= 1 and good:
            counts['clipped_values'] += sum(v < 0 for v in yn)
            y = [max(0., v) for v in yn]; t = dt if h == dt-t else t+h
            counts['accepted'] += 1; age += 1; fresh = True; previous = max(error, 1e-12)
            if rejected:
                factor = min(1., factor)
            if 1/p.get('step_hold', 1) <= factor <= p.get('step_hold', 1):
                factor = 1.
            previous_h = h; h *= factor; rejected = False
        else:
            counts['rejected'] += 1
            if not good:
                counts['positivity_rejections'] += 1; factor = min(.5, factor)
            h *= min(.8, factor); rejected = True
            if age:
                age = 0; fresh = True
            matrix_h = None
    return y, counts


def bdf_coefficients(q, times, h):
    x = [0.] + [(v-h)/h for v in times[:q]]
    d = [-sum(1/v for v in x[1:])]
    for j in range(1, q+1):
        value = 1/x[j]
        for k in range(1, q+1):
            if k != j:
                value *= -x[k]/(x[j]-x[k])
        d.append(value)
    return d


def extrapolate(history, times, points, h):
    result = [0.] * len(history[0])
    for j in range(points):
        weight = 1.
        for k in range(points):
            if k != j:
                weight *= (h-times[k])/(times[j]-times[k])
        result = axpy(result, weight, history[j])
    return result


def bdf(f, jac, y0, dt, p, rtol, atol):
    counts = dict.fromkeys(COUNTERS, 0)
    y = y0[:]; history = [y[:]]; times = [0.]; order = 1
    t = 0.; h = p['h_start']; successes = 0; rejected = False; factor_alpha = None; age = 0

    def estimator(q, yn):
        poly = extrapolate(history, times, q+1, h)
        d0 = bdf_coefficients(q, times, h)[0]
        fac = 1 / (d0 * (1-times[q]/h))
        e = [fac*(v-w) for v, w in zip(yn, poly)]
        return norm(y, yn, e, rtol, atol, p['norm'])

    def matrix(alpha):
        a = [[-v for v in row] for row in jac(y)]
        counts['jacobian_updates'] += 1
        for i in range(len(y)):
            a[i][i] += alpha
        counts['decompositions'] += 1
        return lu(a)

    while t < dt:
        h = min(h, p.get('h_max', 0) or dt, dt-t)
        q = min(order, p['max_order'], len(history))
        d = bdf_coefficients(q, times, h); alpha = d[0]/h
        constant = [0.] * len(y)
        for k in range(1, q+1):
            constant = axpy(constant, -d[k]/h, history[k-1])
        predictor = y[:] if len(history) == 1 else extrapolate(history, times, min(q+1, len(history)), h)
        yn = [max(0., v) for v in predictor]
        updated = not p['reuse'] or factor_alpha is None or age >= p['jacobian_age'] or abs(alpha/factor_alpha-1) > p['reuse_ratio']
        if updated:
            a = matrix(alpha); factor_alpha = alpha; age = 0
        else:
            counts['matrix_reuses'] += 1
        converged = False; last = math.inf
        for retry in range(2):
            for iteration in range(8):
                force = f(yn); counts['function_calls'] += 1
                delta = [v + (c-alpha*w) for v, c, w in zip(force, constant, yn)]
                delta = solve(a, delta); counts['solves'] += 1; counts['nonlinear_iterations'] += 1
                err = norm(y, yn, delta, rtol, atol, p['norm'])
                if not math.isfinite(err):
                    break
                yn = axpy(yn, 1., delta)
                if err < .03:
                    converged = True; break
                if iteration and err > last*.9:
                    break
                last = err
            if converged or updated:
                break
            a = matrix(alpha); factor_alpha = alpha; age = 0; updated = True
            last = math.inf; yn = [max(0., v) for v in predictor]
        counts['internal_steps'] += 1
        error = math.inf
        if converged:
            if len(history) == 1:
                force = f(y); counts['function_calls'] += 1
                e = [.5*(v-w-h*z) for v, w, z in zip(yn, y, force)]
                error = norm(y, yn, e, rtol, atol, p['norm'])
            else:
                error = estimator(q, yn)
        else:
            counts['failed_iterations'] += 1; factor_alpha = None
        good = positive(y, yn, rtol, atol, p['positivity'])
        if converged and good and error <= 1:
            counts['clipped_values'] += sum(v < 0 for v in yn)
            counts['accepted'] += 1; age += 1; successes += 1
            y = [max(0., v) for v in yn]
            history = [y[:]] + history[:6]; times = [0.] + [v-h for v in times[:6]]
            t = dt if h == dt-t else t+h
            ratio = p['safety'] * max(error, 1e-12)**(-1/(q+1)); next_order = q
            if successes >= q+1 and not rejected:
                if q > 1:
                    pred = extrapolate(history, times, q+1, h)
                    low = estimator(q-1, pred)
                    other = p['safety'] * max(low, 1e-12)**(-1/q)
                    if other > ratio*1.2:
                        ratio = other; next_order = q-1
                if q < p['max_order'] and len(history) >= q+3:
                    pred = extrapolate(history, times, q+3, h)
                    high = estimator(q+1, pred)
                    other = p['safety'] * max(high, 1e-12)**(-1/(q+2))
                    if other > ratio*1.1:
                        ratio = other; next_order = q+1
                if next_order != q:
                    successes = 0
            ratio = min(p['growth'], max(p['shrink'], ratio))
            if rejected:
                ratio = min(1., ratio)
            if p['reuse'] and .8 <= ratio <= 1.3:
                ratio = 1.
            order = next_order; h *= ratio; rejected = False
        else:
            counts['rejected'] += 1
            if not good:
                counts['positivity_rejections'] += 1
            ratio = p['safety'] * max(error, 1e-12)**(-1/(q+1)) if math.isfinite(error) else .25
            h *= max(p['shrink'], min(.5, ratio)); rejected = True; successes = 0
            if q > 1:
                order = q-1
    return y, counts
