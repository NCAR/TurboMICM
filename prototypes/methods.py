"""The integrators, in Python, over the library's forcing, Jacobian, and
linear algebra (tmclib): the Rosenbrock-Wanner loop and the linearly
implicit Euler extrapolation, each a mirror of its C loop (src/row.c,
src/extrapolation.c) so that the agreement test can compare them step for
step. The prototypes measure accuracy and count operations; speed comes only
from the C programs.

    solver = RowSolver(lib, tableau, Options(rtol=1e-3, atol=1e-18, ...))
    result = solver.solve(state, y, dt)      # one solve interval, like tmc_solver_solve
    result.y, result.counts, result.status

A solver keeps, like a C state, the step it proposes for the next solve, the
factorization's age and step, and the controller's history.
"""
from dataclasses import dataclass, field

import numpy as np

H_MIN_FRACTION = 1.0e-15
H_START_FRACTION = 1.0e-6
SAME_STEP = 1.0e-12
ERROR_FLOOR = 1.0e-10


@dataclass
class Options:
    rtol: float = 1.0e-6
    atol: float = 1.0e-3
    warm_start: bool = True
    h_start: float = 0.0          # seconds (0: h_start_factor of the interval)
    h_start_factor: float = 0.0   # 0: MICM's 1e-6
    h_max: float = 0.0            # seconds (0: the interval)
    norm: str = "rms"             # rms (MICM), cellmax, max
    controller: str = "micm"      # micm, predictive, h211b
    clamp: bool = False
    quantize: bool = False
    reuse_age: int = 0            # steps a factorization may serve again (0: never)
    reuse_change: float = 0.1     # refactor when the rate constants changed by more than this
    freeze_band: float = 0.0
    max_steps: int = 1000
    order: int = 4                # extrapolation
    sequence: str = "harmonic"    # harmonic (1,2,3,..), shifted (2,3,4,..), even (2,4,6,..)
    safety: float = 0.0           # 0: MICM's 0.9
    fac_min: float = 0.0          # 0: MICM's 0.2
    fac_max: float = 0.0          # 0: MICM's 6
    h211_b: float = 0.0           # 0: 4
    h211_k: float = 0.0           # 0: the estimator order


@dataclass
class Counts:
    steps: int = 0          # attempts
    accepted: int = 0
    rejected: int = 0
    function_calls: int = 0
    jacobian_updates: int = 0
    decompositions: int = 0
    solves: int = 0
    matrix_reuses: int = 0

    def add(self, other):
        for k in self.__dataclass_fields__:
            setattr(self, k, getattr(self, k) + getattr(other, k))


@dataclass
class Result:
    y: np.ndarray
    counts: Counts
    status: str              # converged, nan, inf, too_small, max_steps
    final_time: float


def quantize(h, interval):
    if h >= interval:
        return interval
    n = max(1.0, np.ceil(interval / h - 1.0e-9))
    return interval / n


class Loop:
    """What the ROW and extrapolation loops share: the step-size rules."""

    def __init__(self, lib, options):
        self.lib, self.o = lib, options
        k = lib.loop_constants()
        self.factor_min = options.fac_min if options.fac_min > 0.0 else k.factor_min
        self.factor_max = options.fac_max if options.fac_max > 0.0 else k.factor_max
        self.rejection_factor_decrease = k.rejection_factor_decrease
        self.safety = options.safety if options.safety > 0.0 else k.safety_factor
        self.round_off = k.round_off
        self.h_next = 0.0
        self.lu_h = 0.0
        self.lu_age = 0
        self.err_prev = 0.0
        self.h_prev = 0.0
        self.order = 1.0  # the estimator order, set by the method

    def initial_step(self, interval, h_min, h_max):
        o = self.o
        if o.warm_start and self.h_next > 0.0:
            H = self.h_next
        elif o.h_start != 0.0:
            H = o.h_start
        else:
            H = (o.h_start_factor if o.h_start_factor > 0.0 else H_START_FRACTION) * interval
        H = min(max(h_min, abs(H)), h_max)
        if o.quantize:
            H = quantize(H, interval)
        return H

    def step_factor(self, error, H):
        o = self.o
        fac = self.safety / error ** (1.0 / self.order)
        history = self.err_prev > 0.0 and self.h_prev > 0.0
        if o.controller == "predictive" and history:
            predicted = self.safety * (H / self.h_prev) * (self.err_prev / (error * error)) ** (1.0 / self.order)
            fac = min(fac, predicted)
        elif o.controller == "h211b" and history:
            b = o.h211_b if o.h211_b > 0.0 else 4.0
            k = o.h211_k if o.h211_k > 0.0 else self.order
            fac = self.safety * (error * self.err_prev) ** (-1.0 / (b * k)) * (H / self.h_prev) ** (-1.0 / b)
        return min(max(fac, self.factor_min), self.factor_max)

    def error_norm(self, y, ynew, err, atol):
        scale = atol[None, :] + self.o.rtol * np.maximum(np.abs(y), np.abs(ynew))
        e = np.abs(err) / scale
        if self.o.norm == "rms":
            value = np.sqrt(np.mean(e * e))
        elif self.o.norm == "cellmax":
            value = np.max(np.sqrt(np.mean(e * e, axis=1)))
        else:
            value = np.max(e)
        if np.isnan(value):
            return value
        return max(value, ERROR_FLOOR)

    def accept(self, error, H, interval, h_min, h_max, reject_last, keep_matrix):
        """The new step after an accepted step; records the controller history."""
        o = self.o
        fac = self.step_factor(error, H)
        self.err_prev, self.h_prev = error, H
        Hnew = H * fac
        if reject_last:
            Hnew = min(Hnew, H)
        if keep_matrix and o.freeze_band > 1.0 and H / o.freeze_band < Hnew < H * o.freeze_band:
            Hnew = H
        Hnew = max(h_min, min(Hnew, h_max))
        if o.quantize:
            Hnew = quantize(Hnew, interval)
        return Hnew

    def reject(self, error, H, interval, reject_more):
        fac = self.step_factor(error, H)
        Hnew = H * fac
        if reject_more:
            Hnew = H * self.rejection_factor_decrease
        if self.o.quantize:
            Hnew = quantize(Hnew, interval)
        return Hnew


class RowSolver(Loop):
    def __init__(self, lib, tableau, options):
        super().__init__(lib, options)
        self.t = tableau
        self.order = tableau.order

    def stages(self, state, y, F0, H, counts):
        """K_i as tmc_rosenbrock_stages computes them; the matrix is factored."""
        t = self.t
        K = []
        raw_prev = None
        for stage in range(t.stages):
            base = stage * (stage - 1) // 2
            if stage == 0:
                raw = F0
            elif t.new_f[stage]:
                ys = y.copy()
                for j in range(stage):
                    ys = ys + t.a[base + j] * K[j]
                raw = state.forcing(ys)
                counts.function_calls += 1
            else:
                raw = raw_prev
            rhs = raw.copy()
            for j in range(stage):
                rhs = rhs + (t.c[base + j] / H) * K[j]
            K.append(state.solve(rhs))
            counts.solves += 1
            raw_prev = raw
        ynew = y.copy()
        yerr = np.zeros_like(y)
        for j in range(t.stages):
            ynew = ynew + t.m[j] * K[j]
            yerr = yerr + t.e[j] * K[j]
        return ynew, yerr

    def solve(self, state, y, interval):
        o, t = self.o, self.t
        counts = Counts()
        atol = state.absolute_tolerances()
        h_min = H_MIN_FRACTION * interval
        h_max = interval if o.h_max == 0.0 else min(interval, o.h_max)
        H = self.initial_step(interval, h_min, h_max)
        proposed = H
        time = 0.0
        reject_last = reject_more = False
        status = "running"
        y = np.array(y, dtype=np.float64)
        while time - interval + self.round_off <= 0.0 and status == "running":
            if counts.steps > o.max_steps:
                status = "max_steps"
                break
            if time + 0.1 * H == time or H <= min(self.round_off, h_min):
                status = "too_small"
                break
            H = min(H, abs(interval - time))
            F0 = state.forcing(y)
            counts.function_calls += 1
            reuse = (o.reuse_age > 0 and self.lu_h > 0.0 and abs(self.lu_h - H) <= SAME_STEP * H
                     and self.lu_age < o.reuse_age and state.rate_change() <= o.reuse_change)
            if reuse:
                counts.matrix_reuses += 1
            else:
                state.factor(y, 1.0 / (H * t.gamma))
                counts.jacobian_updates += 1
                counts.decompositions += 1
                self.lu_h, self.lu_age = H, 0
                state.reset_rate_change()
            factored = True
            while True:
                if not factored:
                    state.factor(y, 1.0 / (H * t.gamma))
                    counts.jacobian_updates += 1
                    counts.decompositions += 1
                    self.lu_h, self.lu_age = H, 0
                    state.reset_rate_change()
                factored = False
                ynew, yerr = self.stages(state, y, F0, H, counts)
                error = self.error_norm(y, ynew, yerr, atol)
                counts.steps += 1
                if np.isnan(error):
                    status = "nan"
                    break
                if np.isinf(error):
                    status = "inf"
                    break
                if error < 1.0 or H < h_min:
                    counts.accepted += 1
                    time += H
                    y = np.maximum(ynew, 0.0) if o.clamp else ynew
                    self.lu_age += 1
                    Hnew = self.accept(error, H, interval, h_min, h_max, reject_last, keep_matrix=True)
                    reject_last = reject_more = False
                    proposed = H = Hnew
                    break
                H = self.reject(error, H, interval, reject_more)
                reject_more, reject_last = reject_last, True
                counts.rejected += 1
                self.lu_h = 0.0
        if status == "running":
            status = "converged"
        self.h_next = proposed
        return Result(y, counts, status, time)


def sequence_entry(kind, j):
    if kind == "shifted":
        return j + 2
    if kind == "even":
        return 2 * (j + 1)
    return j + 1


class ExtrapolationSolver(Loop):
    def __init__(self, lib, options):
        super().__init__(lib, options)
        self.k = options.order
        self.order = float(self.k)

    def attempt(self, state, y, F0, H, counts):
        """One attempt: the tableau's rows for the sequence, the last two entries of the last row."""
        o, k = self.o, self.k
        prev = []
        for j in range(k):
            n_j = sequence_entry(o.sequence, j)
            state.factor(y, n_j / H)      # the C loop assembles J once and shifts; same numbers
            counts.decompositions += 1
            yj = y.copy()
            for m in range(n_j):
                rhs = F0 if m == 0 else state.forcing(yj)
                if m > 0:
                    counts.function_calls += 1
                yj = yj + state.solve(rhs)
                counts.solves += 1
            row = [yj]
            for l in range(1, j + 1):
                ratio = n_j / sequence_entry(o.sequence, j - l)
                row.append(row[l - 1] + (row[l - 1] - prev[l - 1]) / (ratio - 1.0))
            prev = row
        return prev[k - 1], prev[k - 1] - prev[k - 2]

    def solve(self, state, y, interval):
        o = self.o
        counts = Counts()
        atol = state.absolute_tolerances()
        h_min = H_MIN_FRACTION * interval
        h_max = interval if o.h_max == 0.0 else min(interval, o.h_max)
        H = self.initial_step(interval, h_min, h_max)
        proposed = H
        time = 0.0
        reject_last = reject_more = False
        status = "running"
        y = np.array(y, dtype=np.float64)
        while time - interval + self.round_off <= 0.0 and status == "running":
            if counts.steps > o.max_steps:
                status = "max_steps"
                break
            if time + 0.1 * H == time or H <= min(self.round_off, h_min):
                status = "too_small"
                break
            H = min(H, abs(interval - time))
            F0 = state.forcing(y)
            counts.function_calls += 1
            counts.jacobian_updates += 1
            while True:
                ynew, yerr = self.attempt(state, y, F0, H, counts)
                error = self.error_norm(y, ynew, yerr, atol)
                counts.steps += 1
                if np.isnan(error):
                    status = "nan"
                    break
                if np.isinf(error):
                    status = "inf"
                    break
                if error < 1.0 or H < h_min:
                    counts.accepted += 1
                    time += H
                    y = np.maximum(ynew, 0.0) if o.clamp else ynew
                    Hnew = self.accept(error, H, interval, h_min, h_max, reject_last, keep_matrix=False)
                    reject_last = reject_more = False
                    proposed = H = Hnew
                    break
                H = self.reject(error, H, interval, reject_more)
                reject_more, reject_last = reject_last, True
                counts.rejected += 1
        if status == "running":
            status = "converged"
        self.h_next = proposed
        self.lu_h = 0.0
        return Result(y, counts, status, time)


def make_solver(lib, method, options):
    if method == "extrap":
        return ExtrapolationSolver(lib, options)
    return RowSolver(lib, lib.tableau(method), options)
