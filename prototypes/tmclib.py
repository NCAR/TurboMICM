"""The C library's pieces for the prototypes, through ctypes (see shim.c).

    lib = Library()                      # build/libtmc_proto.{dylib,so}, from `make proto`
    st = lib.state(n_cells)              # a TS1 state with InitState in every cell
    st.set_conditions(cell, T, P); st.set_variable(cell, v, value); st.set_custom(cell, k, value)
    st.update_rate_constants()           # every rate constant
    st.update_custom_rate_constants()    # only those reading custom parameters (photolysis)
    f = st.forcing(y)                    # F(y), y and f (n_cells, n_vars) arrays
    st.factor(y, alpha); x = st.solve(b) # alpha I - J(y) factored in the state; (alpha I - J)^-1 b
    J = st.dense_jacobian(cell, y)       # dF/dy of one cell, dense
    lib.tableau("rodas3")                # a Rosenbrock parameter set in MICM's form

Row-major (cell, variable) numpy arrays everywhere; the shim converts to the
state's vector ordering.
"""
import ctypes
import os
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ENTRY = HERE.parent


def _library_path():
    for name in ("libtmc_proto.dylib", "libtmc_proto.so"):
        path = ENTRY / "build" / name
        if path.exists():
            return path
    sys.exit("prototypes: build the shared library first: make proto")


@dataclass
class Tableau:
    name: str
    stages: int
    a: np.ndarray      # [15] stage combination coefficients, MICM's packed order
    c: np.ndarray      # [15]
    m: np.ndarray      # [6] solution weights
    e: np.ndarray      # [6] error weights
    gamma: float       # the diagonal gamma
    new_f: list        # per stage: a new forcing evaluation, or the previous stage's
    order: float       # MICM's estimator_of_local_order


@dataclass
class LoopConstants:
    factor_min: float
    factor_max: float
    rejection_factor_decrease: float
    safety_factor: float
    max_steps: int
    round_off: float


class Library:
    def __init__(self, generic=False):
        self._lib = ctypes.CDLL(str(_library_path()))
        L = self._lib
        D, I, P, S = ctypes.c_double, ctypes.c_int, ctypes.c_void_p, ctypes.c_char_p
        dp = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
        L.proto_open.restype = P
        L.proto_open.argtypes = [I]
        L.proto_close.argtypes = [P]
        for name in ("proto_n_vars", "proto_n_custom", "proto_n_reactions", "proto_uses_kernels"):
            getattr(L, name).restype = I
            getattr(L, name).argtypes = [P]
        L.proto_variable_name.restype = S
        L.proto_variable_name.argtypes = [P, I]
        L.proto_custom_label.restype = S
        L.proto_custom_label.argtypes = [P, I]
        L.proto_state_create.restype = P
        L.proto_state_create.argtypes = [P, I]
        L.proto_state_free.argtypes = [P]
        L.proto_set_conditions.argtypes = [P, I, D, D]
        L.proto_set_variable.argtypes = [P, I, I, D]
        L.proto_get_variable.restype = D
        L.proto_get_variable.argtypes = [P, I, I]
        L.proto_set_custom.argtypes = [P, I, I, D]
        L.proto_get_custom.restype = D
        L.proto_get_custom.argtypes = [P, I, I]
        L.proto_rate_constant.restype = D
        L.proto_rate_constant.argtypes = [P, I, I]
        L.proto_set_tolerances.argtypes = [P, D, D]
        L.proto_absolute_tolerance.restype = D
        L.proto_absolute_tolerance.argtypes = [P, I]
        L.proto_update_rate_constants.argtypes = [P]
        L.proto_update_custom_rate_constants.argtypes = [P]
        L.proto_rate_change.restype = D
        L.proto_rate_change.argtypes = [P]
        L.proto_reset_rate_change.argtypes = [P]
        L.proto_get_variables.argtypes = [P, dp]
        L.proto_set_variables.argtypes = [P, dp]
        L.proto_forcing.argtypes = [P, dp, dp]
        L.proto_factor.argtypes = [P, dp, D]
        L.proto_solve.argtypes = [P, dp]
        L.proto_dense_jacobian.argtypes = [P, I, dp, dp]
        L.proto_factor_slot.argtypes = [P, I, dp, D]
        L.proto_solve_slot.argtypes = [P, I, dp]
        L.proto_tableau.restype = I
        L.proto_tableau.argtypes = [S, ctypes.POINTER(I), dp, dp, dp, dp, dp, np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"), ctypes.POINTER(D)]
        L.proto_loop_constants.argtypes = [ctypes.POINTER(D)] * 4 + [ctypes.POINTER(I), ctypes.POINTER(D)]
        self._solver = L.proto_open(1 if generic else 0)
        if not self._solver:
            sys.exit("prototypes: the TS1 solver could not be built")
        self.n_vars = L.proto_n_vars(self._solver)
        self.n_custom = L.proto_n_custom(self._solver)
        self.n_reactions = L.proto_n_reactions(self._solver)
        self.uses_kernels = bool(L.proto_uses_kernels(self._solver))
        self.variable_names = [L.proto_variable_name(self._solver, i).decode() for i in range(self.n_vars)]
        self.custom_labels = [L.proto_custom_label(self._solver, i).decode() for i in range(self.n_custom)]
        self.variable_index = {name: i for i, name in enumerate(self.variable_names)}
        self.custom_index = {label: i for i, label in enumerate(self.custom_labels)}
        self.photolysis = [i for i, label in enumerate(self.custom_labels) if label.startswith("j")]

    def __del__(self):
        if getattr(self, "_solver", None):
            self._lib.proto_close(self._solver)
            self._solver = None

    def state(self, n_cells):
        return State(self, n_cells)

    def tableau(self, name):
        stages = ctypes.c_int()
        order = ctypes.c_double()
        a, c = np.zeros(15), np.zeros(15)
        m, e, gamma = np.zeros(6), np.zeros(6), np.zeros(1)
        new_f = np.zeros(6, dtype=np.int32)
        if not self._lib.proto_tableau(name.encode(), ctypes.byref(stages), a, c, m, e, gamma, new_f, ctypes.byref(order)):
            raise KeyError(f"no parameter set named {name}")
        return Tableau(name, stages.value, a, c, m, e, float(gamma[0]), [bool(x) for x in new_f[: stages.value]], order.value)

    def loop_constants(self):
        values = [ctypes.c_double() for _ in range(4)]
        steps = ctypes.c_int()
        round_off = ctypes.c_double()
        self._lib.proto_loop_constants(*[ctypes.byref(v) for v in values], ctypes.byref(steps), ctypes.byref(round_off))
        return LoopConstants(*[v.value for v in values], steps.value, round_off.value)


class State:
    def __init__(self, lib, n_cells):
        self.lib, self.n_cells, self.n_vars = lib, n_cells, lib.n_vars
        self._state = lib._lib.proto_state_create(lib._solver, n_cells)

    def __del__(self):
        if getattr(self, "_state", None):
            self.lib._lib.proto_state_free(self._state)
            self._state = None

    def set_conditions(self, cell, temperature, pressure):
        self.lib._lib.proto_set_conditions(self._state, cell, temperature, pressure)

    def set_variable(self, cell, v, value):
        self.lib._lib.proto_set_variable(self._state, cell, v, value)

    def set_custom(self, cell, k, value):
        self.lib._lib.proto_set_custom(self._state, cell, k, value)

    def get_custom(self, cell, k):
        return self.lib._lib.proto_get_custom(self._state, cell, k)

    def rate_constant(self, cell, r):
        return self.lib._lib.proto_rate_constant(self._state, cell, r)

    def set_tolerances(self, rtol, atol):
        self.lib._lib.proto_set_tolerances(self._state, rtol, atol)

    def absolute_tolerances(self):
        return np.array([self.lib._lib.proto_absolute_tolerance(self._state, v) for v in range(self.n_vars)])

    def update_rate_constants(self):
        self.lib._lib.proto_update_rate_constants(self._state)

    def update_custom_rate_constants(self):
        self.lib._lib.proto_update_custom_rate_constants(self._state)

    def rate_change(self):
        return self.lib._lib.proto_rate_change(self._state)

    def reset_rate_change(self):
        self.lib._lib.proto_reset_rate_change(self._state)

    def variables(self):
        out = np.zeros((self.n_cells, self.n_vars))
        self.lib._lib.proto_get_variables(self._state, out)
        return out

    def set_variables(self, y):
        self.lib._lib.proto_set_variables(self._state, np.ascontiguousarray(y, dtype=np.float64))

    def forcing(self, y):
        f = np.zeros((self.n_cells, self.n_vars))
        self.lib._lib.proto_forcing(self._state, np.ascontiguousarray(y, dtype=np.float64), f)
        return f

    def factor(self, y, alpha):
        self.lib._lib.proto_factor(self._state, np.ascontiguousarray(y, dtype=np.float64), alpha)

    def solve(self, b):
        x = np.array(b, dtype=np.float64, order="C", copy=True)
        self.lib._lib.proto_solve(self._state, x)
        return x

    def factor_slot(self, slot, y, alpha):
        self.lib._lib.proto_factor_slot(self._state, slot, np.ascontiguousarray(y, dtype=np.float64), alpha)

    def solve_slot(self, slot, b):
        x = np.array(b, dtype=np.float64, order="C", copy=True)
        self.lib._lib.proto_solve_slot(self._state, slot, x)
        return x

    def dense_jacobian(self, cell, y):
        J = np.zeros((self.n_vars, self.n_vars))
        self.lib._lib.proto_dense_jacobian(self._state, cell, np.ascontiguousarray(y, dtype=np.float64), J)
        return J
