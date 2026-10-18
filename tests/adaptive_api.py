"""Bindings for opt-in integrators; inherited ctypes layouts stay unchanged."""
from api import C, D, I, Z, P, SP, Stats, lib, bind


class Adaptive(C.Structure):
    _fields_ = [(k, I) for k in ('method', 'controller', 'norm', 'max_order')] + [
        ('jacobian_age', Z), ('max_steps', Z)] + [(k, D) for k in (
        'h_start h_max safety growth shrink positivity reuse_ratio').split()] + [('reuse', I)] + [
        (k, D) for k in ['step_hold', 'controller_b', 'controller_k']]


class AdaptiveStats(C.Structure):
    _fields_ = [('base', Stats)] + [(k, C.c_uint64) for k in (
        'nonlinear_iterations matrix_reuses positivity_rejections failed_iterations clipped_values').split()] + [('largest_clip', D)]


bind('chem_adaptive_parameters', Adaptive, I)
bind('chem_workspace_create', P, SP)
bind('chem_workspace_free', None, P)
bind('chem_workspace_reset', None, P)
bind('chem_adaptive_solve', I, C.POINTER(Adaptive), D, SP, P, C.POINTER(AdaptiveStats))


def counts(stats):
    result = {k: getattr(stats.base, k) for k, _ in Stats._fields_ if k != 'final_time'}
    result.pop('constraint_init_iterations')
    result.update({k: getattr(stats, k) for k, _ in AdaptiveStats._fields_ if k not in ['base', 'largest_clip']})
    return result
