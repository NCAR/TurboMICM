import math
import struct

CASES = {}


class NotApplicable(Exception):
    pass


def case(*identifiers):
    def register(function):
        for identifier in identifiers:
            if identifier in CASES:
                raise ValueError(f'duplicate case: {identifier}')
            CASES[identifier] = function
        return function
    return register


def near(actual, expected, tolerance):
    if not math.isfinite(actual) or abs(actual-expected) > tolerance:
        raise AssertionError(f'{actual:.17g} != {expected:.17g} within {tolerance:.6g}')


def equal(actual, expected):
    if actual != expected:
        raise AssertionError(f'{actual!r} != {expected!r}')


def real_equal(actual,expected):
    def ordered(value):
        bits = struct.unpack('>Q',struct.pack('>d',float(value)))[0]
        return (1 << 63)-bits if bits >> 63 else bits+(1 << 63)
    if not math.isfinite(actual) or abs(ordered(actual)-ordered(expected)) > 4:
        raise AssertionError(f'{actual:.17g} != {expected:.17g} within 4 ULP')
