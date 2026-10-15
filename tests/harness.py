import math

CASES = {}


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
