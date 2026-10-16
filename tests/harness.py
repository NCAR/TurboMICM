"""Explicit registry: missing cases fail, and N/A needs a case-specific reason."""
from pathlib import Path

CASE_LIST = [line.split()[0] for line in (Path(__file__).parent/'cases.txt').read_text().splitlines()
             if line and not line.startswith('#')]
REGISTRY = {}
NOT_APPLICABLE = {}


def case(*ids):
    def register(function):
        for id_ in ids:
            assert id_ in CASE_LIST, id_
            assert id_ not in REGISTRY, id_
            REGISTRY[id_] = function
        return function
    return register


def group(prefix, names):
    return case(*(prefix+name for name in names.split()))


def na(id_,reason):
    assert id_ in CASE_LIST and id_ not in REGISTRY
    NOT_APPLICABLE[id_] = reason
