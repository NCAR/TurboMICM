"""Stand-in for PyYAML used only for JSON inputs (JSON is a YAML subset).

The entries' tools import yaml unconditionally; no system python3 here has it.
"""
import json


def safe_load(stream):
    if hasattr(stream, "read"):
        stream = stream.read()
    return json.loads(stream)


def safe_dump(data, **kwargs):
    return json.dumps(data)
