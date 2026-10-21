"""Stand-in for PyYAML, so the copied generators can import yaml.

Only JSON inputs are used in this scratch run, and JSON is read with json.loads.
"""
import json


def safe_load(stream):
    text = stream.read() if hasattr(stream, "read") else stream
    return json.loads(text)
