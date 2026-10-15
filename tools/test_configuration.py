"""Configuration validation, normalization and complete TS1-default checks."""
import copy
import csv
import io
import json
from pathlib import Path
import subprocess
import tempfile
import yaml
from configure import read_mechanism, emit_c, AVOGADRO_CM

root = Path(__file__).resolve().parents[1]
checks = 0
base = {"camp-data": [{"type": "CHEM_SPEC", "name": "A"},
                      {"type": "CHEM_SPEC", "name": "B"},
                      {"type": "MECHANISM", "reactions": [{"type": "ARRHENIUS", "A": 2,
                        "reactants": {"A": {"qty": 2}}, "products": {"B": {}}}]}]}
with tempfile.TemporaryDirectory(dir=root / "build", prefix="configuration-") as temp:
    directory = Path(temp)
    def read(document, suffix="json"):
        path = directory / ("fixture." + suffix)
        path.write_text(json.dumps(document) if suffix == "json" else yaml.safe_dump(document))
        return read_mechanism(path)
    def reject(document):
        global checks
        try:
            emit_c(read(document), "fixture")
        except (ValueError, KeyError, TypeError):
            checks += 1
        else:
            raise AssertionError(f"Invalid document accepted: {document}")
    assert read(base) == read(base, "yaml")
    assert read(base)["reactions"][0]["params"]["A"] == 2 * AVOGADRO_CM
    checks += 2
    for value in ({}, None, 1, "record"):
        reject({"camp-data": value})
    for value in ({}, [1], "species.json"):
        reject({"camp-files": value})
    for value in ("", None, 12, "A\0B"):
        doc = copy.deepcopy(base); doc["camp-data"][0]["name"] = value; reject(doc)
    doc = copy.deepcopy(base); doc["camp-data"].insert(0, doc["camp-data"][0]); reject(doc)
    for key, value in (("type", "UNKNOWN"), ("A", float("inf")), ("A", True),
                       ("rate units", "unknown"), ("reactants", {"missing": {}}),
                       ("reactants", {"A": {"qty": 1.5}}), ("reactants", {"A": {"qty": 0}}),
                       ("reactants", "A"), ("reactants", {"A": 3}), ("products", {"missing": {}})):
        doc = copy.deepcopy(base); doc["camp-data"][-1]["reactions"][0][key] = value; reject(doc)
    for value in ({"missing": 1}, {"A": float("nan")}):
        doc = copy.deepcopy(base); doc["camp-data"][-1]["initial state"] = value; reject(doc)
    doc = copy.deepcopy(base); doc["camp-data"][0]["parameterization"] = "unknown"; reject(doc)
    doc = copy.deepcopy(base); doc["camp-data"][0]["absolute tolerance"] = -1; reject(doc)
    for count in (0, 17):
        doc = copy.deepcopy(base); r = doc["camp-data"][-1]["reactions"][0]
        r.update(type="TAYLOR_SERIES", coefficients=[1] * count); reject(doc)
    for n in (-1, 1.5):
        doc = copy.deepcopy(base); r = doc["camp-data"][-1]["reactions"][0]
        r.update(type="BRANCHED", n=n, **{"alkoxy products": {"B": {}}, "nitrate products": {"B": {}}}); reject(doc)
    for kind in ("ARRHENIUS", "TROE", "TERNARY_CHEMICAL_ACTIVATION", "TERNARY_CHEMICAL_ACTIVATION_JPL19",
                 "TUNNELING", "TAYLOR_SERIES", "REVERSIBLE", "USER_DEFINED"):
        doc = copy.deepcopy(base); r = doc["camp-data"][-1]["reactions"][0]; r["type"] = kind
        assert len(read(doc)["reactions"]) == 1
        checks += 1
    doc = copy.deepcopy(base); r = doc["camp-data"][-1]["reactions"][0]
    r.update(type="BRANCHED", n=2, a0=.001, **{"alkoxy products": {"B": {}}, "nitrate products": {"A": {}}})
    branches = read(doc)["reactions"]
    assert [r["params"]["branch"] for r in branches] == [0, 1]
    checks += 1
    doc = copy.deepcopy(base); r = doc["camp-data"][-1]["reactions"][0]
    r.update(reactants=["A", "A"], products=[{"species": "B", "yield": .3}, {"species": "B", "yield": .7}])
    assert read(doc)["reactions"][0]["products"] == [("B", .3), ("B", .7)]
    checks += 1
    for kind in ("PHOTOLYSIS", "FIRST_ORDER_LOSS", "EMISSION"):
        doc = copy.deepcopy(base); r = doc["camp-data"][-1]["reactions"][0]
        r.update(type=kind, species="A", label="custom")
        assert read(doc)["reactions"][0]["type"] == "USER"
        checks += 1
    (directory / "cycle.json").write_text('{"camp-files":["cycle.json"]}')
    try:
        read_mechanism(directory / "cycle.json")
    except ValueError:
        checks += 1
    else:
        raise AssertionError("include cycle accepted")
    for invalid in ("bad-name", "", "x;return 0;"):
        try:
            emit_c(read(base), invalid)
        except ValueError:
            checks += 1
        else:
            raise AssertionError("invalid generated symbol accepted")

output = subprocess.run([root / "bin/test_config", "--defaults"], capture_output=True, text=True, check=True)
actual = list(csv.DictReader(io.StringIO(output.stdout)))
with (root / "tests/data/ts1_defaults.csv").open() as f:
    expected = list(csv.DictReader(f))
assert len(actual) == len(expected) == 1
assert actual[0].keys() == expected[0].keys()
assert all(float(actual[0][key]) == float(value) for key, value in expected[0].items())
print(f"Configuration validation: PASS ({checks} checks; {len(expected[0])} TS1 defaults identical)")
