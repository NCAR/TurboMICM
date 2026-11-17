"""Read the CAMP-style mechanism files used by MICM's tutorial examples."""
from pathlib import Path
import json
import yaml


def read_document(path):
    path = Path(path)
    with path.open(encoding="utf-8") as stream:
        if path.suffix.lower() == ".json":
            return json.load(stream)
        if path.suffix.lower() in (".yaml", ".yml"):
            return yaml.safe_load(stream)
    raise ValueError(f"unsupported configuration format: {path}")


def read_mechanism(path):
    path = Path(path)
    root = read_document(path)
    species = []
    mechanisms = []
    for name in root["camp-files"]:
        for entry in read_document(path.parent / name)["camp-data"]:
            if entry["type"] == "CHEM_SPEC":
                species.append(entry)
            elif entry["type"] == "MECHANISM":
                mechanisms.append(entry)
            else:
                raise ValueError(f"unsupported CAMP entry: {entry['type']}")
    names = [entry["name"] for entry in species]
    if len(names) != len(set(names)):
        raise ValueError("duplicate species name")
    if len(mechanisms) != 1:
        raise ValueError("expected exactly one mechanism")
    return species, mechanisms[0]["reactions"]
