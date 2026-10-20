#!/usr/bin/env python3
"""Exercise mechanism ordering, defaults, and custom parameter indexes."""
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TERMS = {"reactants": [{"name": "A"}], "products": [{"name": "B"}]}
reactions = [
    {"type": "user_defined", **TERMS,
     "parameters": {"label": "user", "scaling_factor": 3.0}},
    {"type": "surface", **TERMS, "surface_species": "A",
     "parameters": {"label": "aerosol", "reaction_probability": 0.74}},
    {"type": "arrhenius", **TERMS,
     "parameters": {"A": 1.25, "B": 2.0, "C": -4.0, "D": 300.0, "E": 1e-8}},
    {"type": "troe", **TERMS, "parameters": {}},
    {"type": "ternary_chemical_activation", **TERMS, "parameters": {}},
    {"type": "ternary_chemical_activation_jpl19", **TERMS, "parameters": {}},
    {"type": "branched_alkoxy", **TERMS,
     "parameters": {"X": 1.2, "Y": 204.3, "a0": 1e-3, "n": 2.0}},
    {"type": "tunneling", **TERMS, "parameters": {"A": 1.2, "B": 2.3, "C": 302.3}},
    {"type": "taylor_series", **TERMS,
     "parameters": {"A": 3e-11, "C": -200.0, "coefficients": [12.5, 0.013, 0.00052]}},
    {"type": "reversible", **TERMS,
     "parameters": {"A": 0.0114, "C": 2300.0, "k_r": 0.32}},
]
data = {
    "species": [
        {"name": "A", "molecular_weight": 0.025, "diffusion": 230.0},
        {"name": "B", "density_parameterized": True, "density_factor": 2.0},
    ],
    "reactions": reactions,
    "defaults": {"rate_parameters": {
        "user": 0.5,
        "aerosol.effective radius [m]": 1e-7,
        "aerosol.particle number concentration [# m-3]": 2.5e6,
    }},
}

with tempfile.TemporaryDirectory() as tmp:
    source = Path(tmp) / "mechanism.json"
    generated = Path(tmp) / "mechanism.c"
    source.write_text(json.dumps(data))
    subprocess.run([sys.executable, str(ROOT / "tools/generate_c.py"),
                    str(source), str(generated)], check=True)
    text = generated.read_text()

block = text.split("static const chem_reaction reactions[] = {\n", 1)[1].split("};", 1)[0]
rows = [line.strip() for line in block.splitlines() if line.strip()]
parsed = []
for line in rows:
    match = re.match(r'\{(CHEM_[A-Z0-9_]+), "([^"]*)".*\{([^}]*)\}, NULL\},', line)
    if not match:
        raise SystemExit(f"unrecognized generated reaction row: {line}")
    kind, label, values = match.groups()
    parsed.append((kind, label, [float(value.strip()) for value in values.split(",")]))

expected_order = [
    "CHEM_ARRHENIUS", "CHEM_TROE", "CHEM_TERNARY_ACTIVATION",
    "CHEM_TERNARY_ACTIVATION_JPL19", "CHEM_BRANCHED_ALKOXY",
    "CHEM_TUNNELING", "CHEM_TAYLOR_SERIES", "CHEM_REVERSIBLE",
    "CHEM_USER_DEFINED", "CHEM_SURFACE",
]
if [row[0] for row in parsed] != expected_order or len(parsed) != len(reactions):
    raise SystemExit("reaction ordering or total count differs from the rate type order")
species_block = text.split("static const chem_species species[] = {\n", 1)[1].split("};", 1)[0]
if not re.search(r'\{"B", 0, 0, 0, 0, 1, 2, NULL, 0, 0\.0\},', species_block):
    raise SystemExit("density parameter metadata was not preserved")
by_kind = {kind: values for kind, _, values in parsed}
if by_kind["CHEM_ARRHENIUS"][:5] != [1.25, 2.0, -4.0, 300.0, 1e-8]:
    raise SystemExit("Arrhenius parameters were not preserved")
if by_kind["CHEM_TROE"][:8] != [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.6, 1.0]:
    raise SystemExit("Troe parameter defaults differ")
if by_kind["CHEM_TERNARY_ACTIVATION"][:8] != [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.6, 1.0]:
    raise SystemExit("ternary activation parameter defaults differ")
if by_kind["CHEM_TERNARY_ACTIVATION_JPL19"][3:14] != [298.0, 1.0, 0.0, 0.0, 298.0,
                                                        1.0, 0.0, 0.0, 298.0, 0.6, 1.0]:
    raise SystemExit("JPL19 parameter defaults differ")
if by_kind["CHEM_BRANCHED_ALKOXY"][4] <= 0.0:
    raise SystemExit("branched derived parameter was not computed")
if by_kind["CHEM_USER_DEFINED"][0] != 3.0 or by_kind["CHEM_USER_DEFINED"][7] != 0.0:
    raise SystemExit("user-defined scaling or custom parameter index differs")
if by_kind["CHEM_SURFACE"][0] != 0.74 or by_kind["CHEM_SURFACE"][7] != 1.0:
    raise SystemExit("surface probability or custom parameter index differs")
labels = text.split("static const char *const custom_labels[] = {\n", 1)[1].split("};", 1)[0]
if re.findall(r'"([^"]+)"', labels) != [
        "user", "aerosol.effective radius [m]",
        "aerosol.particle number concentration [# m-3]"]:
    raise SystemExit("custom parameter labels or indexes differ")

# CAMP tutorial reader results must agree with the corresponding hand-built
# mechanisms, including the concentration-unit conversions done by MICM.
configs = ROOT.parent / "micm" / "test" / "tutorial" / "configs" / "rate_constants_no_user_defined"
with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    no_user_path = tmp / "no-user.json"
    user_path = tmp / "user.json"
    subprocess.run([sys.executable, str(ROOT / "tools/read_config.py"),
                    str(configs / "config.yaml"), str(no_user_path)], check=True)
    user_configs = ROOT.parent / "micm" / "test" / "tutorial" / "configs" / "rate_constants_user_defined"
    subprocess.run([sys.executable, str(ROOT / "tools/read_config.py"),
                    str(user_configs / "config.json"), str(user_path)], check=True)
    no_user = json.loads(no_user_path.read_text())
    user = json.loads(user_path.read_text())
    for source, output in ((no_user_path, tmp / "no-user.c"),
                           (user_path, tmp / "user.c")):
        subprocess.run([sys.executable, str(ROOT / "tools/generate_c.py"),
                        str(source), str(output)], check=True)
        generated_text = output.read_text()
        if generated_text.count("CHEM_") < 7:
            raise SystemExit("tutorial configuration did not generate C reaction definitions")

if no_user["species"] != user["species"] or no_user["reactions"] != user["reactions"][:7]:
    raise SystemExit("by-hand tutorial configurations do not share the same seven base reactions")
if len(no_user["species"]) != 7 or len(no_user["reactions"]) != 7 or len(user["reactions"]) != 10:
    raise SystemExit("rate-constants tutorial configurations have an unexpected mechanism size")
rate_conv = 1.0e-6 * 6.02214076e23
ternary = no_user["reactions"][4]["parameters"]
troe = no_user["reactions"][5]["parameters"]
if ternary["kinf_A"] != 2.6 / rate_conv:
    raise SystemExit("ternary concentration units do not match the hand-built tutorial")
if troe["k0_A"] != 1.8e-30 * rate_conv**2 or troe["kinf_A"] != 2.8e-11 * rate_conv:
    raise SystemExit("Troe concentration units do not match the hand-built tutorial")
if [reaction["type"] for reaction in user["reactions"][7:]] != [
        "USER_DEFINED", "USER_DEFINED", "USER_DEFINED"]:
    raise SystemExit("photolysis, emission, and loss entries were not normalized as custom rates")
