#!/usr/bin/env python3
"""Run a copied entry generator on a CAMP configuration and check the emitted
reaction terms, by species name, against the configuration.

usage: python3 -I e5_check.py RUN CONFIG_DIR OUT_DIR

RUN is one of the nine entry names; the copied tools are in ../RUN/tools.
Writes the generated C to OUT_DIR and prints a report.  Independent of the
entries' own code: the configuration is parsed here with json only.
"""
import json
import re
import runpy
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORK = HERE.parent


# --------------------------------------------------------------------------
# Reference: the configuration as written (JSON only), species in input order,
# reactions expanded the way MICM's tutorials build them by hand.
def reference(config_dir):
    cfg = json.loads((config_dir / "config.json").read_text())
    entries = []
    for name in cfg.get("camp-files", []):
        entries.extend(json.loads((config_dir / name).read_text())["camp-data"])
    entries.extend(cfg.get("camp-data", []))
    species = [e["name"] for e in entries if e["type"] == "CHEM_SPEC"]
    reactions = []

    def terms(mapping, key, default):
        out = []
        for name, opt in (mapping or {}).items():
            opt = opt or {}
            if key == "qty":
                out.extend([name] * int(opt.get("qty", 1)))
            else:
                out.append((name, float(opt.get("yield", 1))))
        return out

    for e in entries:
        if e["type"] != "MECHANISM":
            continue
        for r in e["reactions"]:
            t = r["type"]
            if t == "BRANCHED":
                rs = terms(r.get("reactants"), "qty", 1)
                reactions.append((t + "_ALKOXY", rs, terms(r.get("alkoxy products"), "yield", 1)))
                reactions.append((t + "_NITRATE", rs, terms(r.get("nitrate products"), "yield", 1)))
            elif t == "SURFACE":
                reactions.append((t, [r["gas-phase reactant"]], terms(r.get("gas-phase products"), "yield", 1)))
            elif t == "EMISSION":
                reactions.append((t, [], [(r["species"], 1.0)]))
            elif t == "FIRST_ORDER_LOSS":
                reactions.append((t, [r["species"]], []))
            else:
                reactions.append((t, terms(r.get("reactants"), "qty", 1), terms(r.get("products"), "yield", 1)))
    return species, reactions


def signature(rs, ps):
    return (tuple(sorted(rs)), tuple(sorted((n, round(y, 12)) for n, y in ps)))


# --------------------------------------------------------------------------
# Running a copied tool as a script, with the yaml stand-in and the tools
# directory placed on sys.path explicitly (python3 -I adds neither).
def run_script(run, script, argv):
    tools = WORK / run / "tools"
    for p in (str(HERE), str(tools)):
        if p not in sys.path:
            sys.path.insert(0, p)
    sys.argv = [script] + [str(a) for a in argv]
    runpy.run_path(str(tools / script), run_name="__main__")


def generate(run, config_dir, out_dir):
    cfg = config_dir / "config.json"
    out = out_dir / f"{run}-{config_dir.name}.c"
    if run == "sol-high":
        run_script(run, "configure.py", [cfg, out, "--symbol", "x"])
    elif run == "sol-xhigh":
        run_script(run, "configure.py", [cfg, out, "--symbol", "x"])
    elif run == "sol-max":
        run_script(run, "configure.py", [cfg, out, "--name", "x"])
    elif run == "astra-high":
        run_script(run, "configure.py", [cfg, out, "--symbol", "x"])
    elif run == "astra-xhigh":
        run_script(run, "generate.py", [cfg, out, "--symbol", "x"])
    elif run == "astra-max":
        run_script(run, "configure.py", [cfg, out, "--name", "x"])
    elif run == "luna-max":
        run_script(run, "generate_c.py", [cfg, out])
    elif run == "luna-xhigh":
        run_script(run, "generate.py", [cfg, out])
    elif run == "luna-high":
        normalized = out_dir / f"{run}-{config_dir.name}.json"
        run_script(run, "read_config.py", [cfg, normalized])
        run_script(run, "generate_c.py", [normalized, out])
    else:
        raise SystemExit("unknown run " + run)
    return out.read_text()


# --------------------------------------------------------------------------
# Parsers of the emitted C.  Each returns (species_table, reactions) where
# species_table is the emitted species name list (in the emitted order) and
# reactions is a list of (type_text, reactant_indices, [(product_index, yield)],
# params_text).  kernel(text, species_table) returns extra findings for the
# two entries that emit kernels.
def nums(s):
    return [float(x) for x in re.findall(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?", s)]


def parse_sol_high(text):
    species = re.findall(r'"([^"]+)"', re.search(r"species_names\[\] = \{(.*?)\};", text).group(1))
    arrays = {m.group(1): m.group(2) for m in re.finditer(r"static const (?:int|ChemProduct) (\w+)\[\] = \{(.*?)\};", text)}
    reactions = []
    for m in re.finditer(r"^\{(\d+),(\d+),(reactants_\d+),(products_\d+),\{\.type=CHEM_(\w+),\.p=\{([^}]*)\}", text, re.M):
        nr, np_, ra, pa, typ, params = m.groups()
        nr, np_ = int(nr), int(np_)
        rs = [int(x) for x in re.findall(r"-?\d+", arrays[ra])][:nr] if nr else []
        ps = [(int(a), float(b)) for a, b in re.findall(r"\{(\d+),([^}]+)\}", arrays[pa])][:np_] if np_ else []
        reactions.append((typ, rs, ps, params))
    return species, reactions


def parse_sol_xhigh(text):
    species = re.findall(r'\{\.name="([^"]+)"', text)
    arrays = {m.group(1): m.group(2) for m in re.finditer(r"static const (?:size_t|ck_product) (\w+)\[\] = \{(.*?)\};", text)}
    reactions = []
    for m in re.finditer(r"\{\.nreactants=(\d+), \.reactants=(\S+), \.nproducts=(\d+), \.products=(\S+), \.rate=\{(.*?)\}\},", text):
        nr, ra, np_, pa, params = m.groups()
        nr, np_ = int(nr), int(np_)
        rs = [int(x) for x in re.findall(r"\d+", arrays[ra])] if nr else []
        ps = [(int(a), float(b)) for a, b in re.findall(r"\{(\d+),([^}]+)\}", arrays[pa])] if np_ else []
        typ = re.search(r"\.type=(\w+)", params).group(1)
        reactions.append((typ, rs, ps, params))
    return species, reactions


def parse_sol_max(text):
    species = re.findall(r'chem_model_species\(m, "([^"]+)"', text)
    reactions = []
    for m in re.finditer(r"ChemRate rate = chem_rate\((CHEM_\w+)\);(.*?)const size_t reactants\[\] = \{([^}]*)\}, products\[\] = \{([^}]*)\};\s*const double yields\[\] = \{([^}]*)\};\s*if \(chem_model_reaction\(m, reactants, (\d+), products, yields, (\d+), rate\)\)", text, re.S):
        typ, params, ra, pa, ya, nr, np_ = m.groups()
        nr, np_ = int(nr), int(np_)
        rs = [int(x) for x in re.findall(r"\d+", ra)][:nr] if nr else []
        pi = [int(x) for x in re.findall(r"\d+", pa)][:np_] if np_ else []
        ps = list(zip(pi, nums(ya)[:np_])) if np_ else []
        reactions.append((typ, rs, ps, params))
    return species, reactions


def kernel_sol_max(text, species_ids_per_reaction):
    """Kernel forcing blocks use permuted variable indices; recover the map
    variable -> species id from each block and require one consistent bijection."""
    blocks = re.findall(r"double rate = k\[(\d+)\*s\+c\];((?:\s*(?:rate \*= y\[\d+\*s\+c\];|out\[\d+\*s\+c\] [-+]= [^;]*;))+)", text)
    mapping, conflicts = {}, []
    forcing_blocks = [(int(ri), body) for ri, body in blocks if ri is not None]
    # The first len(reactions) blocks with distinct ri in increasing order belong to the force kernel
    # (the jacobian kernel repeats reaction indices per column); take only blocks before "_jacobian".
    jac_start = text.find("_jacobian_0(")
    force_text = text[:jac_start] if jac_start > 0 else text
    blocks = re.findall(r"double rate = k\[(\d+)\*s\+c\];((?:\s*(?:rate \*= y\[\d+\*s\+c\];|out\[\d+\*s\+c\] [-+]= [^;]*;))+)", force_text)
    for ri, body in blocks:
        ri = int(ri)
        rs_var = [int(v) for v in re.findall(r"rate \*= y\[(\d+)\*s\+c\]", body)]
        ps_var = [int(v) for v in re.findall(r"out\[(\d+)\*s\+c\] \+= ", body)]
        rs_id, ps_id = species_ids_per_reaction[ri]
        if len(rs_var) != len(rs_id) or len(ps_var) != len(ps_id):
            conflicts.append(f"reaction {ri}: kernel term counts {len(rs_var)}/{len(ps_var)} vs descriptor {len(rs_id)}/{len(ps_id)}")
            continue
        for v, s in list(zip(rs_var, rs_id)) + list(zip(ps_var, ps_id)):
            if mapping.setdefault(v, s) != s:
                conflicts.append(f"reaction {ri}: variable {v} maps to species {mapping[v]} and {s}")
    bijective = len(set(mapping.values())) == len(mapping)
    return mapping, bijective, conflicts, len(blocks)


def parse_astra_high(text):
    species = re.findall(r'\{\.name="([^"]+)"', text)
    reactions = []
    rate_lines = {int(m.group(1)): (m.group(2), m.group(3)) for m in re.finditer(r"reactions\[(\d+)\]\.rate = chem_rate_default\((CHEM_\w+)\);((?:\s*reactions\[\d+\]\.rate\.[^\n]*\n)*)", text)}
    rarr = {int(m.group(1)): [int(x) for x in re.findall(r"\d+", m.group(2))] for m in re.finditer(r"const size_t reactants_(\d+)\[\] = \{([^}]*)\};", text)}
    parr = {int(m.group(1)): [(int(a), float(b)) for a, b in re.findall(r"\{(\d+), ([^}]+)\}", m.group(2))] for m in re.finditer(r"const ChemProduct products_(\d+)\[\] = \{(.*?)\};", text)}
    for i in sorted(rate_lines):
        typ, params = rate_lines[i]
        reactions.append((typ, rarr.get(i, []), parr.get(i, []), params))
    return species, reactions


def parse_astra_xhigh(text):
    species = re.findall(r'\{\.name="([^"]+)"', text)
    rarr = {int(m.group(1)): m.group(2) for m in re.finditer(r"static const size_t r(\d+)\[\]=\{([^}]*)\};", text)}
    parr = {int(m.group(1)): m.group(2) for m in re.finditer(r"static const chem_product p(\d+)\[\]=\{(.*?)\};", text)}
    reactions = []
    for m in re.finditer(r"^\s*\{(\d+),(\d+),r(\d+),p(\d+),\{\.type=(CHEM_\w+),(.*?)\}\},", text, re.M):
        nr, np_, ri, pi, typ, params = m.groups()
        nr, np_ = int(nr), int(np_)
        rs = [int(x) for x in re.findall(r"\d+", rarr[int(ri)])][:nr] if nr else []
        ps = [(int(a), float(b)) for a, b in re.findall(r"\{(\d+),([^}]+)\}", parr[int(pi)])][:np_] if np_ else []
        reactions.append((typ, rs, ps, params))
    return species, reactions


def kernel_astra_xhigh(text, species_ids_per_reaction):
    start = text.find("static void forcing(")
    end = text.find("static void jacobian(")
    force_text = text[start:end] if start >= 0 and end > start else ""
    blocks = re.findall(r"double r=rates\[(\d+)\*cells\+c\];((?:\s*(?:r\*=y\[\d+\*cells\+c\];|out\[\d+\*cells\+c\][-+]=[^;]*;))+)", force_text)
    mapping, conflicts = {}, []
    for ri, body in blocks:
        ri = int(ri)
        rs_var = [int(v) for v in re.findall(r"r\*=y\[(\d+)\*cells\+c\]", body)]
        ps_var = [int(v) for v in re.findall(r"out\[(\d+)\*cells\+c\]\+=", body)]
        rs_id, ps_id = species_ids_per_reaction[ri]
        if len(rs_var) != len(rs_id) or len(ps_var) != len(ps_id):
            conflicts.append(f"reaction {ri}: kernel term counts {len(rs_var)}/{len(ps_var)} vs descriptor {len(rs_id)}/{len(ps_id)}")
            continue
        for v, s in list(zip(rs_var, rs_id)) + list(zip(ps_var, ps_id)):
            if mapping.setdefault(v, s) != s:
                conflicts.append(f"reaction {ri}: variable {v} maps to species {mapping[v]} and {s}")
    bijective = len(set(mapping.values())) == len(mapping)
    return mapping, bijective, conflicts, len(blocks)


def parse_astra_max(text):
    species = re.findall(r'ch_model_add_species\(m,"([^"]+)"', text)
    reactions = []
    for m in re.finditer(r"ch_model_add_reaction\(m,(\d+),(NULL|\(int\[\]\)\{[^}]*\}),(\d+),(NULL|\(int\[\]\)\{[^}]*\}),(NULL|\(double\[\]\)\{[^}]*\}),&\(ch_rate\)\{\.type=(\d+),\.p=\{([^}]*)\}", text):
        nr, ra, np_, pa, ya, typ, params = m.groups()
        rs = [int(x) for x in re.findall(r"\d+", ra)] if ra != "NULL" else []
        pi = [int(x) for x in re.findall(r"\d+", pa)] if pa != "NULL" else []
        ys = nums(ya) if ya != "NULL" else []
        reactions.append(("type" + typ, rs, list(zip(pi, ys)), params))
    return species, reactions


def parse_luna_max(text):
    species = re.findall(r'\{\.name="([^"]+)"', text)
    flat_r = [int(x) for x in re.findall(r"\d+", re.search(r"chem_reactant_data\[\] = \{([^}]*)\}", text).group(1))]
    flat_p = [int(x) for x in re.findall(r"\d+", re.search(r"chem_product_data\[\] = \{([^}]*)\}", text).group(1))]
    flat_y = nums(re.search(r"chem_yield_data\[\] = \{([^}]*)\}", text).group(1))
    reactions = []
    for m in re.finditer(r"\{\.kind=(CHEM_\w+), \.reactant_offset=(\d+), \.reactant_count=(\d+), \.product_offset=(\d+), \.product_count=(\d+),(.*?)\.p=\{([^}]*)\}\}", text):
        typ, ro, nr, po, np_, _, params = m.groups()
        ro, nr, po, np_ = map(int, (ro, nr, po, np_))
        rs = flat_r[ro:ro + nr]
        ps = list(zip(flat_p[po:po + np_], flat_y[po:po + np_]))
        reactions.append((typ, rs, ps, params))
    return species, reactions


def parse_luna_xhigh(text):
    species = re.findall(r'^\s*\{ "([^"]+)",', text, re.M)
    block = lambda name: re.search(name + r"\[\d+\] = \{(.*?)\};", text, re.S).group(1)
    flat_r = [(int(a), float(b)) for a, b in re.findall(r"\{ (\d+)u, ([^}]+) \}", block("generated_reactants"))]
    flat_p = [(int(a), float(b)) for a, b in re.findall(r"\{ (\d+)u, ([^}]+) \}", block("generated_products"))]
    reactions = []
    for m in re.finditer(r"^\s*\{ (\d+)u, (\d+)u, (\d+)u, (\d+)u, \{ \.kind = (CHEM_RATE_\w+), (.*?)\.p = \{([^}]*)\}", text, re.M):
        ro, nr, po, np_, typ, _, params = m.groups()
        ro, nr, po, np_ = map(int, (ro, nr, po, np_))
        rs = [s for s, _ in flat_r[ro:ro + nr]]
        ps = flat_p[po:po + np_]
        reactions.append((typ, rs, ps, params))
    return species, reactions


def parse_luna_high(text):
    sp_block = re.search(r"static const chem_species species\[\] = \{(.*?)\n\};", text, re.S).group(1)
    species = re.findall(r'^\s*\{"([^"]+)",', sp_block, re.M)
    terms = {}
    for m in re.finditer(r"static const chem_term (rxn_\d+_[rp])\[\] = \{(.*?)\};", text, re.S):
        terms[m.group(1)] = [(int(a), float(b)) for a, b in re.findall(r"\{(\d+), ([^}]+)\}", m.group(2))]
    reactions = []
    rx_block = re.search(r"static const chem_reaction reactions\[\] = \{(.*?)\n\};", text, re.S).group(1)
    for m in re.finditer(r'^\s*\{(CHEM_\w+), "([^"]*)", (rxn_\d+_r|NULL), (\d+), (rxn_\d+_p|NULL), (\d+), ([^,]+), \{([^}]*)\}, NULL\},', rx_block, re.M):
        typ, label, ra, nr, pa, np_, surf, params = m.groups()
        rs = [s for s, _ in terms[ra]] if ra != "NULL" else []
        ps = terms[pa] if pa != "NULL" else []
        reactions.append((typ, rs, ps, params))
    return species, reactions


PARSERS = {
    "sol-high": parse_sol_high, "sol-xhigh": parse_sol_xhigh, "sol-max": parse_sol_max,
    "astra-high": parse_astra_high, "astra-xhigh": parse_astra_xhigh, "astra-max": parse_astra_max,
    "luna-max": parse_luna_max, "luna-xhigh": parse_luna_xhigh, "luna-high": parse_luna_high,
}
KERNELS = {"sol-max": kernel_sol_max, "astra-xhigh": kernel_astra_xhigh}

# Where each entry stores the pre-exponential (A) of an Arrhenius row and the
# k0_A / kinf_A of a Troe row, for the units probe.
def units_probe(run, reactions, species):
    out = {}
    for typ, rs, ps, params in reactions:
        names = [species[i] for i in rs]
        if sorted(names) == ["A", "B"]:  # bimolecular Arrhenius
            out["Arrhenius A+B A"] = first_param(run, params, "A")
        if names == ["A"] and ("TROE" in typ.upper() or (run == "astra-max" and typ == "type1")):
            out["Troe A->B k0_A"] = first_param(run, params, "k0_A")
            out["Troe A->B kinf_A"] = first_param(run, params, "kinf_A")
    return out


def first_param(run, params, key):
    if run in ("sol-xhigh",):
        return float(re.search(r"\." + key + r"=([^,}]+)", params).group(1))
    if run == "sol-max":
        m = re.search(r"rate\." + key + r" = ([^;]+);", params)
        return float(m.group(1)) if m else None
    if run == "astra-high":
        m = re.search(r"\." + key + r" = ([^;]+);", params)
        return float(m.group(1)) if m else None
    if run == "astra-xhigh":
        m = re.search(r"\." + key + r"=([^,}]+)", params)
        return float(m.group(1)) if m else None
    values = nums(params)
    slot = {"sol-high": {"A": 0, "k0_A": 0, "kinf_A": 3}, "astra-max": {"A": 0, "k0_A": 0, "kinf_A": 3},
            "luna-max": {"A": 0, "k0_A": 0, "kinf_A": 3}, "luna-xhigh": {"A": 0, "k0_A": 5, "kinf_A": 8},
            "luna-high": {"A": 0, "k0_A": 0, "kinf_A": 3}}[run][key]
    return values[slot] if slot < len(values) else None


def main():
    run, config_dir, out_dir = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
    out_dir.mkdir(parents=True, exist_ok=True)
    text = generate(run, config_dir, out_dir)
    species, reactions = PARSERS[run](text)
    ref_species, ref_reactions = reference(config_dir)
    print(f"== {run} on {config_dir.name}: input species {ref_species}")
    print(f"   emitted species table {species}")
    if sorted(species) != sorted(ref_species):
        print("   !! emitted species set differs from the configuration")
    ref_sigs = Counter(signature(rs, ps) for _, rs, ps in ref_reactions)
    got = []
    bad = 0
    for typ, rs, ps, _ in reactions:
        names_r = [species[i] if i < len(species) else f"?{i}" for i in rs]
        names_p = [(species[i] if i < len(species) else f"?{i}", y) for i, y in ps]
        sig = signature(names_r, names_p)
        ok = ref_sigs.get(sig, 0) > 0
        if ok:
            ref_sigs[sig] -= 1
        else:
            bad += 1
        got.append((typ, names_r, names_p, ok))
    for typ, nr, np_, ok in got:
        mark = "ok " if ok else "BAD"
        print(f"   {mark} {typ:28s} {' + '.join(nr) or '(none)':14s} -> {', '.join(f'{n}x{y:g}' for n, y in np_) or '(none)'}")
    missing = [s for s, c in ref_sigs.items() if c > 0]
    print(f"   result: {len(reactions)} emitted reactions, {len(ref_reactions)} in configuration, {bad} with terms not in the configuration, {len(missing)} configuration reactions unmatched")
    if run in KERNELS:
        ids = {i: (rs, [p for p, _ in ps]) for i, (_, rs, ps, _) in enumerate(reactions)}
        mapping, bijective, conflicts, nblocks = KERNELS[run](text, ids)
        print(f"   kernel check: {nblocks} forcing blocks; variable->species map {dict(sorted(mapping.items()))}; bijective={bijective}; conflicts={conflicts or 'none'}")
    if config_dir.name == "units":
        print("   units probe:", units_probe(run, reactions, species))
    print()


if __name__ == "__main__":
    main()
