"""Mechanism configurations and their compiled form.

Reads MICM's configuration format (CAMP-style YAML or JSON, as in MICM's
test/tutorial/configs) and compiles a mechanism the way the C library's
builder does: the same reaction order, species order, Jacobian pattern, LU
fill, and fingerprint. tools/generate.py turns the result into C.

Format notes (see DECISIONS.md for the extensions):
  - `camp-files` lists further files, relative to the file that names them.
  - `camp-data` holds CHEM_SPEC entries (species with properties) and
    MECHANISM entries (reactions).
  - Reactants and products are maps {name: {qty, yield}} as in CAMP, or
    ordered lists (an extension): ["A", "B"] and [{"species": "C", "yield": 0.5}].
  - Rate parameters are in molecule cm-3 units, converted as MICM's tutorial
    programs do by hand, unless the mechanism says
    "rate constant units": "mol m-3" (an extension used by TS1).
  - A species with "third body": true (an extension) is parameterized on air
    density and holds no state variable.
  - "initial concentration [mol m-3]" on a species, and top-level
    "initial conditions" and "initial rate parameters" (extensions) record TS1's
    initial state.
"""
from __future__ import annotations

import json
import struct
from dataclasses import dataclass, field
from pathlib import Path

AVOGADRO = 6.02214076e23
BOLTZMANN = 1.380649e-23
MOLES_M3_TO_MOLECULES_CM3 = 1.0e-6 * AVOGADRO

# Rate constant types in the C enum's order, which is also MICM's sort order.
RATE_TYPES = [
    "ARRHENIUS", "TROE", "TERNARY_CHEMICAL_ACTIVATION", "TERNARY_JPL19", "BRANCHED",
    "TUNNELING", "TAYLOR_SERIES", "REVERSIBLE", "USER_DEFINED", "SURFACE", "LAMBDA",
]
RATE_TYPE_INDEX = {t: i for i, t in enumerate(RATE_TYPES)}


@dataclass
class Species:
    name: str
    properties: dict = field(default_factory=dict)  # key -> float | str | int | bool
    third_body: bool = False
    initial_concentration: float | None = None

    @property
    def diffusion_coefficient(self):
        return self.properties.get("diffusion coefficient [m2 s-1]")

    @property
    def molecular_weight(self):
        return self.properties.get("molecular weight [kg mol-1]")


@dataclass
class Reaction:
    type: str                          # one of RATE_TYPES
    reactants: list                    # ordered species names, repeats allowed
    products: list                     # ordered (name, yield)
    params: dict = field(default_factory=dict)  # numeric parameters by C field name
    label: str = ""                    # user-defined / surface / lambda label
    scaling_factor: float = 1.0
    reaction_probability: float = 1.0
    branch: str = "ALKOXY"
    surface_species: str = ""


@dataclass
class Mechanism:
    name: str
    species: list                       # declaration order
    reactions: list                     # declaration order
    temperature: float | None = None
    pressure: float | None = None
    rate_parameters: dict = field(default_factory=dict)

    def species_by_name(self):
        return {s.name: s for s in self.species}


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------

def _load(path: Path):
    text = path.read_text()
    if path.suffix.lower() in (".yaml", ".yml"):
        import yaml
        return yaml.safe_load(text)
    return json.loads(text)


def _collect(path: Path, seen=None):
    """All camp-data entries reachable from a configuration file, in order."""
    path = path.resolve()
    seen = seen if seen is not None else set()
    if path in seen:
        return [], {}
    seen.add(path)
    doc = _load(path)
    entries, extras = [], {}
    for key, value in doc.items():
        if key == "camp-files":
            for sub in value:
                e, x = _collect(path.parent / sub, seen)
                entries += e
                extras.update(x)
        elif key == "camp-data":
            entries += value
        else:
            extras[key] = value
    return entries, extras


def _species_list(value):
    """Reactants or products from a CAMP map or an ordered list; returns [(name, qty_or_yield)]."""
    out = []
    if isinstance(value, dict):
        for name, props in value.items():
            props = props or {}
            qty = props.get("qty", 1)
            if "yield" in props:
                out.append((name, float(props["yield"])))
            else:
                out += [(name, 1.0)] * int(qty)
    elif isinstance(value, list):
        for item in value:
            if isinstance(item, str):
                out.append((item, 1.0))
            else:
                out.append((item["species"], float(item.get("yield", 1.0))))
    elif value is None:
        pass
    else:
        raise ValueError(f"bad species list: {value!r}")
    return out


def _conv(n_reactants, order_shift, si):
    """Unit conversion factor for a pre-exponential: (1e-6 N_A)^(n_reactants + shift)."""
    if si:
        return 1.0
    return MOLES_M3_TO_MOLECULES_CM3 ** (n_reactants + order_shift)


def _reaction(entry, species, si):
    """One reaction entry to a Reaction, with MICM's unit conversions."""
    rtype = entry["type"]
    g = entry.get
    if rtype in ("PHOTOLYSIS", "USER_DEFINED"):
        reactants = [n for n, _ in _species_list(g("reactants"))]
        products = _species_list(g("products"))
        return Reaction("USER_DEFINED", reactants, products, label=g("MUSICA name", ""),
                        scaling_factor=float(g("scaling factor", 1.0)))
    if rtype == "EMISSION":
        return Reaction("USER_DEFINED", [], [(entry["species"], 1.0)], label=g("MUSICA name", ""),
                        scaling_factor=float(g("scaling factor", 1.0)))
    if rtype == "FIRST_ORDER_LOSS":
        return Reaction("USER_DEFINED", [entry["species"]], [], label=g("MUSICA name", ""),
                        scaling_factor=float(g("scaling factor", 1.0)))
    if rtype == "SURFACE":
        reactant = g("gas-phase reactant") or _species_list(g("reactants"))[0][0]
        products = _species_list(g("gas-phase products", g("products")))
        return Reaction("SURFACE", [reactant], products, label=g("MUSICA name", ""),
                        reaction_probability=float(g("reaction probability", 1.0)), surface_species=reactant)
    reactants = [n for n, _ in _species_list(g("reactants"))]
    products = _species_list(g("products"))
    # reactants that hold state (not parameterized) set the reaction order
    n_state = sum(1 for n in reactants if not species[n].third_body)
    if rtype == "ARRHENIUS":
        c = float(g("C", 0.0))
        if "Ea" in entry:
            c = -float(entry["Ea"]) / BOLTZMANN
        p = {"A": float(g("A", 1.0)) * _conv(n_state, -1, si), "B": float(g("B", 0.0)), "C": c,
             "D": float(g("D", 300.0)), "E": float(g("E", 0.0))}
        return Reaction("ARRHENIUS", reactants, products, p)
    if rtype in ("TROE", "TERNARY_CHEMICAL_ACTIVATION"):
        shift = 0 if rtype == "TROE" else -1
        p = {"k0_A": float(g("k0_A", 1.0)) * _conv(n_state, shift, si), "k0_B": float(g("k0_B", 0.0)),
             "k0_C": float(g("k0_C", 0.0)), "kinf_A": float(g("kinf_A", 1.0)) * _conv(n_state, shift - 1, si),
             "kinf_B": float(g("kinf_B", 0.0)), "kinf_C": float(g("kinf_C", 0.0)),
             "Fc": float(g("Fc", 0.6)), "N": float(g("N", 1.0))}
        return Reaction(rtype, reactants, products, p)
    if rtype == "TERNARY_JPL19":
        p = {k: float(g(k, d)) for k, d in (("k0_A", 1.0), ("k0_B", 0.0), ("k0_C", 0.0), ("k0_D", 298.0),
                                             ("kinf_A", 1.0), ("kinf_B", 0.0), ("kinf_C", 0.0), ("kinf_D", 298.0),
                                             ("kint_A", 1.0), ("kint_B", 0.0), ("kint_C", 0.0), ("kint_D", 298.0),
                                             ("Fc", 0.6), ("N", 1.0))}
        p["k0_A"] *= _conv(n_state, -1, si)
        p["kinf_A"] *= _conv(n_state, -2, si)
        p["kint_A"] *= _conv(n_state, -1, si)
        return Reaction(rtype, reactants, products, p)
    if rtype == "BRANCHED":
        # one entry gives the alkoxy and the nitrate branches, two reactions
        raise ValueError("BRANCHED is expanded by the caller")
    if rtype == "TUNNELING":
        p = {"A": float(g("A", 1.0)) * _conv(n_state, -1, si), "B": float(g("B", 0.0)), "C": float(g("C", 0.0))}
        return Reaction(rtype, reactants, products, p)
    if rtype == "TAYLOR_SERIES":
        coefficients = [float(x) for x in g("coefficients", [1.0])]
        p = {"A": float(g("A", 1.0)) * _conv(n_state, -1, si), "B": float(g("B", 0.0)), "C": float(g("C", 0.0)),
             "D": float(g("D", 300.0)), "E": float(g("E", 0.0)), "coefficients": coefficients}
        return Reaction(rtype, reactants, products, p)
    if rtype == "REVERSIBLE":
        p = {"A": float(g("A", 1.0)) * _conv(n_state, -1, si), "C": float(g("C", 0.0)), "k_r": float(g("k_r", 0.0))}
        return Reaction(rtype, reactants, products, p)
    raise ValueError(f"unknown reaction type {rtype!r}")


def read(path) -> Mechanism:
    """Read a configuration file (YAML or JSON) into a Mechanism."""
    path = Path(path)
    entries, extras = _collect(path)
    species = []
    reactions = []
    name = path.parent.name
    si = False
    for entry in entries:
        if entry.get("type") == "CHEM_SPEC":
            props = {k: v for k, v in entry.items() if k not in ("name", "type", "third body", "initial concentration [mol m-3]")}
            s = Species(entry["name"], props, bool(entry.get("third body", False)),
                        entry.get("initial concentration [mol m-3]"))
            species.append(s)
    by_name = {s.name: s for s in species}
    for entry in entries:
        if entry.get("type") != "MECHANISM":
            continue
        name = entry.get("name", name)
        si = entry.get("rate constant units", "molecule cm-3") == "mol m-3"
        for r in entry.get("reactions", []):
            if r["type"] == "BRANCHED":
                reactants = [n for n, _ in _species_list(r.get("reactants"))]
                n_state = sum(1 for n in reactants if not by_name[n].third_body)
                p = {"X": float(r["X"]) * _conv(n_state, -1, si), "Y": float(r["Y"]), "a0": float(r["a0"]), "n": int(r["n"])}
                reactions.append(Reaction("BRANCHED", reactants, _species_list(r.get("alkoxy products")), dict(p), branch="ALKOXY"))
                reactions.append(Reaction("BRANCHED", reactants, _species_list(r.get("nitrate products")), dict(p), branch="NITRATE"))
            else:
                reactions.append(_reaction(r, by_name, si))
    for r in reactions:
        for n in r.reactants + [p for p, _ in r.products]:
            if n not in by_name:
                raise ValueError(f"reaction references unknown species {n!r}")
    conditions = extras.get("initial conditions", {}) or {}
    return Mechanism(name, species, reactions, conditions.get("temperature"), conditions.get("pressure"),
                     dict(extras.get("initial rate parameters", {}) or {}))


# ---------------------------------------------------------------------------
# Compilation, mirroring src/builder.c
# ---------------------------------------------------------------------------

M64 = (1 << 64) - 1


def markowitz(n, pattern):
    """MICM's DiagonalMarkowitzReorder on a set of (row, col) edges; returns perm[new] = old.
    Degrees of zero wrap around as MICM's std::size_t arithmetic does."""
    out = [set() for _ in range(n)]
    inc = [set() for _ in range(n)]
    for i, j in pattern:
        out[i].add(j)
        inc[j].add(i)
    row_deg = [len(out[v]) for v in range(n)]
    col_deg = [len(inc[v]) for v in range(n)]
    alive = [True] * n
    perm = []
    for _ in range(n):
        pivot, best = -1, 0
        for v in range(n):
            if not alive[v]:
                continue
            cost = (((row_deg[v] - 1) & M64) * ((col_deg[v] - 1) & M64)) & M64
            if pivot < 0 or cost < best:
                best, pivot = cost, v
        perm.append(pivot)
        alive[pivot] = False
        cols = [c for c in out[pivot] if c != pivot and alive[c]]
        ins = [i for i in inc[pivot] if i != pivot and alive[i]]
        for i in ins:
            for c in cols:
                if c not in out[i]:
                    out[i].add(c)
                    row_deg[i] += 1
                    inc[c].add(i)
                    col_deg[c] += 1
        for c in cols:
            inc[c].discard(pivot)
            col_deg[c] -= 1
        for i in ins:
            out[i].discard(pivot)
            row_deg[i] -= 1
    return perm


def lu_fill(n, pattern):
    """MICM's GetLUMatrix: the pattern closed under elimination fill."""
    p = set(pattern)
    rows = {i: set() for i in range(n)}
    cols = {j: set() for j in range(n)}
    for i, j in p:
        rows[i].add(j)
        cols[j].add(i)
    for i in range(n):
        right = [k for k in rows[i] if k > i]
        below = [j for j in cols[i] if j > i]
        for k in right:
            for j in below:
                if (j, k) not in p:
                    p.add((j, k))
                    rows[j].add(k)
                    cols[k].add(j)
    return p


@dataclass
class Compiled:
    mechanism: Mechanism
    reactions: list          # sorted by rate type (stable)
    variables: list          # state variable names in solver order
    index: dict              # name -> solver index
    processes: list          # per reaction: (reactant ids, [(product id, yield)])
    jobs: list               # per Jacobian job: (process, independent, dependent ids, [(product id, yield)])
    jac_elements: list       # sorted (row, col) of the Jacobian, diagonal included
    lu_elements: list        # sorted (row, col) of the LU pattern (CSR numbering)
    lu_index: dict           # (row, col) -> element id
    fingerprint: int
    custom_labels: list      # custom rate parameter labels in order


def compile_mechanism(m: Mechanism, reorder=True) -> Compiled:
    by_name = m.species_by_name()
    reactions = sorted(m.reactions, key=lambda r: RATE_TYPE_INDEX[r.type])  # Python's sort is stable
    declared = [s.name for s in m.species if not s.third_body]
    dindex = {s: i for i, s in enumerate(declared)}
    rids = [[dindex[n] for n in r.reactants if not by_name[n].third_body] for r in reactions]
    pids = [[(dindex[n], y) for n, y in r.products if not by_name[n].third_body] for r in reactions]
    n = len(declared)
    perm = list(range(n))
    if reorder:
        pattern = set()
        for rs, ps in zip(rids, pids):
            for ind in rs:
                for d in rs:
                    pattern.add((d, ind))
                for d, _ in ps:
                    pattern.add((d, ind))
        perm = markowitz(n, pattern)
    new_of = {old: k for k, old in enumerate(perm)}
    variables = [declared[old] for old in perm]
    index = {v: k for k, v in enumerate(variables)}
    processes = [([new_of[i] for i in rs], [(new_of[i], y) for i, y in ps]) for rs, ps in zip(rids, pids)]
    # Jacobian jobs: per process, per reactant occurrence, stably sorted by independent variable
    raw = []
    for r, (rs, ps) in enumerate(processes):
        for ind in rs:
            raw.append((ind, r))
    raw.sort(key=lambda t: t[0])
    jobs = []
    for ind, r in raw:
        rs, ps = processes[r]
        found = False
        dep = []
        for x in rs:
            if x == ind and not found:
                found = True
                continue
            dep.append(x)
        jobs.append((r, ind, dep, list(ps)))
    jac = set()
    for rs, ps in processes:
        for ind in rs:
            for d in rs:
                jac.add((d, ind))
            for d, _ in ps:
                jac.add((d, ind))
    for i in range(n):
        jac.add((i, i))
    lu = sorted(lu_fill(n, jac))
    lu_index = {e: k for k, e in enumerate(lu)}
    custom_labels = []
    for r in reactions:
        if r.type == "USER_DEFINED":
            custom_labels.append(r.label)
        elif r.type == "SURFACE":
            custom_labels += [f"{r.label}.effective radius [m]", f"{r.label}.particle number concentration [# m-3]"]
    return Compiled(m, reactions, variables, index, processes, jobs, sorted(jac), lu, lu_index,
                    fingerprint(variables, reactions, processes, lu), custom_labels)


def fingerprint(variables, reactions, processes, lu):
    """FNV-1a over the same byte sequence as src/builder.c."""
    h = 14695981039346656037

    def bytes_(b):
        nonlocal h
        for x in b:
            h ^= x
            h = (h * 1099511628211) & M64

    def text(s):
        bytes_(s.encode() + b"\0")

    def integer(v):
        text(str(v))

    def double(v):
        bytes_(struct.pack("<d", float(v)))

    integer(len(variables))
    for v in variables:
        text(v)
    integer(len(reactions))
    for r, (rs, ps) in zip(reactions, processes):
        integer(RATE_TYPE_INDEX[r.type])
        integer(len(rs))
        for x in rs:
            integer(x)
        integer(len(ps))
        for x, y in ps:
            integer(x)
            double(y)
    integer(len(lu))
    for i, j in lu:
        integer(i)
        integer(j)
    return h


if __name__ == "__main__":
    import sys
    for arg in sys.argv[1:]:
        mech = read(arg)
        c = compile_mechanism(mech)
        print(f"{arg}: {len(mech.species)} species ({len(c.variables)} variables), {len(mech.reactions)} reactions, "
              f"Jacobian nnz {len(c.jac_elements)}, LU nnz {len(c.lu_elements)}, fingerprint {c.fingerprint:016x}")
