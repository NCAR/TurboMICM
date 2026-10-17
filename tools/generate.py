#!/usr/bin/env python3
"""Generate the C for a mechanism: its species and reactions as library calls,
its recorded initial state, and specialized kernels (forcing, Jacobian, sparse
LU, triangular solves) with the species and matrix indices as literals.

    generate.py CONFIG OUTPUT_STEM --name NAME [--width W] [--chunk N]
                [--no-kernels] [--reciprocal] [--table-forcing] [--table-jacobian]

writes OUTPUT_STEM.c and OUTPUT_STEM.h. The kernels compute exactly what the
library's generic kernels compute, in the same order (see DESIGN.md), for
states whose vector width is W; the library checks the width and a
fingerprint of the mechanism when the kernels are attached.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mechanism as mech  # noqa: E402


def c_double(x):
    """A C literal that round-trips the double."""
    r = repr(float(x))
    if r in ("inf", "-inf", "nan"):
        raise ValueError(f"non-finite constant {x}")
    if "e" not in r and "." not in r:
        r += ".0"
    return r


def c_string(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


class Emitter:
    def __init__(self):
        self.lines = []

    def __call__(self, line=""):
        self.lines.append(line)

    def text(self):
        return "\n".join(self.lines) + "\n"


# ---------------------------------------------------------------------------
# Mechanism tables: species, reactions, initial state
# ---------------------------------------------------------------------------

def emit_mechanism(out, name, m: mech.Mechanism, c: mech.Compiled):
    phase = "gas"
    out(f"tmc_system {name}_system(void)")
    out("{")
    out(f"  tmc_phase phase = tmc_phase_create({c_string(phase)});")
    out("  tmc_species s;")
    for s in m.species:
        out(f"  s = tmc_species_create({c_string(s.name)});")
        for key, value in s.properties.items():
            if isinstance(value, bool):
                out(f"  tmc_species_set_bool(&s, {c_string(key)}, {'true' if value else 'false'});")
            elif isinstance(value, int):
                out(f"  tmc_species_set_int(&s, {c_string(key)}, {value});")
            elif isinstance(value, float):
                out(f"  tmc_species_set_double(&s, {c_string(key)}, {c_double(value)});")
            elif isinstance(value, str):
                out(f"  tmc_species_set_string(&s, {c_string(key)}, {c_string(value)});")
        if s.third_body:
            out("  tmc_species_set_third_body(&s);")
        out("  tmc_phase_add_species(&phase, &s);")
    out("  tmc_system system = tmc_system_create(&phase);")
    out("  tmc_phase_free(&phase);")
    out("  return system;")
    out("}")
    out()
    out(f"int {name}_reactions(tmc_reaction **reactions)")
    out("{")
    out(f"  tmc_system system = {name}_system();")
    out("  const tmc_phase *phase = &system.gas_phase;")
    out(f"  tmc_reaction *list = malloc({max(len(m.reactions), 1)} * sizeof *list);")
    out("  tmc_reaction *r = list;")
    for rx in m.reactions:  # declared order; the builder sorts
        out("  *r = tmc_reaction_create();")
        for sp in rx.reactants:
            out(f"  tmc_reaction_add_reactant(r, &tmc_phase_find(phase, {c_string(sp)})->species);")
        for sp, y in rx.products:
            out(f"  tmc_reaction_add_product(r, &tmc_phase_find(phase, {c_string(sp)})->species, {c_double(y)});")
        p = rx.params
        if rx.type == "ARRHENIUS":
            out(f"  tmc_reaction_set_rate_constant(r, tmc_rate_arrhenius((tmc_arrhenius_parameters){{ {c_double(p['A'])}, {c_double(p['B'])}, {c_double(p['C'])}, {c_double(p['D'])}, {c_double(p['E'])} }}));")
        elif rx.type in ("TROE", "TERNARY_CHEMICAL_ACTIVATION"):
            fn = "tmc_rate_troe" if rx.type == "TROE" else "tmc_rate_ternary_chemical_activation"
            vals = ", ".join(c_double(p[k]) for k in ("k0_A", "k0_B", "k0_C", "kinf_A", "kinf_B", "kinf_C", "Fc", "N"))
            out(f"  tmc_reaction_set_rate_constant(r, {fn}((tmc_falloff_parameters){{ {vals} }}));")
        elif rx.type == "TERNARY_JPL19":
            vals = ", ".join(c_double(p[k]) for k in ("k0_A", "k0_B", "k0_C", "k0_D", "kinf_A", "kinf_B", "kinf_C", "kinf_D", "kint_A", "kint_B", "kint_C", "kint_D", "Fc", "N"))
            out(f"  tmc_reaction_set_rate_constant(r, tmc_rate_ternary_jpl19((tmc_ternary_jpl19_parameters){{ {vals} }}));")
        elif rx.type == "BRANCHED":
            out(f"  tmc_reaction_set_rate_constant(r, tmc_rate_branched((tmc_branched_parameters){{ TMC_BRANCH_{rx.branch}, {c_double(p['X'])}, {c_double(p['Y'])}, {c_double(p['a0'])}, {int(p['n'])} }}));")
        elif rx.type == "TUNNELING":
            out(f"  tmc_reaction_set_rate_constant(r, tmc_rate_tunneling((tmc_tunneling_parameters){{ {c_double(p['A'])}, {c_double(p['B'])}, {c_double(p['C'])} }}));")
        elif rx.type == "TAYLOR_SERIES":
            coefs = p["coefficients"]
            out("  { tmc_taylor_series_parameters t = tmc_taylor_series_defaults();")
            out(f"    t.A = {c_double(p['A'])}; t.B = {c_double(p['B'])}; t.C = {c_double(p['C'])}; t.D = {c_double(p['D'])}; t.E = {c_double(p['E'])};")
            out(f"    t.n_coefficients = {len(coefs)};")
            for i, v in enumerate(coefs):
                out(f"    t.coefficients[{i}] = {c_double(v)};")
            out("    tmc_reaction_set_rate_constant(r, tmc_rate_taylor_series(t)); }")
        elif rx.type == "REVERSIBLE":
            out(f"  tmc_reaction_set_rate_constant(r, tmc_rate_reversible((tmc_reversible_parameters){{ {c_double(p['A'])}, {c_double(p['C'])}, {c_double(p['k_r'])} }}));")
        elif rx.type == "USER_DEFINED":
            out(f"  tmc_reaction_set_rate_constant(r, tmc_rate_user_defined({c_string(rx.label)}, {c_double(rx.scaling_factor)}));")
        elif rx.type == "SURFACE":
            out(f"  tmc_reaction_set_rate_constant(r, tmc_rate_surface({c_string(rx.label)}, tmc_phase_find(phase, {c_string(rx.surface_species)}), {c_double(rx.reaction_probability)}));")
        else:
            raise ValueError(f"cannot emit rate type {rx.type}")
        out("  tmc_reaction_set_phase(r, phase);")
        out("  ++r;")
    out("  tmc_system_free(&system);")
    out("  *reactions = list;")
    out(f"  return {len(m.reactions)};")
    out("}")
    out()
    out(f"int {name}_initial_state(tmc_state *state)")
    out("{")
    out("  int count = 0;")
    for s in m.species:
        if s.initial_concentration is not None:
            out(f"  count += tmc_state_set_concentration(state, {c_string(s.name)}, {c_double(s.initial_concentration)}, NULL) == TMC_OK;")
    for label, value in m.rate_parameters.items():
        out(f"  count += tmc_state_set_custom_rate_parameter(state, {c_string(label)}, {c_double(value)}, NULL) == TMC_OK;")
    if m.temperature is not None and m.pressure is not None:
        out("  { tmc_conditions *cond = tmc_state_conditions(state);")
        out("    for (int i = 0; i < tmc_state_num_cells(state); ++i) {")
        out(f"      cond[i].temperature = {c_double(m.temperature)}; cond[i].pressure = {c_double(m.pressure)};")
        out("      tmc_conditions_calculate_ideal_air_density(&cond[i]); } }")
    out("  return count;")
    out("}")
    out()


# ---------------------------------------------------------------------------
# Kernels
# ---------------------------------------------------------------------------

def emit_kernels(out, name, c: mech.Compiled, width, chunk, reciprocal, table_forcing, table_jacobian):
    n = len(c.variables)
    nnz = len(c.lu_elements)
    nr = len(c.reactions)
    W = "TMC_WIDTH"
    lu = c.lu_index

    def Y(v):
        return f"y[{v}*{W}+c]"

    def F(v):
        return f"f[{v}*{W}+c]"

    def RC(r):
        return f"rc[{r}*{W}+c]"

    def A(e):
        return f"a[{lu[e]}*{W}+c]"

    def X(v):
        return f"x[{v}*{W}+c]"

    lane = f"for (int c = 0; c < {W}; ++c)"
    out(f"#if TMC_WIDTH != {width}")
    out(f'#error "{name} kernels were generated for TMC_WIDTH == {width}"')
    out("#endif")
    out(f"#define {name.upper()}_NV {n}")
    out(f"#define {name.upper()}_NR {nr}")
    out(f"#define {name.upper()}_NNZ {nnz}")
    out()

    # ---- forcing: one straight-line pass per block over all reactions
    touched = set()
    out("/* forcing f = F(y): rate = k * prod(reactants); reactants lose it, products gain yield * rate */")
    funcs = []
    for start in range(0, nr, chunk * 4):
        fn = f"{name}_forcing_{start}"
        funcs.append(fn)
        out(f"static void {fn}(const double *restrict y, const double *restrict rc, double *restrict f, int nb)")
        out("{")
        out("  for (int b = 0; b < nb; ++b, y += " + f"{n}*{W}, rc += {nr}*{W}, f += {n}*{W})")
        out("  {")
        out(f"    double rate[{W}];")
        if start == 0:
            for v in range(n):
                if not any(v in rs or any(pv == v for pv, _ in ps) for rs, ps in c.processes):
                    out(f"    {lane} {F(v)} = 0.0;")
        for r in range(start, min(nr, start + chunk * 4)):
            rs, ps = c.processes[r]
            out(f"    {lane} rate[c] = {RC(r)};")
            for v in rs:
                out(f"    {lane} rate[c] *= {Y(v)};")
            for v in rs:
                op = "-=" if v in touched else "= -"
                touched.add(v)
                out(f"    {lane} {F(v)} {op} rate[c];")
            for v, yld in ps:
                if v in touched:
                    out(f"    {lane} {F(v)} += {c_double(yld)} * rate[c];")
                else:
                    out(f"    {lane} {F(v)} = {c_double(yld)} * rate[c];")
                touched.add(v)
        out("  }")
        out("}")
        out()
    out(f"static void {name}_forcing(const double *restrict y, const double *restrict rc, double *restrict f, int nb)")
    out("{")
    for fn in funcs:
        out(f"  {fn}(y, rc, f, nb);")
    out("}")
    out()

    # ---- Jacobian: alpha I - J in the LU pattern, MICM's job order, first touch assigns
    out("/* alpha I - dF/dy into the LU matrix: for each (independent variable, reaction) job,")
    out("   d = k * prod(other reactants); reactant rows gain d, product rows lose yield * d */")
    touched = set()
    diag = {i: lu[(i, i)] for i in range(n)}
    funcs = []
    njobs = len(c.jobs)
    for start in range(0, njobs, chunk * 4):
        fn = f"{name}_jacobian_{start}"
        funcs.append(fn)
        out(f"static void {fn}(const double *restrict y, const double *restrict rc, double *restrict a, int nb)")
        out("{")
        out("  for (int b = 0; b < nb; ++b, y += " + f"{n}*{W}, rc += {nr}*{W}, a += {nnz}*{W})")
        out("  {")
        out(f"    double d[{W}];")
        if start == 0:
            contributed = set()
            for proc, ind, dep, ps in c.jobs:
                for v in dep:
                    contributed.add(lu[(v, ind)])
                contributed.add(lu[(ind, ind)])
                for v, _ in ps:
                    contributed.add(lu[(v, ind)])
            for e in range(nnz):
                if e not in contributed:
                    out(f"    {lane} a[{e}*{W}+c] = 0.0;")
        for q in range(start, min(njobs, start + chunk * 4)):
            proc, ind, dep, ps = c.jobs[q]
            out(f"    {lane} d[c] = {RC(proc)};")
            for v in dep:
                out(f"    {lane} d[c] *= {Y(v)};")
            for v in dep + [ind]:
                e = lu[(v, ind)]
                op = "+=" if e in touched else "="
                touched.add(e)
                out(f"    {lane} a[{e}*{W}+c] {op} d[c];")
            for v, yld in ps:
                e = lu[(v, ind)]
                if e in touched:
                    out(f"    {lane} a[{e}*{W}+c] -= {c_double(yld)} * d[c];")
                else:
                    out(f"    {lane} a[{e}*{W}+c] = -({c_double(yld)} * d[c]);")
                touched.add(e)
        out("  }")
        out("}")
        out()
    out(f"static void {name}_jacobian(const double *restrict y, const double *restrict rc, double alpha, double *restrict a, int nb)")
    out("{")
    for fn in funcs:
        out(f"  {fn}(y, rc, a, nb);")
    out(f"  for (int b = 0; b < nb; ++b, a += {nnz}*{W})")
    out("  {")
    for i in range(n):
        out(f"    {lane} a[{diag[i]}*{W}+c] += alpha;")
    out("  }")
    out("}")
    out()

    # ---- LU: MICM's Mozart in-place elimination, pivot by pivot, blocks in lock step
    rows = {i: sorted(j for (r, j) in c.lu_elements if r == i) for i in range(n)}
    cols = {j: sorted(i for (i, cc) in c.lu_elements if cc == j) for j in range(n)}
    out("/* in-place LU (MICM's Mozart order): for each pivot, scale the column below by the")
    out("   pivot's reciprocal, then update the trailing submatrix; blocks run in lock step */")
    funcs = []
    for start in range(0, n, chunk):
        fn = f"{name}_lu_{start}"
        funcs.append(fn)
        out(f"static void {fn}(double *restrict a0, double *restrict r0, int nb)")
        out("{")
        for i in range(start, min(n, start + chunk)):
            below = [j for j in cols[i] if j > i]
            right = [k for k in rows[i] if k > i]
            out(f"  for (int b = 0; b < nb; ++b) {{ /* pivot {i} */")
            out(f"    double *restrict a = a0 + (size_t)b * ({nnz}*{W}); double *restrict inv = r0 + (size_t)b * ({n}*{W}) + {i}*{W};")
            out(f"    {lane} inv[c] = 1.0 / {A((i, i))};")
            for j in below:
                out(f"    {lane} {A((j, i))} *= inv[c];")
            for j in below:
                for k in right:
                    out(f"    {lane} {A((j, k))} -= {A((j, i))} * {A((i, k))};")
            out("  }")
        out("}")
        out()
    out(f"static void {name}_lu(double *restrict a, double *restrict recip, int nb)")
    out("{")
    for fn in funcs:
        out(f"  {fn}(a, recip, nb);")
    out("}")
    out()

    # ---- solves: forward (unit lower) then backward, blocks in lock step
    out("/* L U x = b: forward substitution with unit diagonal, then backward substitution */")
    funcs = []
    for start in range(0, n, chunk * 2):
        fn = f"{name}_fwd_{start}"
        funcs.append(fn)
        out(f"static void {fn}(double *restrict x0, const double *restrict a0, int nb)")
        out("{")
        for i in range(start, min(n, start + chunk * 2)):
            js = [j for j in rows[i] if j < i]
            if not js:
                continue
            out(f"  for (int b = 0; b < nb; ++b) {{ double t[{W}]; double *restrict x = x0 + (size_t)b * ({n}*{W}); const double *restrict a = a0 + (size_t)b * ({nnz}*{W});")
            out(f"    {lane} t[c] = {X(i)};")
            for j in js:
                out(f"    {lane} t[c] -= {A((i, j))} * {X(j)};")
            out(f"    {lane} {X(i)} = t[c];")
            out("  }")
        out("}")
        out()
    for start in range(0, n, chunk * 2):
        fn = f"{name}_bwd_{start}"
        funcs.append(fn)
        out(f"static void {fn}(double *restrict x0, const double *restrict a0, const double *restrict r0, int nb)")
        out("{")
        if not reciprocal:
            out("  (void)r0;")
        for r in range(start, min(n, start + chunk * 2)):
            i = n - 1 - r
            js = [j for j in rows[i] if j > i]
            recip_decl = f" const double *restrict recip = r0 + (size_t)b * ({n}*{W});" if reciprocal else ""
            out(f"  for (int b = 0; b < nb; ++b) {{ double t[{W}]; double *restrict x = x0 + (size_t)b * ({n}*{W}); const double *restrict a = a0 + (size_t)b * ({nnz}*{W});{recip_decl}")
            out(f"    {lane} t[c] = {X(i)};")
            for j in js:
                out(f"    {lane} t[c] -= {A((i, j))} * {X(j)};")
            if reciprocal:
                out(f"    {lane} {X(i)} = t[c] * recip[{i}*{W}+c];")
            else:
                out(f"    {lane} {X(i)} = t[c] / {A((i, i))};")
            out("  }")
        out("}")
        out()
    out(f"static void {name}_solve(double *restrict x, const double *restrict a, const double *restrict recip, int nb)")
    out("{")
    for fn in funcs:
        if "_fwd_" in fn:
            out(f"  {fn}(x, a, nb);")
        else:
            out(f"  {fn}(x, a, recip, nb);")
    out("  (void)recip;")
    out("}")
    out()
    out(f"static const tmc_kernels {name}_kernel_table = {{ {width}, {c.fingerprint}ULL, {name}_forcing, {name}_jacobian, {name}_lu, {name}_solve }};")
    out(f"const tmc_kernels *{name}_kernels(void)")
    out("{")
    out(f"  return &{name}_kernel_table;")
    out("}")
    out()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config")
    ap.add_argument("output_stem")
    ap.add_argument("--name", required=True, help="C identifier prefix")
    ap.add_argument("--width", type=int, default=2, help="vector width the kernels are compiled for")
    ap.add_argument("--chunk", type=int, default=16, help="pivots per generated LU function")
    ap.add_argument("--no-kernels", action="store_true", help="emit only the mechanism tables")
    ap.add_argument("--reciprocal", action="store_true", help="back substitution multiplies by stored pivot reciprocals")
    ap.add_argument("--table-forcing", action="store_true")
    ap.add_argument("--table-jacobian", action="store_true")
    args = ap.parse_args()

    m = mech.read(args.config)
    c = mech.compile_mechanism(m)
    name = args.name
    out = Emitter()
    out(f"/* Generated by tools/generate.py from {Path(args.config).name}: {len(m.species)} species, {len(m.reactions)} reactions. Do not edit. */")
    out(f'#include "{Path(args.output_stem).name}.h"')
    out("#include <stddef.h>")
    out("#include <stdlib.h>")
    out()
    emit_mechanism(out, name, m, c)
    if args.no_kernels:
        out(f"const tmc_kernels *{name}_kernels(void)\n{{\n  return NULL;\n}}")
    else:
        emit_kernels(out, name, c, args.width, args.chunk, args.reciprocal, args.table_forcing, args.table_jacobian)
    Path(args.output_stem + ".c").write_text(out.text())
    guard = f"{name.upper()}_GENERATED_H"
    header = Emitter()
    header(f"/* Generated by tools/generate.py. Mechanism {m.name}: {len(m.species)} species, {len(m.reactions)} reactions. */")
    header(f"#ifndef {guard}")
    header(f"#define {guard}")
    header('#include "tmc.h"')
    header(f"#define {name.upper()}_NUM_SPECIES {len(m.species)}")
    header(f"#define {name.upper()}_NUM_VARIABLES {len(c.variables)}")
    header(f"#define {name.upper()}_NUM_REACTIONS {len(m.reactions)}")
    header(f"#define {name.upper()}_FINGERPRINT {c.fingerprint}ULL")
    header(f"tmc_system {name}_system(void);")
    header(f"int {name}_reactions(tmc_reaction **reactions);")
    header(f"int {name}_initial_state(tmc_state *state);")
    header(f"const tmc_kernels *{name}_kernels(void);")
    header("#endif")
    Path(args.output_stem + ".h").write_text(header.text())
    print(f"{args.output_stem}.c: {len(c.variables)} variables, {len(c.reactions)} reactions, LU nnz {len(c.lu_elements)}, "
          f"{len(c.jobs)} Jacobian jobs, fingerprint {c.fingerprint:016x}")


if __name__ == "__main__":
    main()
