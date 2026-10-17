#!/usr/bin/env python3
"""Generate the C for a mechanism: its species and reactions as library calls,
its recorded initial state, and specialized kernels (forcing, Jacobian, sparse
LU, triangular solves) with the species and matrix indices as literals.

    generate.py CONFIG OUTPUT_STEM --name NAME [--width W] [--chunk N]
                [--no-kernels] [--mode fast|exact]

writes OUTPUT_STEM.c and OUTPUT_STEM.h. The kernels compute exactly what the
library's generic kernels compute, in the same order (see DESIGN.md), for
states whose vector width is W; the library checks the width and a
fingerprint of the mechanism when the kernels are attached.
"""
import argparse
import sys

INTERLEAVE = 4  # independent accumulations emitted side by side
SPLIT_MIN = 12  # fast mode splits sums with at least this many terms
PANEL = 4  # pivots per LU panel
ROW_BUDGET = 64  # solve terms per block loop body
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
    out("  (void)state;")
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

def emit_kernels(final_out, name, c: mech.Compiled, width, chunk, fast):
    out = Emitter()
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

    def block_ptr(base, stride, table):
        return f"{base} + (size_t)b * ({stride})"

    def K(value):
        """C literal for a yield or its negative."""
        return c_double(value)

    out(f"#if TMC_WIDTH != {width}")
    out(f'#error "{name} kernels were generated for TMC_WIDTH == {width}"')
    out("#endif")
    out(f"#define {name.upper()}_NV {n}")
    out(f"#define {name.upper()}_NR {nr}")
    out(f"#define {name.upper()}_NNZ {nnz}")
    out()

    # ---- forcing, gather form: the reaction rates once per block, then each species sums
    # its terms in MICM's accumulation order (reaction by reaction; reactants, then products)
    out("/* forcing f = F(y). Pass 1: rate_r = k_r * prod(reactants). Pass 2: each species sums")
    out("   -rate over the reactions consuming it and yield * rate over those producing it. */")
    out(f"static _Thread_local double {name}_rates[{nr}*{W}];")
    terms = {v: [] for v in range(n)}  # species -> [(reaction, yield or None)], in MICM's order
    for r, (rs, ps) in enumerate(c.processes):
        for v in rs:
            terms[v].append((r, None))
        for v, yld in ps:
            terms[v].append((r, yld))
    funcs = []
    for start in range(0, nr, chunk * 4):
        fn = f"{name}_rates_{start}"
        funcs.append(fn)
        out(f"static void {fn}(const double *restrict y, const double *restrict rc, double *restrict rate)")
        out("{")
        for r in range(start, min(nr, start + chunk * 4)):
            rs, ps = c.processes[r]
            out(f"  {{ double t[{W}];")
            out(f"    {lane} t[c] = {RC(r)};")
            for v in rs:
                out(f"    {lane} t[c] *= {Y(v)};")
            out(f"    {lane} rate[{r}*{W}+c] = t[c]; }}")
        out("}")
        out()
    # Species are emitted longest-sum first, INTERLEAVE at a time with their terms
    # interleaved, so the accumulation chains of the dense species overlap.
    order = sorted(range(n), key=lambda v: -len(terms[v]))
    groups = [order[i:i + INTERLEAVE] for i in range(0, n, INTERLEAVE)]
    per_fn = max(1, (chunk * 4) // INTERLEAVE)
    for start in range(0, len(groups), per_fn):
        fn = f"{name}_gather_{start}"
        funcs.append(fn)
        out(f"static void {fn}(const double *restrict rate, double *restrict f)")
        out("{")
        for grp in groups[start:start + per_fn]:
            out("  { " + " ".join(f"double t{k}[{W}];" for k in range(len(grp))))
            longest = max(len(terms[v]) for v in grp)
            for step in range(longest):
                for k, v in enumerate(grp):
                    if step >= len(terms[v]):
                        continue
                    r, yld = terms[v][step]
                    first = step == 0
                    if yld is None:
                        out(f"    {lane} t{k}[c] {'= -' if first else '-='} rate[{r}*{W}+c];")
                    else:
                        out(f"    {lane} t{k}[c] {'=' if first else '+='} {K(yld)} * rate[{r}*{W}+c];")
            for k, v in enumerate(grp):
                if not terms[v]:
                    out(f"    {lane} t{k}[c] = 0.0;")
                out(f"    {lane} {F(v)} = t{k}[c];")
            out("  }")
        out("}")
        out()
    forcing_funcs = list(funcs)
    out(f"static void {name}_forcing(const double *restrict y, const double *restrict rc, double *restrict f, int nb)")
    out("{")
    out(f"  double *restrict rate = {name}_rates;")
    out(f"  for (int b = 0; b < nb; ++b, y += {n}*{W}, rc += {nr}*{W}, f += {n}*{W})")
    out("  {")
    for fn in funcs:
        if "_rates_" in fn:
            out(f"    {fn}(y, rc, rate);")
        else:
            out(f"    {fn}(rate, f);")
    out("  }")
    out("}")
    out()

    # ---- Jacobian, gather form: the job partials once per block, then each LU element sums
    # its terms in MICM's order (jobs by independent variable, then process; within a job the
    # dependent reactants, the independent variable's diagonal, then the products)
    out("/* alpha I - dF/dy into the LU matrix. Pass 1: d_q = k * prod(other reactants) for each")
    out("   (independent variable, reaction) job. Pass 2: each element sums +d over the jobs whose")
    out("   reactant rows it holds and -yield * d over product rows; diagonals add alpha last. */")
    njobs = len(c.jobs)
    out(f"static _Thread_local double {name}_partials[{njobs}*{W}];")
    eterms = {e: [] for e in range(nnz)}  # element -> [(job, yield or None)]
    for q, (proc, ind, dep, ps) in enumerate(c.jobs):
        for v in dep:
            eterms[lu[(v, ind)]].append((q, None))
        eterms[lu[(ind, ind)]].append((q, None))
        for v, yld in ps:
            eterms[lu[(v, ind)]].append((q, yld))
    funcs = []
    for start in range(0, njobs, chunk * 4):
        fn = f"{name}_partials_{start}"
        funcs.append(fn)
        out(f"static void {fn}(const double *restrict y, const double *restrict rc, double *restrict d)")
        out("{")
        for q in range(start, min(njobs, start + chunk * 4)):
            proc, ind, dep, ps = c.jobs[q]
            out(f"  {{ double t[{W}];")
            out(f"    {lane} t[c] = {RC(proc)};")
            for v in dep:
                out(f"    {lane} t[c] *= {Y(v)};")
            out(f"    {lane} d[{q}*{W}+c] = t[c]; }}")
        out("}")
        out()
    diag_elements = {lu[(i, i)] for i in range(n)}
    order = sorted(range(nnz), key=lambda e: -len(eterms[e]))
    groups = [order[i:i + INTERLEAVE] for i in range(0, nnz, INTERLEAVE)]
    per_fn = max(1, (chunk * 8) // INTERLEAVE)
    for start in range(0, len(groups), per_fn):
        fn = f"{name}_assemble_{start}"
        funcs.append(fn)
        out(f"static void {fn}(const double *restrict d, double alpha, double *restrict a)")
        out("{")
        out("  (void)alpha; (void)d;")
        for grp in groups[start:start + per_fn]:
            out("  { " + " ".join(f"double t{k}[{W}];" for k in range(len(grp))))
            longest = max(len(eterms[e]) for e in grp)
            for step in range(longest):
                for k, e in enumerate(grp):
                    if step >= len(eterms[e]):
                        continue
                    q, yld = eterms[e][step]
                    first = step == 0
                    if yld is None:
                        out(f"    {lane} t{k}[c] {'=' if first else '+='} d[{q}*{W}+c];")
                    elif first:
                        out(f"    {lane} t{k}[c] = {K(-yld)} * d[{q}*{W}+c];")
                    else:
                        out(f"    {lane} t{k}[c] -= {K(yld)} * d[{q}*{W}+c];")
            for k, e in enumerate(grp):
                if not eterms[e]:
                    out(f"    {lane} t{k}[c] = 0.0;")
                if e in diag_elements:
                    out(f"    {lane} t{k}[c] += alpha;")
                out(f"    {lane} a[{e}*{W}+c] = t{k}[c];")
            out("  }")
        out("}")
        out()
    jacobian_funcs = list(funcs)
    out(f"static void {name}_jacobian(const double *restrict y, const double *restrict rc, double alpha, double *restrict a, int nb)")
    out("{")
    out(f"  double *restrict d = {name}_partials;")
    out(f"  for (int b = 0; b < nb; ++b, y += {n}*{W}, rc += {nr}*{W}, a += {nnz}*{W})")
    out("  {")
    for fn in funcs:
        if "_partials_" in fn:
            out(f"    {fn}(y, rc, d);")
        else:
            out(f"    {fn}(d, alpha, a);")
    out("  }")
    out("}")
    out()

    # ---- forcing and Jacobian of the same y in one pass: the rates come from the
    # partials (the last reactant's partial times that reactant, MICM's product order),
    # so the reactions are read once
    out("/* forcing and alpha I - dF/dy of the same y: the partials once, the rates from the")
    out("   last reactant's partial, then the two gathers */")
    last_job = {}  # reaction -> job index of its last reactant occurrence
    for q, (proc, ind, dep, ps) in enumerate(c.jobs):
        rs, _ = c.processes[proc]
        # jobs of a reaction are in reactant order within MICM's stable sort by variable;
        # the job whose dependents are all but the last occurrence is the last reactant's
        if dep == rs[:-1]:
            last_job[proc] = (q, rs[-1])
    rate_funcs = []
    for start in range(0, nr, chunk * 4):
        fn = f"{name}_rates_from_partials_{start}"
        rate_funcs.append(fn)
        out(f"static void {fn}(const double *restrict y, const double *restrict rc, const double *restrict d, double *restrict rate)")
        out("{")
        out("  (void)y; (void)rc; (void)d;")
        for r in range(start, min(nr, start + chunk * 4)):
            if r in last_job:
                q, v = last_job[r]
                out(f"  {lane} rate[{r}*{W}+c] = d[{q}*{W}+c] * {Y(v)};")
            else:
                out(f"  {lane} rate[{r}*{W}+c] = {RC(r)};")
        out("}")
        out()
    out(f"static void {name}_forcing_jacobian(const double *restrict y, const double *restrict rc, double alpha, double *restrict f, double *restrict a, int nb)")
    out("{")
    out(f"  double *restrict d = {name}_partials;")
    out(f"  double *restrict rate = {name}_rates;")
    out(f"  for (int b = 0; b < nb; ++b, y += {n}*{W}, rc += {nr}*{W}, f += {n}*{W}, a += {nnz}*{W})")
    out("  {")
    for fn in jacobian_funcs:
        if "_partials_" in fn:
            out(f"    {fn}(y, rc, d);")
    for fn in rate_funcs:
        out(f"    {fn}(y, rc, d, rate);")
    for fn in forcing_funcs:
        if "_gather_" in fn:
            out(f"    {fn}(rate, f);")
    for fn in jacobian_funcs:
        if "_assemble_" in fn:
            out(f"    {fn}(d, alpha, a);")
    out("  }")
    out("}")
    out()

    # ---- LU: MICM's Mozart in-place elimination in panels of PANEL pivots. Within a
    # panel the pivots run one by one on the panel's own rows and columns; the rest of
    # the matrix then takes the panel's updates target by target, each target loaded and
    # stored once with its terms in pivot order, so every element sees exactly the
    # update sequence of the pivot-by-pivot algorithm. Blocks run in lock step.
    rows = {i: sorted(j for (r, j) in c.lu_elements if r == i) for i in range(n)}
    cols = {j: sorted(i for (i, cc) in c.lu_elements if cc == j) for j in range(n)}
    out("/* in-place LU (MICM's Mozart order), pivots in panels: scale the column below each")
    out("   pivot by its reciprocal, update the panel's rows and columns, then apply the panel's")
    out("   rank-" + str(PANEL) + " update to the trailing matrix; blocks run in lock step */")
    funcs = []
    panels = [(i0, min(n, i0 + PANEL)) for i0 in range(0, n, PANEL)]
    per_fn = max(1, chunk // PANEL)
    for start in range(0, len(panels), per_fn):
        fn = f"{name}_lu_{panels[start][0]}"
        funcs.append(fn)
        out(f"static void {fn}(double *restrict a0, double *restrict r0, int nb)")
        out("{")
        for i0, i1 in panels[start:start + per_fn]:
            out(f"  for (int b = 0; b < nb; ++b) {{ /* pivots {i0}-{i1 - 1} */")
            out(f"    double *restrict a = a0 + (size_t)b * ({nnz}*{W}); double *restrict recip = r0 + (size_t)b * ({n}*{W});")
            for p in range(i0, i1):
                below = [j for j in cols[p] if j > p]
                right = [k for k in rows[p] if k > p]
                out(f"    {{ double inv[{W}];")
                out(f"      {lane} inv[c] = 1.0 / {A((p, p))};")
                out(f"      {lane} recip[{p}*{W}+c] = inv[c];")
                for j in below:
                    out(f"      {lane} {A((j, p))} *= inv[c];")
                for j in below:
                    for k in right:
                        if j < i1 or k < i1:
                            out(f"      {lane} {A((j, k))} -= {A((j, p))} * {A((p, k))};")
                out("    }")
            # trailing update, row by row
            trailing_rows = sorted({j for p in range(i0, i1) for j in cols[p] if j >= i1})
            for j in trailing_rows:
                lp = [p for p in range(i0, i1) if j in cols[p]]
                targets = []
                for k in rows[j]:
                    if k < i1:
                        continue
                    terms = [p for p in lp if k in rows[p]]
                    if terms:
                        targets.append((k, terms))
                if not targets:
                    continue
                out("    { " + " ".join(f"double l{p - i0}[{W}];" for p in lp))
                for p in lp:
                    out(f"      {lane} l{p - i0}[c] = {A((j, p))};")
                for k, terms in targets:
                    out(f"      {{ double t[{W}];")
                    out(f"        {lane} t[c] = {A((j, k))};")
                    for p in terms:
                        out(f"        {lane} t[c] -= l{p - i0}[c] * {A((p, k))};")
                    out(f"        {lane} {A((j, k))} = t[c]; }}")
                out("    }")
            out("  }")
        out("}")
        out()
    out(f"static void {name}_lu(double *restrict a, double *restrict recip, int nb)")
    out("{")
    for fn in funcs:
        out(f"  {fn}(a, recip, nb);")
    out("}")
    out()

    # ---- solves: forward (unit lower) then backward, blocks in lock step. In fast
    # mode a long row's sum is split over INTERLEAVE independent partial sums and the
    # back substitution multiplies by the stored pivot reciprocal; exact mode keeps
    # MICM's single sum and division.
    out("/* L U x = b: forward substitution with unit diagonal, then backward substitution */")
    funcs = []

    def emit_row_sum(i, js, final):
        """t = x_i - sum_j a_ij x_j, then `final` (a format with {t}) stores the result."""
        if fast and len(js) >= SPLIT_MIN:
            P = INTERLEAVE
            out("    " + " ".join(f"double t{k}[{W}];" for k in range(P)))
            out(f"    {lane} t0[c] = {X(i)};")
            for k in range(1, P):
                out(f"    {lane} t{k}[c] = 0.0;")
            for m, j in enumerate(js):
                k = m % P
                if k == 0:
                    out(f"    {lane} t0[c] -= {A((i, j))} * {X(j)};")
                else:
                    out(f"    {lane} t{k}[c] += {A((i, j))} * {X(j)};")
            for k in range(1, P):
                out(f"    {lane} t0[c] -= t{k}[c];")
            out("    " + lane + " " + final.format(t="t0[c]"))
        else:
            out(f"    double t[{W}];")
            out(f"    {lane} t[c] = {X(i)};")
            for j in js:
                out(f"    {lane} t[c] -= {A((i, j))} * {X(j)};")
            out("    " + lane + " " + final.format(t="t[c]"))

    # Rows share a block loop until the body holds about ROW_BUDGET terms: short rows
    # would otherwise spend half their instructions on loop control, while a long
    # row alone keeps its body short enough for the blocks' chains to overlap.
    def row_groups(order, terms_of):
        groups, cur, cost = [], [], 0
        for i in order:
            w = max(1, len(terms_of(i)))
            if cur and cost + w > ROW_BUDGET:
                groups.append(cur)
                cur, cost = [], 0
            cur.append(i)
            cost += w
        if cur:
            groups.append(cur)
        return groups

    fw_groups = row_groups(range(n), lambda i: [j for j in rows[i] if j < i])
    bw_groups = row_groups(range(n - 1, -1, -1), lambda i: [j for j in rows[i] if j > i])
    per_fn = 8
    for start in range(0, len(fw_groups), per_fn):
        if not any(j < i for grp in fw_groups[start:start + per_fn] for i in grp for j in rows[i]):
            continue
        fn = f"{name}_fwd_{start}"
        funcs.append(fn)
        out(f"static void {fn}(double *restrict x0, const double *restrict a0, int nb)")
        out("{")
        for grp in fw_groups[start:start + per_fn]:
            if not any(j < i for i in grp for j in rows[i]):
                continue
            out(f"  for (int b = 0; b < nb; ++b) {{ double *restrict x = {block_ptr('x0', f'{n}*{W}', 'xb')}; const double *restrict a = {block_ptr('a0', f'{nnz}*{W}', 'ab')};")
            for i in grp:
                js = [j for j in rows[i] if j < i]
                if not js:
                    continue
                out("    {")
                emit_row_sum(i, js, f"{X(i)} = {{t}};")
                out("    }")
            out("  }")
        out("}")
        out()
    for start in range(0, len(bw_groups), per_fn):
        fn = f"{name}_bwd_{start}"
        funcs.append(fn)
        out(f"static void {fn}(double *restrict x0, const double *restrict a0, const double *restrict r0, int nb)")
        out("{")
        out("  (void)a0; (void)r0;")
        for grp in bw_groups[start:start + per_fn]:
            recip_decl = f" const double *restrict recip = {block_ptr('r0', f'{n}*{W}', 'rb')};" if fast else ""
            uses_a = (not fast) or any(j > i for i in grp for j in rows[i])
            a_decl = f" const double *restrict a = {block_ptr('a0', f'{nnz}*{W}', 'ab')};" if uses_a else ""
            out(f"  for (int b = 0; b < nb; ++b) {{ double *restrict x = {block_ptr('x0', f'{n}*{W}', 'xb')};{a_decl}{recip_decl}")
            for i in grp:
                js = [j for j in rows[i] if j > i]
                out("    {")
                if fast:
                    emit_row_sum(i, js, f"{X(i)} = {{t}} * recip[{i}*{W}+c];")
                else:
                    emit_row_sum(i, js, f"{X(i)} = {{t}} / {A((i, i))};")
                out("    }")
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
    out(f"static const tmc_kernels {name}_kernel_table = {{ {width}, {c.fingerprint}ULL, {name}_forcing, {name}_jacobian, {name}_lu, {name}_solve, {name}_forcing_jacobian }};")
    out(f"const tmc_kernels *{name}_kernels(void)")
    out("{")
    out(f"  return &{name}_kernel_table;")
    out("}")
    out()
    final_out.lines += out.lines


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config")
    ap.add_argument("output_stem")
    ap.add_argument("--name", required=True, help="C identifier prefix")
    ap.add_argument("--width", type=int, default=2, help="vector width the kernels are compiled for")
    ap.add_argument("--chunk", type=int, default=16, help="pivots per generated LU function")
    ap.add_argument("--no-kernels", action="store_true", help="emit only the mechanism tables")
    ap.add_argument("--panel", type=int, default=4, help="pivots per LU panel")
    ap.add_argument("--row-budget", type=int, default=64, help="solve terms per block loop body")
    ap.add_argument("--mode", choices=("fast", "exact"), default="fast",
                    help="exact: MICM's operation order bit for bit; fast: split long sums, pivot reciprocals")
    args = ap.parse_args()

    global PANEL, ROW_BUDGET
    PANEL = max(1, args.panel)
    ROW_BUDGET = max(1, args.row_budget)
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
        emit_kernels(out, name, c, args.width, args.chunk, args.mode == "fast")
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
