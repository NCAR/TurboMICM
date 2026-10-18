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
PANEL = 64  # pivots per LU panel
ROW_BUDGET = 64  # solve terms per block loop body
ORDERING_ENUM = {"markowitz": "TMC_ORDERING_MARKOWITZ", "minimum-fill": "TMC_ORDERING_MINIMUM_FILL"}
LANE_PRAGMA = True  # each lane loop an OpenMP simd loop (Derecho default; --no-lane-pragma for plain loops)
SOLVE_PAIRS = False  # the solves take two blocks per loop iteration, statements interleaved (--solve-pairs)
ASSEMBLE_IN_LU = True  # alpha I - J assembled inside the LU at each element's first touch (--no-assemble-in-lu)
NBMAX = 8  # blocks whose partials the assembling LU holds at once
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
            out(f"  count += tmc_state_fill_concentration(state, {c_string(s.name)}, {c_double(s.initial_concentration)}, NULL) == TMC_OK;")
    for label, value in m.rate_parameters.items():
        out(f"  count += tmc_state_fill_custom_rate_parameter(state, {c_string(label)}, {c_double(value)}, NULL) == TMC_OK;")
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

def emit_kernels(final_out, name, c: mech.Compiled, width, chunk, fast, ordering):
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

    # each lane loop is one vector of the width: by default an OpenMP simd loop,
    # which GCC's loop vectorizer takes whole instead of its SLP vectorizer
    # rebuilding vectors from the unrolled scalars (Derecho tuning: half the
    # instructions, no shuffles, few spills); --no-lane-pragma for plain loops
    lane = ('_Pragma("omp simd") ' if LANE_PRAGMA else "") + f"for (int c = 0; c < {W}; ++c)"

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
    out(f"static _Thread_local _Alignas(64) double {name}_rates[{nr}*{W}];")
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
    # the forcing at y + a k, the stage combination formed block by block into a
    # thread-local buffer (one fused multiply-add per element, as the library's pass)
    out(f"static _Thread_local _Alignas(64) double {name}_ytmp[{n}*{W}];")
    out(f"static void {name}_forcing_at(const double *restrict y, const double *restrict k, double a, const double *restrict rc, double *restrict f, int nb)")
    out("{")
    out(f"  double *restrict rate = {name}_rates;")
    out(f"  double *restrict yt = {name}_ytmp;")
    out(f"  for (int b = 0; b < nb; ++b, y += {n}*{W}, k += {n}*{W}, rc += {nr}*{W}, f += {n}*{W})")
    out("  {")
    out(f"    for (int i = 0; i < {n}*{W}; ++i)")
    out("      yt[i] = y[i] + a * k[i];")
    for fn in funcs:
        if "_rates_" in fn:
            out(f"    {fn}(yt, rc, rate);")
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
    out(f"static _Thread_local _Alignas(64) double {name}_partials[{njobs}*{W}];")
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
                        # the rounded product: a bare product followed by another term would
                        # leave the compiler two ways to contract the pair (Derecho: GCC's choice
                        # differed between this pass and the assembling LU)
                        out(f"    {lane} t{k}[c] = fma({K(-yld)}, d[{q}*{W}+c], 0.0);")
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

    # ---- LU: MICM's Mozart in-place elimination in panels of PANEL pivots. Every element
    # is finished target by target, loaded and stored once per panel with its terms in
    # pivot order, so it sees exactly the update sequence of the pivot-by-pivot algorithm:
    # inside the panel at its own pivot (the diagonal, the column below it with the
    # scaling, the row to its right), outside it in the panel's rank-PANEL trailing
    # update. Blocks run in lock step.
    rows = {i: sorted(j for (r, j) in c.lu_elements if r == i) for i in range(n)}
    cols = {j: sorted(i for (i, cc) in c.lu_elements if cc == j) for j in range(n)}
    row_set = {i: set(rows[i]) for i in range(n)}
    col_set = {j: set(cols[j]) for j in range(n)}
    out("/* in-place LU (MICM's Mozart order), pivots in panels of " + str(PANEL) + ": at each pivot its")
    out("   diagonal, its column (updates from the panel's earlier pivots, then the scaling by the")
    out("   reciprocal), and its row are finished target by target; the trailing matrix then takes")
    out("   the panel's update; blocks run in lock step */")
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
            # Within the panel every element is finished target by target: at pivot p the
            # diagonal, then the column below it (updates from the panel's earlier pivots in
            # pivot order, then the scaling by the reciprocal) and the row to its right. The
            # sources, column p above the diagonal and row p left of it, were finished at
            # their own pivots, so each element sees exactly the pivot-by-pivot sequence.
            for p in range(i0, i1):
                lp = [q for q in rows[p] if i0 <= q < p]  # L entries (p, q) of row p inside the panel
                up = [q for q in cols[p] if i0 <= q < p]  # U entries (q, p) of column p inside the panel
                below = [j for j in cols[p] if j > p]
                right = [k for k in rows[p] if k > p]
                diag_terms = [q for q in lp if q in col_set[p]]
                out(f"    {{ /* pivot {p} */")
                out(f"      double inv[{W}];")
                if diag_terms:
                    out(f"      {{ double d[{W}];")
                    out(f"        {lane} d[c] = {A((p, p))};")
                    for q in diag_terms:
                        out(f"        {lane} d[c] -= {A((p, q))} * {A((q, p))};")
                    out(f"        {lane} {A((p, p))} = d[c];")
                    out(f"        {lane} inv[c] = 1.0 / d[c]; }}")
                else:
                    out(f"      {lane} inv[c] = 1.0 / {A((p, p))};")
                out(f"      {lane} recip[{p}*{W}+c] = inv[c];")
                if below:
                    col_targets = [(j, [q for q in up if q in row_set[j]]) for j in below]
                    used = [q for q in up if any(q in terms for _, terms in col_targets)]
                    out("      { " + " ".join(f"double u{q - i0}[{W}];" for q in used))
                    for q in used:
                        out(f"        {lane} u{q - i0}[c] = {A((q, p))};")
                    for j, terms in col_targets:
                        out(f"        {{ double t[{W}];")
                        out(f"          {lane} t[c] = {A((j, p))};")
                        for q in terms:
                            out(f"          {lane} t[c] -= {A((j, q))} * u{q - i0}[c];")
                        out(f"          {lane} {A((j, p))} = t[c] * inv[c]; }}")
                    out("      }")
                row_targets = [(k, [q for q in lp if k in row_set[q]]) for k in right]
                row_targets = [(k, terms) for k, terms in row_targets if terms]
                if row_targets:
                    used = [q for q in lp if any(q in terms for _, terms in row_targets)]
                    out("      { " + " ".join(f"double l{q - i0}[{W}];" for q in used))
                    for q in used:
                        out(f"        {lane} l{q - i0}[c] = {A((p, q))};")
                    for k, terms in row_targets:
                        out(f"        {{ double t[{W}];")
                        out(f"          {lane} t[c] = {A((p, k))};")
                        for q in terms:
                            out(f"          {lane} t[c] -= l{q - i0}[c] * {A((q, k))};")
                        out(f"          {lane} {A((p, k))} = t[c]; }}")
                    out("      }")
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

    # ---- LU with the Jacobian assembled at first touch (a Derecho tuning). The same
    # elimination, but an element's starting value is its sum of partials (alpha last on
    # the diagonal) formed where the LU first reads it, so the assembled matrix is never
    # stored and reloaded; the result is bit for bit that of assembly followed by the LU.
    # The partials of up to NBMAX blocks sit in one thread-local buffer, so the blocks of
    # a chunk still run in lock step through each panel.
    lua_funcs = []
    if ASSEMBLE_IN_LU:
        out(f"#define {name.upper()}_NBMAX {NBMAX}")
        out(f"static _Thread_local _Alignas(64) double {name}_partials_all[{NBMAX}*{njobs}*{W}];")
        touched = set()

        def assembled(e, tn):
            """tn[c] = element e's assembled value: its terms in MICM's order, alpha last on the diagonal."""
            lines = []
            terms = eterms[e]
            if not terms:
                lines.append(f"{lane} {tn}[c] = 0.0;")
            for step, (q, yld) in enumerate(terms):
                first = step == 0
                if yld is None:
                    lines.append(f"{lane} {tn}[c] {'=' if first else '+='} d[{q}*{W}+c];")
                elif first:
                    lines.append(f"{lane} {tn}[c] = fma({K(-yld)}, d[{q}*{W}+c], 0.0);")
                else:
                    lines.append(f"{lane} {tn}[c] -= {K(yld)} * d[{q}*{W}+c];")
            if e in diag_elements:
                lines.append(f"{lane} {tn}[c] += alpha;")
            return lines

        def start_value(e, tn, indent):
            """tn = the element's current value: assembled at its first touch, loaded
            afterwards. Returns whether tn is a literal zero (a fill element's first
            touch), whose first update must then be the rounded product fma(-l, u, 0.0),
            as the stored zero gave: 0.0 - l*u would fold into an unrounded product."""
            if e in touched:
                out(f"{indent}{lane} {tn}[c] = a[{e}*{W}+c];")
                return False
            touched.add(e)
            for line in assembled(e, tn):
                out(indent + line)
            return not eterms[e]

        def update(tn, left, right, zero_first, indent):
            """tn -= left * right, or the rounded product itself when tn is a literal zero."""
            if zero_first:
                out(f"{indent}{lane} {tn}[c] = fma(-{left}, {right}, 0.0);")
            else:
                out(f"{indent}{lane} {tn}[c] -= {left} * {right};")

        out("/* the LU with alpha I - J assembled on the way: each element's partial sums are formed")
        out("   where the elimination first reads it, so the matrix is never stored and reloaded */")
        for start in range(0, len(panels), per_fn):
            fn = f"{name}_lua_{panels[start][0]}"
            lua_funcs.append(fn)
            out(f"static void {fn}(double *restrict a0, double *restrict r0, const double *restrict d0, double alpha, int nb)")
            out("{")
            out("  (void)alpha;")
            for i0, i1 in panels[start:start + per_fn]:
                out(f"  for (int b = 0; b < nb; ++b) {{ /* pivots {i0}-{i1 - 1} */")
                out(f"    double *restrict a = a0 + (size_t)b * ({nnz}*{W}); double *restrict recip = r0 + (size_t)b * ({n}*{W}); const double *restrict d = d0 + (size_t)b * ({njobs}*{W}); (void)d;")
                for p in range(i0, i1):
                    lp = [q for q in rows[p] if i0 <= q < p]
                    up = [q for q in cols[p] if i0 <= q < p]
                    below = [j for j in cols[p] if j > p]
                    right = [k for k in rows[p] if k > p]
                    diag_terms = [q for q in lp if q in col_set[p]]
                    out(f"    {{ /* pivot {p} */")
                    out(f"      double inv[{W}];")
                    out(f"      {{ double dg[{W}];")
                    start_value(lu[(p, p)], "dg", "        ")
                    for q in diag_terms:
                        out(f"        {lane} dg[c] -= {A((p, q))} * {A((q, p))};")
                    out(f"        {lane} {A((p, p))} = dg[c];")
                    out(f"        {lane} inv[c] = 1.0 / dg[c]; }}")
                    out(f"      {lane} recip[{p}*{W}+c] = inv[c];")
                    if below:
                        col_targets = [(j, [q for q in up if q in row_set[j]]) for j in below]
                        used = [q for q in up if any(q in terms for _, terms in col_targets)]
                        out("      { " + " ".join(f"double u{q - i0}[{W}];" for q in used))
                        for q in used:
                            out(f"        {lane} u{q - i0}[c] = {A((q, p))};")
                        for j, terms in col_targets:
                            out(f"        {{ double t[{W}];")
                            zero = start_value(lu[(j, p)], "t", "          ")
                            for idx, q in enumerate(terms):
                                update("t", A((j, q)), f"u{q - i0}[c]", zero and idx == 0, "          ")
                            out(f"          {lane} {A((j, p))} = t[c] * inv[c]; }}")
                        out("      }")
                    # every element right of the diagonal is finished here: with its terms, or
                    # stored assembled if nothing has touched it yet
                    row_targets = [(k, [q for q in lp if k in row_set[q]]) for k in right]
                    row_targets = [(k, terms) for k, terms in row_targets if terms or lu[(p, k)] not in touched]
                    if row_targets:
                        used = [q for q in lp if any(q in terms for _, terms in row_targets)]
                        out("      { " + " ".join(f"double l{q - i0}[{W}];" for q in used))
                        for q in used:
                            out(f"        {lane} l{q - i0}[c] = {A((p, q))};")
                        for k, terms in row_targets:
                            out(f"        {{ double t[{W}];")
                            zero = start_value(lu[(p, k)], "t", "          ")
                            for idx, q in enumerate(terms):
                                update("t", f"l{q - i0}[c]", A((q, k)), zero and idx == 0, "          ")
                            out(f"          {lane} {A((p, k))} = t[c]; }}")
                        out("      }")
                    out("    }")
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
                        zero = start_value(lu[(j, k)], "t", "        ")
                        for idx, p in enumerate(terms):
                            update("t", f"l{p - i0}[c]", A((p, k)), zero and idx == 0, "        ")
                        out(f"        {lane} {A((j, k))} = t[c]; }}")
                    out("    }")
                out("  }")
            out("}")
            out()
        assert touched == set(range(nnz)), "every LU element must be assembled somewhere"
        partials_funcs = [fn for fn in jacobian_funcs if "_partials_" in fn]
        gather_funcs = [fn for fn in forcing_funcs if "_gather_" in fn]
        # alpha I - J assembled and factored, up to NBMAX blocks at a time
        out(f"static void {name}_jacobian_lu(const double *restrict y, const double *restrict rc, double alpha, double *restrict a, double *restrict recip, int nb)")
        out("{")
        out(f"  double *restrict d_all = {name}_partials_all;")
        out(f"  for (int b0 = 0; b0 < nb; b0 += {NBMAX})")
        out("  {")
        out(f"    const int nbc = nb - b0 < {NBMAX} ? nb - b0 : {NBMAX};")
        out("    for (int b = 0; b < nbc; ++b)")
        out("    {")
        out(f"      const double *yb = y + (size_t)(b0 + b) * ({n}*{W}); const double *rcb = rc + (size_t)(b0 + b) * ({nr}*{W}); double *db = d_all + (size_t)b * ({njobs}*{W});")
        for fn in partials_funcs:
            out(f"      {fn}(yb, rcb, db);")
        out("    }")
        for fn in lua_funcs:
            out(f"    {fn}(a + (size_t)b0 * ({nnz}*{W}), recip + (size_t)b0 * ({n}*{W}), d_all, alpha, nbc);")
        out("  }")
        out("}")
        out()
        # the same with the forcing of y, for the start of a step
        out(f"static void {name}_forcing_jacobian_lu(const double *restrict y, const double *restrict rc, double alpha, double *restrict f, double *restrict a, double *restrict recip, int nb)")
        out("{")
        out(f"  double *restrict d_all = {name}_partials_all;")
        out(f"  double *restrict rate = {name}_rates;")
        out(f"  for (int b0 = 0; b0 < nb; b0 += {NBMAX})")
        out("  {")
        out(f"    const int nbc = nb - b0 < {NBMAX} ? nb - b0 : {NBMAX};")
        out("    for (int b = 0; b < nbc; ++b)")
        out("    {")
        out(f"      const double *yb = y + (size_t)(b0 + b) * ({n}*{W}); const double *rcb = rc + (size_t)(b0 + b) * ({nr}*{W}); double *db = d_all + (size_t)b * ({njobs}*{W}); double *fb = f + (size_t)(b0 + b) * ({n}*{W});")
        for fn in partials_funcs:
            out(f"      {fn}(yb, rcb, db);")
        for fn in rate_funcs:
            out(f"      {fn}(yb, rcb, db, rate);")
        for fn in gather_funcs:
            out(f"      {fn}(rate, fb);")
        out("    }")
        for fn in lua_funcs:
            out(f"    {fn}(a + (size_t)b0 * ({nnz}*{W}), recip + (size_t)b0 * ({n}*{W}), d_all, alpha, nbc);")
        out("  }")
        out("}")
        out()

    # ---- solves: forward (unit lower) then backward, blocks in lock step. In fast
    # mode a long row's sum is split over INTERLEAVE independent partial sums and the
    # back substitution multiplies by the stored pivot reciprocal; exact mode keeps
    # MICM's single sum and division.
    out("/* L U x = b: forward substitution with unit diagonal, then backward substitution */")
    funcs = []

    def row_sum_lines(i, js, final, X_, A_, tn, first=None):
        """The statements for t = x_i - sum_j a_ij x_j and `final` (a format with {t}),
        on the block whose vectors X_ and A_ address, with temporaries named tn;
        `first` addresses the initial value x_i (X_ unless the right-hand side is
        a separate array)."""
        first = first or X_
        lines = []
        if fast and len(js) >= SPLIT_MIN:
            P = INTERLEAVE
            lines.append("    " + " ".join(f"double {tn}{k}[{W}];" for k in range(P)))
            lines.append(f"    {lane} {tn}0[c] = {first(i)};")
            for k in range(1, P):
                lines.append(f"    {lane} {tn}{k}[c] = 0.0;")
            for m, j in enumerate(js):
                k = m % P
                if k == 0:
                    lines.append(f"    {lane} {tn}0[c] -= {A_((i, j))} * {X_(j)};")
                else:
                    lines.append(f"    {lane} {tn}{k}[c] += {A_((i, j))} * {X_(j)};")
            for k in range(1, P):
                lines.append(f"    {lane} {tn}0[c] -= {tn}{k}[c];")
            lines.append("    " + lane + " " + final.format(t=f"{tn}0[c]"))
        else:
            lines.append(f"    double {tn}[{W}];")
            lines.append(f"    {lane} {tn}[c] = {first(i)};")
            for j in js:
                lines.append(f"    {lane} {tn}[c] -= {A_((i, j))} * {X_(j)};")
            lines.append("    " + lane + " " + final.format(t=f"{tn}[c]"))
        return lines

    def X2(v):
        return f"x2[{v}*{W}+c]"

    def RHS(v):
        return f"rhs[{v}*{W}+c]"

    def RHS2(v):
        return f"rhs2[{v}*{W}+c]"

    def A2(e):
        return f"a2[{lu[e]}*{W}+c]"

    def emit_row(i, js, final1, final2, first1=None, first2=None):
        """One row on one block, or (--solve-pairs, final2 given) on two blocks with
        the statements interleaved, so their chains run side by side."""
        out("    {")
        if final2 is not None:
            for s1, s2 in zip(row_sum_lines(i, js, final1, X, A, "t", first1), row_sum_lines(i, js, final2, X2, A2, "u", first2)):
                out(s1)
                out(s2)
        else:
            for s1 in row_sum_lines(i, js, final1, X, A, "t", first1):
                out(s1)
        out("    }")

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

    # the block loops of a row group: with --solve-pairs, two blocks per iteration
    # and then the odd block alone; otherwise every block alone
    def block_loops(with_recip, with_rhs=False):
        x_decl = f"double *restrict x = {block_ptr('x0', f'{n}*{W}', 'xb')}; const double *restrict a = {block_ptr('a0', f'{nnz}*{W}', 'ab')};"
        r_decl = f" const double *restrict recip = {block_ptr('r0', f'{n}*{W}', 'rb')};" if with_recip else ""
        b_decl = f" const double *restrict rhs = {block_ptr('b0', f'{n}*{W}', 'bb')};" if with_rhs else ""
        if SOLVE_PAIRS:
            x2_decl = f" double *restrict x2 = x + {n}*{W}; const double *restrict a2 = a + {nnz}*{W};"
            r2_decl = f" const double *restrict recip2 = recip + {n}*{W};" if with_recip else ""
            b2_decl = f" const double *restrict rhs2 = rhs + {n}*{W};" if with_rhs else ""
            yield f"  for (int b = 0; b + 1 < nb; b += 2) {{ {x_decl}{r_decl}{b_decl}{x2_decl}{r2_decl}{b2_decl}", True
            yield f"  if (nb & 1) {{ const int b = nb - 1; {x_decl}{r_decl}{b_decl}", False
        else:
            yield f"  for (int b = 0; b < nb; ++b) {{ {x_decl}{r_decl}{b_decl}", False

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
            for header, pairs in block_loops(False):
                out(header)
                for i in grp:
                    js = [j for j in rows[i] if j < i]
                    if not js:
                        continue
                    emit_row(i, js, f"{X(i)} = {{t}};", f"{X2(i)} = {{t}};" if pairs else None)
                out("  }")
        out("}")
        out()
    # The same forward substitution reading its right-hand side b from a separate
    # array (a Derecho tuning): the stage's solve runs straight from the forcing
    # without the library copying it into K first. Rows without lower terms copy
    # their value. The in-place functions above stay for solves in place: a call
    # with b == x would break the restrict promises.
    fwdb_funcs = []
    for start in range(0, len(fw_groups), per_fn):
        fn = f"{name}_fwdb_{start}"
        fwdb_funcs.append(fn)
        out(f"static void {fn}(double *restrict x0, const double *restrict b0, const double *restrict a0, int nb)")
        out("{")
        out("  (void)a0;")
        for grp in fw_groups[start:start + per_fn]:
            for header, pairs in block_loops(False, True):
                out(header)
                if not any(j < i for i in grp for j in rows[i]):
                    out("    (void)a;" + (" (void)a2;" if pairs else ""))
                for i in grp:
                    js = [j for j in rows[i] if j < i]
                    if not js:
                        out(f"    {lane} {X(i)} = {RHS(i)};")
                        if pairs:
                            out(f"    {lane} {X2(i)} = {RHS2(i)};")
                        continue
                    emit_row(i, js, f"{X(i)} = {{t}};", f"{X2(i)} = {{t}};" if pairs else None, RHS, RHS2)
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
            for header, pairs in block_loops(fast):
                out(header)
                if not ((not fast) or any(j > i for i in grp for j in rows[i])):
                    out("    (void)a;" + (" (void)a2;" if pairs else ""))
                for i in grp:
                    js = [j for j in rows[i] if j > i]
                    if fast:
                        emit_row(i, js, f"{X(i)} = {{t}} * recip[{i}*{W}+c];", f"{X2(i)} = {{t}} * recip2[{i}*{W}+c];" if pairs else None)
                    else:
                        emit_row(i, js, f"{X(i)} = {{t}} / {A((i, i))};", f"{X2(i)} = {{t}} / {A2((i, i))};" if pairs else None)
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
    out(f"static void {name}_solve_from(double *restrict x, const double *restrict b, const double *restrict a, const double *restrict recip, int nb)")
    out("{")
    for fn in fwdb_funcs:
        out(f"  {fn}(x, b, a, nb);")
    for fn in funcs:
        if "_bwd_" in fn:
            out(f"  {fn}(x, a, recip, nb);")
    out("  (void)recip;")
    out("}")
    out()
    jl = f"{name}_jacobian_lu, {name}_forcing_jacobian_lu" if ASSEMBLE_IN_LU else "NULL, NULL"
    out(f"static const tmc_kernels {name}_kernel_table = {{ {width}, {c.fingerprint}ULL, {ORDERING_ENUM[ordering]}, {name}_forcing, {name}_jacobian, {name}_lu, {name}_solve, {name}_forcing_jacobian, {name}_forcing_at, {name}_solve_from, {jl} }};")
    out(f"const tmc_kernels *{name}_kernels(void)")
    out("{")
    out(f"  return &{name}_kernel_table;")
    out("}")
    out()
    final_out.lines += out.lines


def main():
    global PANEL, ROW_BUDGET, INTERLEAVE, LANE_PRAGMA, SPLIT_MIN, SOLVE_PAIRS, ASSEMBLE_IN_LU
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config")
    ap.add_argument("output_stem")
    ap.add_argument("--name", required=True, help="C identifier prefix")
    ap.add_argument("--width", type=int, default=2, help="vector width the kernels are compiled for")
    ap.add_argument("--chunk", type=int, default=16, help="pivots per generated LU function")
    ap.add_argument("--no-kernels", action="store_true", help="emit only the mechanism tables")
    ap.add_argument("--panel", type=int, default=PANEL, help="pivots per LU panel")
    ap.add_argument("--row-budget", type=int, default=ROW_BUDGET, help="solve terms per block loop body")
    ap.add_argument("--interleave", type=int, default=INTERLEAVE, help="independent accumulations emitted side by side")
    ap.add_argument("--split-min", type=int, default=SPLIT_MIN, help="fast mode splits solve rows with at least this many terms")
    ap.add_argument("--solve-pairs", action="store_true", help="the solves take two blocks per loop iteration, statements interleaved")
    ap.add_argument("--assemble-in-lu", dest="assemble_in_lu", action="store_true", default=True,
                    help="assemble alpha I - J inside the LU at each element's first touch (the default)")
    ap.add_argument("--no-assemble-in-lu", dest="assemble_in_lu", action="store_false", help="assemble, store, then factor")
    ap.add_argument("--ordering", choices=("markowitz", "minimum-fill"), default="markowitz",
                    help="species ordering: MICM's diagonal Markowitz, or greedy minimum fill (fewer LU elements)")
    ap.add_argument("--mode", choices=("fast", "exact"), default="fast",
                    help="exact: MICM's operation order bit for bit; fast: split long sums, pivot reciprocals")
    ap.add_argument("--lane-pragma", dest="lane_pragma", action="store_true", default=True,
                    help="emit each lane loop as an OpenMP simd loop (the default)")
    ap.add_argument("--no-lane-pragma", dest="lane_pragma", action="store_false", help="plain lane loops")
    args = ap.parse_args()

    PANEL = max(1, args.panel)
    ROW_BUDGET = max(1, args.row_budget)
    INTERLEAVE = max(1, args.interleave)
    SPLIT_MIN = max(2, args.split_min)
    SOLVE_PAIRS = args.solve_pairs
    ASSEMBLE_IN_LU = args.assemble_in_lu
    LANE_PRAGMA = args.lane_pragma
    m = mech.read(args.config)
    c = mech.compile_mechanism(m, ordering=args.ordering)
    name = args.name
    out = Emitter()
    out(f"/* Generated by tools/generate.py from {Path(args.config).name}: {len(m.species)} species, {len(m.reactions)} reactions. Do not edit. */")
    out(f'#include "{Path(args.output_stem).name}.h"')
    out("#include <math.h>")
    out("#include <stddef.h>")
    out("#include <stdlib.h>")
    out()
    emit_mechanism(out, name, m, c)
    if args.no_kernels:
        out(f"const tmc_kernels *{name}_kernels(void)\n{{\n  return NULL;\n}}")
    else:
        emit_kernels(out, name, c, args.width, args.chunk, args.mode == "fast", args.ordering)
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
    header(f"#define {name.upper()}_ORDERING {ORDERING_ENUM[args.ordering]}")
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
