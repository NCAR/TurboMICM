"""Generate direct C chemistry and sparse solves from a CAMP configuration.

Normalization and stable ordering are shared with configure.py. All generated
indices refer to that same descriptor and filled CSR layout. The emitter never
changes accumulation order within a cell, even when it fuses cell loops.
"""
import argparse
from pathlib import Path

from configure import normalize, reorder, emit_c, number


def symbolic(mech):
    """Mirror the generic library's symbolic fill and column-major Mozart plan."""
    n = len(mech['species'])
    pattern = {(i, i) for i in range(n)}
    for r in mech['reactions']:
        pattern.update((row, col) for col in r['reactants']
                       for row in r['reactants'] + [i for i, _ in r['products']])
    for i in range(n):
        upper = [k for k in range(i + 1, n) if (i, k) in pattern]
        lower = [j for j in range(i + 1, n) if (j, i) in pattern]
        pattern.update((j, k) for k in upper for j in lower)
    return {ij: k for k, ij in enumerate(sorted(pattern))}


def at(array, index):
    return f'{array}[{index} * nc + c]'


class Emitter:
    def __init__(self, batch, width=0, tile=4):
        self.lines = []
        self.batch = batch
        self.width = width
        self.tile = tile

    def add(self, *lines):
        self.lines.extend(lines)

    def state_width(self, fallback):
        if self.width:
            self.add(f'    if (s->cells != {self.width}) {{ {fallback}; return; }}',
                     f'    const size_t nc = {self.width};')
        else:
            self.add('    const size_t nc = s->cells;')

    def matrix_width(self, fallback, returns=False):
        if self.width:
            result = 'return ' if returns else ''
            tail = '' if returns else ' return;'
            self.add(f'    if (cells != {self.width}) {{ {result}{fallback};{tail} }}',
                     f'    const size_t nc = {self.width};')
        else:
            self.add('    const size_t nc = cells;')

    def begin_tiles(self):
        if self.tile:
            self.add(f'    for (size_t base = 0; base < nc; base += {self.tile}) {{')
            end = f'base + {self.tile}' if self.width and self.width % self.tile == 0 else f'(nc - base < {self.tile} ? nc : base + {self.tile})'
            self.add(f'    const size_t end = {end};', '    (void)end;')

    def end_tiles(self):
        if self.tile:
            self.add('    }')

    @property
    def cell_range(self):
        return 'size_t c = base; c < end; ++c' if self.tile else 'size_t c = 0; c < nc; ++c'

    def loops(self, jobs):
        """Fuse a configurable number of scalar jobs in each SIMD cell loop."""
        for start in range(0, len(jobs), self.batch):
            self.add('    #pragma omp simd', f'    for ({self.cell_range}) {{')
            for job in jobs[start:start + self.batch]:
                self.add('        {', *('            ' + line for line in job), '        }')
            self.add('    }')


def chemistry(e, mech, index):
    n = len(mech['species'])
    e.add('static void sp_forcing(const ChemState *s, const double *restrict y, double *restrict f)',
          '{')
    e.state_width('chem_forcing(s, y, f)')
    e.add('    const double *restrict rates = s->rates;', '    (void)y; (void)rates;',
          f'    memset(f, 0, {n} * nc * sizeof(double));')
    e.begin_tiles()
    jobs = []
    for ri, r in enumerate(mech['reactions']):
        job = [f'double v = {at("rates", ri)};']
        job += [f'v *= {at("y", i)};' for i in r['reactants']]
        job += [f'{at("f", i)} -= v;' for i in r['reactants']]
        job += [f'{at("f", i)} += {number(v)} * v;' for i, v in r['products']]
        jobs.append(job)
    e.loops(jobs)
    e.end_tiles()
    e.add('}', 'static void sp_jacobian(const ChemState *s, const double *restrict y, double *restrict a)',
          '{')
    e.state_width('chem_negative_jacobian(s, y, a)')
    e.add('    const double *restrict rates = s->rates;', '    (void)y; (void)rates;',
          f'    memset(a, 0, {len(index)} * nc * sizeof(double));')
    e.begin_tiles()
    jobs = []
    # A repeated independent reactant skips its FIRST occurrence each time, as
    # in mechanism.c. Do not combine repeats or cancel self-products.
    for col in range(n):
        for ri, r in enumerate(mech['reactions']):
            rr = r['reactants']
            if col not in rr:
                continue
            skip = rr.index(col)
            for _ in range(rr.count(col)):
                others = [i for j, i in enumerate(rr) if j != skip]
                job = [f'double v = {at("rates", ri)};']
                job += [f'v *= {at("y", i)};' for i in others]
                job += [f'{at("a", index[i, col])} += v;' for i in others]
                job += [f'{at("a", index[col, col])} += v;']
                job += [f'{at("a", index[i, col])} -= {number(v)} * v;' for i, v in r['products']]
                jobs.append(job)
    e.loops(jobs)
    e.end_tiles()
    e.add('}')


def linear(e, mech, index, reciprocal):
    n = len(mech['species'])
    e.add('static void sp_alpha(const ChemMechanism *m, size_t cells, double *restrict a, double alpha)',
          '{', '    (void)m;')
    e.matrix_width('chem_alpha_minus_jacobian(m, cells, a, alpha)')
    e.begin_tiles()
    e.loops([[f'{at("a", index[i, i])} += alpha;'] for i in range(n)])
    e.end_tiles()
    e.add('}', 'static ChemStatus sp_factor(const ChemMechanism *m, size_t cells, double *restrict a)',
          '{', '    (void)m;')
    e.matrix_width('chem_factor(m, cells, a)', returns=True)
    e.begin_tiles()
    for i in range(n):
        pivot = at('a', index[i, i])
        e.add(f'    for ({e.cell_range}) if ({pivot} == 0 || !isfinite({pivot})) return CHEM_SINGULAR;')
        if reciprocal:
            e.loops([[f'{pivot} = 1.0 / {pivot};']])
        jobs = []
        for j in range(i + 1, n):
            if (j, i) in index:
                jobs.append([f'{at("a", index[j, i])} *= ' + (pivot if reciprocal else f'1.0 / {pivot}') + ';'])
        for k in range(i + 1, n):
            if (i, k) in index:
                for j in range(i + 1, n):
                    if (j, i) in index:
                        jobs.append([f'{at("a", index[j, k])} -= {at("a", index[j, i])} * {at("a", index[i, k])};'])
        e.loops(jobs)
    e.end_tiles()
    e.add('    return CHEM_OK;', '}',
          'static void sp_backsolve(const ChemMechanism *m, size_t cells, const double *restrict a, double *restrict x)',
          '{', '    (void)m;')
    e.matrix_width('chem_linear_solve(m, cells, a, x)')
    e.begin_tiles()
    jobs = []
    for i in range(n):
        terms = [f'v -= {at("a", index[i, j])} * {at("x", j)};' for j in range(i) if (i, j) in index]
        if terms:
            jobs.append([f'double v = {at("x", i)};', *terms, f'{at("x", i)} = v;'])
    for i in reversed(range(n)):
        terms = [f'v -= {at("a", index[i, j])} * {at("x", j)};' for j in range(i + 1, n) if (i, j) in index]
        operator = '*' if reciprocal else '/'
        jobs.append([f'double v = {at("x", i)};', *terms, f'{at("x", i)} = v {operator} {at("a", index[i, i])};'])
    e.loops(jobs)
    e.end_tiles()
    e.add('}')


def rate_formula(r):
    """C expressions retain the library's grouping and fixed-parameter formulas."""
    p = [number(v) for v in r['values']]
    kind = r['type']
    arr = lambda a, b, c, d: f'{a} * exp({c} / t) * pow(t / {d}, {b})'
    lines = []
    if kind in ('ARRHENIUS', 'TAYLOR'):
        expr = arr(*p[:4]) + f' * (1.0 + {p[4]} * pressure)'
        if kind == 'TAYLOR':
            lines += ['double poly = 0, power = 1;']
            for coef in p[5:]:
                lines += [f'poly += {coef} * power;', 'power *= t;']
            expr = 'poly * ' + expr
    elif kind in ('TROE', 'TERNARY'):
        lines += [f'double k0 = {arr(*p[:3], "300.0")};',
                  f'double kinf = {arr(*p[3:6], "300.0")};',
                  'double ratio = k0 * air / kinf;']
        numerator = 'air' if kind == 'TROE' else '1.0'
        expr = f'k0 * {numerator} / (1.0 + ratio) * pow({p[6]}, {p[7]} / ({p[7]} + pow(log10(ratio), 2.0)))'
    elif kind == 'TERNARY_JPL19':
        lines += [f'double k0 = {arr(*p[:4])};', f'double kinf = {arr(*p[4:8])};',
                  f'double kint = {arr(*p[8:12])};', 'double k0m = k0 * air, ratio = k0m / kinf;',
                  f'double kf = k0m / (1.0 + ratio) * pow({p[12]}, {p[13]} / ({p[13]} + pow(log10(ratio), 2.0)));']
        expr = 'kf + kint * (1.0 - kf / kinf)'
    elif kind == 'TUNNELING':
        expr = f'{p[0]} * exp(-{p[1]} / t + {p[2]} / (t * t * t))'
    elif kind in ('BRANCHED_ALKOXY', 'BRANCHED_NITRATE'):
        lines += [f'double k0 = 2e-22 * 6.02214076e23 * 1e-6 * exp({p[3]});',
                  'double ar = k0 * (2.45e19 / 6.02214076e23 * 1e6);',
                  'double br = .43 * pow(293.0 / 298.0, -8.0);',
                  'double ref = ar / (1.0 + ar / br) * pow(.41, 1.0 / (1.0 + pow(log10(ar / br), 2.0)));',
                  f'double z = ref * (1.0 - {p[2]}) / {p[2]};',
                  'double a = k0 * air, b = .43 * pow(t / 298.0, -8.0);',
                  'double av = a / (1.0 + a / b) * pow(.41, 1.0 / (1.0 + pow(log10(a / b), 2.0)));',
                  f'double pre = {p[0]} * exp(-{p[1]} / t);']
        expr = 'pre * (z / (z + av))' if kind == 'BRANCHED_ALKOXY' else 'pre * (av / (av + z))'
    elif kind == 'REVERSIBLE':
        expr = f'{p[0]} * exp({p[1]} / t) * {p[2]}'
    elif kind == 'USER':
        expr = f'{at("custom", r["parameter"])} * {p[0]}'
    elif kind == 'SURFACE':
        lines += [f'double radius = {at("custom", r["parameter"])}, count = {at("custom", r["parameter"] + 1)};',
                  f'double speed = sqrt((8.0 * (1.380649e-23 * 6.02214076e23) / (3.14159265358979323846 * {p[1]})) * t);']
        expr = f'4.0 * count * 3.14159265358979323846 * radius * radius / (radius / {p[0]} + 4.0 / (speed * {p[2]}))'
    else:
        raise ValueError(f'no generated formula for {kind}')
    lines += [f'double v = {expr};']
    lines += ['v *= air;' for _ in range(r['third_body_count'])]
    for _, multiplier in r['multipliers']:
        a, b, c, d = [number(multiplier.get(k, 0)) for k in ('constant', 'temperature', 'pressure', 'air density')]
        lines += [f'v *= {a} + {b} * t + {c} * pressure + {d} * air;']
    return lines


def rates(e, mech, prefix):
    e.add(f'ChemStatus {prefix}_update_rates(ChemState *s)', '{',
          f'    if (!s || !s->mechanism || s->mechanism->data != &{prefix}_data) return CHEM_INVALID;',
          *( [f'    if (s->cells != {e.width}) return chem_update_rates(s);', f'    const size_t nc = {e.width};'] if e.width else ['    const size_t nc = s->cells;']),
          '    const double *restrict custom = s->custom;',
          '    double *restrict rates = s->rates;', '    (void)custom;',
          '    for (size_t c = 0; c < nc; ++c)',
          '        if (!(s->temperature[c] > 0) || !(s->pressure[c] >= 0)',
          '            || !isfinite(s->temperature[c]) || !isfinite(s->pressure[c])) return CHEM_INVALID;')
    e.begin_tiles()
    jobs = []
    for i, r in enumerate(mech['reactions']):
        jobs.append(['double t = s->temperature[c], pressure = s->pressure[c], air = s->air_density[c];',
                     '(void)t; (void)pressure; (void)air;', *rate_formula(r), f'{at("rates", i)} = v;'])
    e.loops(jobs)
    e.end_tiles()
    e.add(f'    for (size_t i = 0; i < {len(mech["reactions"])} * nc; ++i)',
          '        if (!isfinite(rates[i])) return CHEM_NONFINITE;', '    return CHEM_OK;', '}')


def generate(mech, prefix='specialized', batch=1, width=0, reciprocal=True, fused_stages=True, tile=4):
    e = Emitter(batch, width, tile)
    e.add(emit_c(mech, prefix + '_data'), '#include "internal.h"', '#include <math.h>', '#include <string.h>')
    index = symbolic(mech)
    chemistry(e, mech, index)
    linear(e, mech, index, reciprocal)
    rates(e, mech, prefix)
    if fused_stages:
        e.add('#define CHEM_FUSED_STAGES 1')
    e.add(f'#define CHEM_SOLVE {prefix}_solve', '#define CHEM_FORCING sp_forcing',
          '#define CHEM_JACOBIAN sp_jacobian', '#define CHEM_ALPHA sp_alpha',
          '#define CHEM_FACTOR(m,nc,a) sp_factor(m,nc,a)',
          '#define CHEM_BACKSOLVE(m,nc,a,x) sp_backsolve(m,nc,a,x)',
          f'#define CHEM_SPECIES(s) {len(mech["species"])}', '#define CHEM_CONSTRAINED(m) 0',
          f'#define CHEM_VALID_STATE(s) ((s)->mechanism->data == &{prefix}_data)',
          '#include "rosenbrock_impl.h"')
    return '\n'.join(e.lines) + '\n'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('config')
    p.add_argument('output')
    p.add_argument('--prefix', default='specialized')
    p.add_argument('--no-reorder', action='store_true')
    p.add_argument('--batch', type=int, default=1, help='positive number of operations fused per cell loop')
    p.add_argument('--fixed-width', type=int, default=0, help='constant cell width, with generic kernels for other widths; zero means dynamic')
    p.add_argument('--no-reciprocal-pivots', action='store_true', help='retain ordinary U diagonals instead of their reciprocals')
    p.add_argument('--no-fused-stages', action='store_true', help='use the generic stage vector loops')
    p.add_argument('--tile', type=int, default=4, help='cells processed through each whole kernel at a time; zero disables tiling')
    a = p.parse_args()
    if not a.prefix.isidentifier() or not a.prefix.isascii():
        p.error('prefix must be an ASCII C identifier')
    if a.tile < 0:
        p.error('tile must be nonnegative')
    if a.fixed_width < 0:
        p.error('fixed-width must be nonnegative')
    if a.batch < 1:
        p.error('batch must be positive')
    mech = normalize(a.config)
    if not a.no_reorder:
        mech = reorder(mech)
    Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    Path(a.output).write_text(generate(mech, a.prefix, a.batch, a.fixed_width, not a.no_reciprocal_pivots, not a.no_fused_stages, a.tile))


if __name__ == '__main__':
    main()
