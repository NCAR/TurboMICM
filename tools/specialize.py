"""Compile a CAMP configuration into fixed-index, cell-vectorized C11 kernels.

Parsing and units are shared with configure.py. This module owns symbolic
ordering/fill and arithmetic emission; solver control lives in the C template.
No chemistry is inferred from mechanism names or from initial concentrations.
"""
import argparse
from dataclasses import dataclass
import pathlib
import re

import configure
from configure import ConfigurationError, number


PRIORITY = {name: i for i, group in enumerate([
    ['CK_ARRHENIUS'], ['CK_TROE'], ['CK_TERNARY'], ['CK_TERNARY_JPL19'],
    ['CK_BRANCH_ALKOXY', 'CK_BRANCH_NITRATE'], ['CK_TUNNELING'],
    ['CK_TAYLOR'], ['CK_REVERSIBLE'], ['CK_USER'], ['CK_SURFACE']
]) for name in group}


def markowitz(pattern):
    """Reference diagonal Markowitz, with stable original-index tie breaks.

    The reference initially omits artificial diagonal entries. Its unsigned
    isolated-row/column costs sort after ordinary costs. Represent that order
    without depending on the host Python integer width.
    """
    a = [set(row) for row in pattern]
    alive = set(range(len(a)))
    order = []
    while alive:
        def cost(v):
            rows = len(a[v] & alive)
            cols = sum(v in a[j] for j in alive)
            value = (rows - 1) * (cols - 1)
            return (value < 0, value, v)
        pivot = min(alive, key=cost)
        order.append(pivot)
        alive.remove(pivot)
        for row in alive:
            if pivot in a[row]:
                a[row].update(a[pivot] & alive)
    return order


@dataclass
class Plan:
    variables: list
    reactions: list
    rows: list
    indices: dict
    derivatives: list
    pivots: list
    ncustom: int

    @property
    def n(self):
        return len(self.variables)

    @property
    def nnz(self):
        return len(self.indices)


def bind(mechanism):
    variables = [i for i, s in enumerate(mechanism.species) if not s.get('parameterization')]
    if not variables:
        raise ConfigurationError('At least one integrated species is required')
    old = {sp: i for i, sp in enumerate(variables)}
    pattern = [set() for _ in variables]
    for r in mechanism.reactions:
        for col in r['reactants']:
            if col not in old:
                continue
            for row in r['reactants'] + [sp for sp, _ in r['products']]:
                if row in old:
                    pattern[old[row]].add(old[col])
    variables = [variables[i] for i in markowitz(pattern)]
    mapping = {sp: i for i, sp in enumerate(variables)}
    reactions = []
    custom = 0
    for r in sorted(mechanism.reactions, key=lambda r: PRIORITY[r['rate']['type']]):
        reactions.append(dict(
            reactants=[mapping[sp] for sp in r['reactants'] if sp in mapping],
            products=[(mapping[sp], y) for sp, y in r['products'] if sp in mapping],
            air=sum(sp not in mapping for sp in r['reactants']),
            rate=r['rate'], custom=custom))
        custom += {'CK_USER': 1, 'CK_SURFACE': 2}.get(r['rate']['type'], 0)
    n = len(variables)
    pattern = [{i} for i in range(n)]
    for r in reactions:
        for col in r['reactants']:
            for row in r['reactants'] + [sp for sp, _ in r['products']]:
                pattern[row].add(col)
    for i in range(n):
        for row in range(i + 1, n):
            if i in pattern[row]:
                pattern[row].update(col for col in pattern[i] if col > i)
    rows = [sorted(row) for row in pattern]
    indices = {pair: i for i, pair in enumerate((row, col) for row in range(n) for col in rows[row])}
    derivatives = []
    for col in range(n):
        for ri, r in enumerate(reactions):
            # Removing the first occurrence matches the generic schedule even
            # when a species appears repeatedly among the reactants.
            for _ in range(r['reactants'].count(col)):
                other = r['reactants'].copy()
                other.remove(col)
                derivatives.append((ri, other,
                                    [indices[row, col] for row in other + [col]],
                                    [(indices[row, col], y) for row, y in r['products']]))
    pivots = []
    for i in range(n):
        lower = [row for row in range(i + 1, n) if i in pattern[row]]
        divisions = [indices[row, i] for row in lower]
        updates = [(indices[row, col], indices[row, i], indices[i, col])
                   for col in rows[i] if col > i for row in lower]
        pivots.append((indices[i, i], divisions, updates))
    return Plan(variables, reactions, rows, indices, derivatives, pivots, custom)


def rate_body(rate):
    """Emit exactly the reference formula, substituting fixed parameters.

    Leave libm calls in C so rounding follows the selected C compiler/libm.
    In particular, do not evaluate transcendental expressions with Python.
    """
    q = {key: number(value) for key, value in rate.items()
         if key not in ('type', 'label', 'coefficients')}
    kind = rate['type']
    def arr(prefix='', denominator=None):
        a, b, c, d = [q[prefix + key] for key in ('A', 'B', 'C', 'D')]
        return f'{a} * exp({c} / t) * pow(t / {denominator or d}, {b})'
    base = arr() + f' * (1.0 + {q["E"]} * pressure)'
    if kind in ('CK_ARRHENIUS', 'CK_TAYLOR'):
        if kind == 'CK_ARRHENIUS':
            return [f'double k = {base};']
        lines = ['double poly = 0, power = 1;']
        for coefficient in rate['coefficients']:
            lines += [f'poly += {number(coefficient)} * power;', 'power *= t;']
        return lines + [f'double k = poly * {base};']
    if kind in ('CK_TROE', 'CK_TERNARY', 'CK_TERNARY_JPL19'):
        jpl = kind == 'CK_TERNARY_JPL19'
        lines = [f'double k0 = {arr("k0_", None if jpl else "300.0")};',
                 f'double kinf = {arr("kinf_", None if jpl else "300.0")};',
                 'double ratio = k0 * air / kinf;']
        broad = f'pow({q["Fc"]}, {q["N"]} / ({q["N"]} + pow(log10(ratio), 2.0)))'
        if jpl:
            lines += [f'double kint = {arr("kint_")};',
                      f'double kf = k0 * air / (1.0 + ratio) * {broad};',
                      'double k = kf + kint * (1.0 - kf / kinf);']
        else:
            lines += [f'double k = k0 * {"air" if kind == "CK_TROE" else "1.0"} / (1.0 + ratio) * {broad};']
        return lines
    if kind == 'CK_TUNNELING':
        return [f'double k = {q["A"]} * exp(-({q["B"]}) / t + {q["C"]} / (t * t * t));']
    if kind in ('CK_BRANCH_ALKOXY', 'CK_BRANCH_NITRATE'):
        # Prepared branch constants are calculated by the C model constructor
        # with the reference operations, outside the timed integration.
        frac = 'z / (z + val)' if kind == 'CK_BRANCH_ALKOXY' else 'val / (val + z)'
        return ['double a = branch_k0 * air, b = 0.43 * pow(t / 298.0, -8.0);',
                'double val = a / (1.0 + a / b) * pow(0.41, 1.0 / (1.0 + pow(log10(a / b), 2.0)));',
                f'double pre = {q["X"]} * exp(-({q["Y"]}) / t);',
                f'double k = pre * ({frac});']
    if kind == 'CK_REVERSIBLE':
        return [f'double k = {q["A"]} * exp({q["C"]} / t) * {q["reverse"]};']
    if kind == 'CK_USER':
        return [f'double k = custom0 * {q["scaling"]};']
    if kind == 'CK_SURFACE':
        return [f'double speed_factor = 8.0 * CK_GAS_CONSTANT / (CK_PI * {q["molecular_weight"]});',
                'double speed = sqrt(speed_factor * t), radius = custom0, count = custom1;',
                f'double k = 4.0 * count * CK_PI * radius * radius / (radius / {q["diffusion"]} + 4.0 / (speed * {q["probability"]}));']
    raise ConfigurationError('Unsupported specialized rate ' + kind)


class Emitter:
    def __init__(self, symbol, batch, split, width):
        self.symbol, self.batch, self.split, self.width = symbol, batch, split, width
        self.lines = []

    def add(self, *lines):
        self.lines.extend(lines)

    def function(self, name, arguments, units, lane_batch=1):
        """Split large kernels into bounded helper functions.

        Each unit is a list of scalar statements. lane_batch adjacent units
        share one SIMD cell loop. Execution order inside each cell is retained.
        """
        helpers = []
        for start in range(0, len(units), self.split):
            helper = f'{self.symbol}_{name}_{len(helpers)}'
            helpers.append(helper)
            versions = ['dynamic', 'fixed'] if self.width else ['dynamic']
            for version in versions:
                signature = arguments
                if version == 'fixed':
                    signature = signature.replace('size_t cells', 'size_t unused_cells', 1)
                self.add(f'static void {helper}_{version}({signature})', '{')
                if version == 'fixed':
                    self.add('  (void)unused_cells;', f'  const size_t cells = {self.width};')
                for begin in range(start, min(start + self.split, len(units)), lane_batch):
                    self.add('#pragma omp simd', '  for (size_t c = 0; c < cells; ++c)', '  {')
                    for unit in units[begin:min(begin + lane_batch, start + self.split, len(units))]:
                        self.add('    {', *('      ' + line for line in unit), '    }')
                    self.add('  }')
                self.add('}')
            names = [arg.split()[-1] for arg in arguments.split(', ')]
            call = ', '.join(names)
            self.add(f'static void {helper}({arguments})', '{')
            if self.width:
                self.add(f'  if (cells == {self.width}) {helper}_fixed({call});', '  else')
            self.add(f'    {helper}_dynamic({call});', '}')
        return helpers


def emit(mechanism, symbol, batch=8, split=32, algebra_batch=1, width=32):
    plan = bind(mechanism)
    e = Emitter(symbol, batch, split, width)
    e.add(configure.emit(mechanism, symbol + '_mechanism'), '#include "internal.h"',
          '#include "chem_specialized.h"')
    ix = lambda array, i: f'{array}[{i} * cells + c]'
    arguments = 'size_t cells, const double *restrict y, const double *restrict rates, double *restrict out'
    units = []
    for ri, r in enumerate(plan.reactions):
        unit = [f'double rate = {ix("rates", ri)};']
        unit += [f'rate *= {ix("y", j)};' for j in r['reactants']]
        unit += [f'{ix("out", j)} -= rate;' for j in r['reactants']]
        unit += [f'{ix("out", j)} += {number(y)} * rate;' for j, y in r['products']]
        units.append(unit)
    helpers = e.function('forcing', arguments, units, batch)
    e.add(f'static void {symbol}_forcing(const ck_state *s, const double *restrict y, double *restrict out)',
          '{', f'  memset(out, 0, {plan.n} * s->cells * sizeof(double));')
    e.add(*(f'  {h}(s->cells, y, s->rates, out);' for h in helpers), '}')
    units = []
    for ri, reactants, targets, products in plan.derivatives:
        unit = [f'double rate = {ix("rates", ri)};']
        unit += [f'rate *= {ix("y", j)};' for j in reactants]
        unit += [f'{ix("out", j)} += rate;' for j in targets]
        unit += [f'{ix("out", j)} -= {number(y)} * rate;' for j, y in products]
        units.append(unit)
    helpers = e.function('jacobian', arguments, units, batch)
    e.add(f'static void {symbol}_jacobian(const ck_state *s, const double *restrict y, double *restrict out)',
          '{', f'  memset(out, 0, {plan.nnz} * s->cells * sizeof(double));')
    e.add(*(f'  {h}(s->cells, y, s->rates, out);' for h in helpers), '}')
    units = []
    for diagonal, divisions, updates in plan.pivots:
        units += [[f'{ix("a", j)} *= 1.0 / {ix("a", diagonal)};'] for j in divisions]
        units += [[f'{ix("a", target)} -= {ix("a", left)} * {ix("a", right)};']
                  for target, left, right in updates]
    helpers = e.function('factor', 'size_t cells, double *restrict a', units, algebra_batch)
    e.add(f'static void {symbol}_factor(size_t cells, double *restrict a)', '{')
    if helpers:
        e.add('  static const ck_factor_kernel kernels[] = {' + ','.join(helpers) + '};',
              f'  ck_run_factor_kernels(cells, a, {len(helpers)}, kernels);')
    else:
        e.add('  (void)cells; (void)a;')
    e.add('}')
    units = []
    for forward, order in [(True, range(plan.n)), (False, reversed(range(plan.n)))]:
        for row in order:
            columns = [col for col in plan.rows[row] if (col < row if forward else col > row)]
            units += [[f'{ix("x", row)} -= {ix("a", plan.indices[row, col])} * {ix("x", col)};']
                      for col in columns]
            if not forward:
                units.append([f'{ix("x", row)} /= {ix("a", plan.indices[row, row])};'])
    helpers = e.function('linear', 'size_t cells, const double *restrict a, double *restrict x', units, algebra_batch)
    e.add(f'static void {symbol}_linear(size_t cells, const double *restrict a, double *restrict x)', '{',
          '  static const ck_solve_kernel kernels[] = {' + ','.join(helpers) + '};',
          f'  ck_run_solve_kernels(cells, a, x, {len(helpers)}, kernels);', '}')
    e.add(f'static ck_status {symbol}_update(ck_state *s)', '{', '  if (!s) return CK_INVALID;',
          '  const size_t cells = s->cells;')
    for ri, r in enumerate(plan.reactions):
        e.add('  {')
        kind = r['rate']['type']
        if kind in ('CK_BRANCH_ALKOXY', 'CK_BRANCH_NITRATE'):
            e.add(f'    double branch_k0 = s->model->rates[{ri}].k0, z = s->model->rates[{ri}].z;')
        # Suppress unused environmental inputs by emitting only those referenced.
        body = rate_body(r['rate'])
        if r['air']:
            body += ['k *= air;'] * r['air']
        text = '\n'.join(body)
        e.add('    for (size_t c = 0; c < cells; ++c)', '    {')
        for name, source in [('t', 's->conditions[c].temperature'),
                             ('pressure', 's->conditions[c].pressure'),
                             ('air', 's->conditions[c].air_density'),
                             ('custom0', ix('s->custom', r['custom'])),
                             ('custom1', ix('s->custom', r['custom'] + 1))]:
            if re.search(r'\b' + name + r'\b', text):
                e.add(f'      double {name} = {source};')
        e.add(*('      ' + line for line in body), f'      {ix("s->rates", ri)} = k;', '    }', '  }')
    e.add('  return CK_OK;', '}')
    # Template macros are local to this translation unit.
    e.add(f'#define SP_PREFIX {symbol}', f'#define SP_N {plan.n}', f'#define SP_NNZ {plan.nnz}')
    e.add('static const size_t sp_diagonal[] = {' + ','.join(str(plan.indices[i, i]) for i in range(plan.n)) + '};',
          '#include "specialized_solver.h"',
          f'const ck_specialized {symbol} = {{ &{symbol}_mechanism, {symbol}_update, {symbol}_forcing,',
          f'  {symbol}_jacobian, {symbol}_factor, {symbol}_linear, {symbol}_solve }};')
    return '\n'.join(e.lines) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('configuration', type=pathlib.Path)
    parser.add_argument('output', type=pathlib.Path)
    parser.add_argument('--symbol', default='specialized')
    parser.add_argument('--batch', type=int, default=8, help='chemistry units per cell loop')
    parser.add_argument('--split', type=int, default=32, help='units per helper function')
    parser.add_argument('--algebra-batch', type=int, default=1, help='LU/solve operations per cell loop')
    parser.add_argument('--width', type=int, default=32, help='fixed-stride fast path; 0 disables it')
    args = parser.parse_args()
    if min(args.batch, args.split, args.algebra_batch, args.width) < 1 or args.width < 0:
        parser.error('batch, split and algebra-batch must be positive')
    content = emit(configure.read(args.configuration), args.symbol, args.batch, args.split, args.algebra_batch, args.width)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(content)


if __name__ == '__main__':
    main()
