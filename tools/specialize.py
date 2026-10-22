"""Generate a complete mechanism-specialized, portable Rosenbrock-3 solver.

The reader owns units and format validation; Plan owns symbolic algebra; this
module lowers ordered arithmetic to independent-cell SIMD loops. No mechanism
name, reaction count, or chemical identity is special-cased.
"""
import argparse
from pathlib import Path
from configure import read, number
from generate import emit as description
from specialize_rates import formula, USER_TYPES
from symbolic import Plan
from split_kernels import split_kernels
from tiled import fixed_kernels, wrappers, transform_solver


class Emitter:
    def __init__(self):
        self.lines = []

    def add(self, *lines):
        self.lines.extend(lines)

    def loop(self, statements):
        self.add('    #pragma omp simd', '    for (size_t c=0; c<cells; ++c) {')
        self.lines.extend('        ' + s for s in statements)
        self.add('    }')

    def batches(self, operations, count):
        for first in range(0, len(operations), count):
            self.loop([s for op in operations[first:first+count] for s in op])

    def kernel(self, name):
        self.add(f'static void {name}(const chem_model *model, size_t cells,',
                 '    const double *restrict rates, const double *restrict y, double *restrict out)',
                 '{', '    (void)model; (void)rates; (void)y; (void)out; (void)cells;')


def chemistry(e, config, plan, batch):
    e.kernel('forcing')
    operations = []
    for r, rx in enumerate(config['reactions']):
        rs = [plan.variable[s] for s in rx['reactants'] if s in plan.variable]
        ps = [(plan.variable[s], v) for s, v in rx['products'] if s in plan.variable]
        op = ['{', f'double r = rates[{r}*cells+c];']
        op += [f'r *= y[{s}*cells+c];' for s in rs]
        op += [f'out[{s}*cells+c] -= r;' for s in rs]
        op += [f'out[{s}*cells+c] += {number(v)}*r;' for s, v in ps]
        operations.append(op + ['}'])
    e.batches(operations, batch)
    e.add('}')
    e.kernel('jacobian')
    operations = []
    for col, r, rs, ps in plan.jobs:
        op = ['{', f'double r = rates[{r}*cells+c];']
        op += [f'r *= y[{s}*cells+c];' for s in rs]
        op += [f'out[{plan.flat[s, col]}*cells+c] += r;' for s in rs+[col]]
        op += [f'out[{plan.flat[s, col]}*cells+c] -= {number(v)}*r;' for s, v in ps]
        operations.append(op + ['}'])
    e.batches(operations, batch)
    e.add('}')


def algebra(e, plan, batch, reciprocal=False, lu="mozart"):
    e.add('static void add_alpha(size_t cells, double *restrict a, double alpha)', '{')
    e.batches([[f'a[{plan.flat[i, i]}*cells+c] += alpha;'] for i in range(len(plan.variables))], batch)
    e.add('}', 'static void factor(size_t cells, double *restrict a, double *restrict scratch)',
          '{', '    (void)scratch; (void)cells; (void)a;')
    if lu == 'doolittle':
        n = len(plan.variables)
        for i in range(n):
            operations = []
            for k in range(i, n):
                if (i, k) not in plan.pattern:
                    continue
                target = plan.flat[i, k]
                op = ['{', f'double v = a[{target}*cells+c];']
                op += [f'v -= a[{plan.flat[i, j]}*cells+c]*a[{plan.flat[j, k]}*cells+c];'
                       for j in range(i) if (i, j) in plan.pattern and (j, k) in plan.pattern]
                op += [f'a[{target}*cells+c] = v;', '}']
                operations.append(op)
            e.batches(operations, batch)
            diagonal = plan.flat[i, i]
            if reciprocal:
                e.loop([f'a[{diagonal}*cells+c] = 1.0/a[{diagonal}*cells+c];'])
            operations = []
            for k in range(i+1, n):
                if (k, i) not in plan.pattern:
                    continue
                target = plan.flat[k, i]
                op = ['{', f'double v = a[{target}*cells+c];']
                op += [f'v -= a[{plan.flat[k, j]}*cells+c]*a[{plan.flat[j, i]}*cells+c];'
                       for j in range(i) if (k, j) in plan.pattern and (j, i) in plan.pattern]
                op += [f'a[{target}*cells+c] = v{"*" if reciprocal else "/"}a[{diagonal}*cells+c];', '}']
                operations.append(op)
            e.batches(operations, batch)
    else:
        for diagonal, multipliers, updates in plan.elimination():
            if not multipliers and not reciprocal:
                continue
            # A pivot's reciprocal is shared by all its column multipliers.
            e.loop([f'double inv = 1.0/a[{diagonal}*cells+c];'] +
                   ([f'a[{diagonal}*cells+c] = inv;'] if reciprocal else []) +
                   [f'a[{entry}*cells+c] *= inv;' for entry in multipliers])
            e.batches([[f'a[{target}*cells+c] -= a[{left}*cells+c]*a[{right}*cells+c];']
                       for target, left, right in updates], batch)
    e.add('}', 'static void linear(size_t cells, const double *restrict a, double *restrict x)', '{')
    n = len(plan.variables)
    for i in range(n):
        statements = [f'v -= a[{plan.flat[i, j]}*cells+c]*x[{j}*cells+c];'
                      for j in range(i) if (i, j) in plan.pattern]
        if statements:
            e.loop([f'double v = x[{i}*cells+c];'] + statements + [f'x[{i}*cells+c] = v;'])
    for i in reversed(range(n)):
        statements = [f'v -= a[{plan.flat[i, j]}*cells+c]*x[{j}*cells+c];'
                      for j in range(i+1, n) if (i, j) in plan.pattern]
        e.loop([f'double v = x[{i}*cells+c];'] + statements +
               [f'x[{i}*cells+c] = v{"*" if reciprocal else "/"}a[{plan.flat[i, i]}*cells+c];'])
    e.add('}')


def rates(e, config):
    e.add('static void calculate_rates(chem_state *s)', '{',
          '    const size_t cells = s->cells;',
          '    const double *restrict custom = s->custom;',
          '    double *restrict rates = s->rates;', '    (void)custom; (void)rates; (void)cells;')
    offset = 0
    for r, rx in enumerate(config['reactions']):
        setup, body = formula(rx, config['species'], offset)
        e.add('    {', *['    ' + line for line in setup])
        e.loop(['double T = s->conditions[c].temperature;',
                'double P = s->conditions[c].pressure;',
                'double M = s->conditions[c].air_density;',
                '(void)T; (void)P; (void)M;'] + body + [f'rates[{r}*cells+c] = k;'])
        e.add('    }')
        offset += 2 if rx['type'] == 'SURFACE' else int(rx['type'] in USER_TYPES)
    e.add('}')


def emit(config, symbol, chemistry_batch=8, algebra_batch=8, function_lines=256, tile_size=4, flat=False, direct_calls=False, division=False, lu="doolittle"):
    plan = Plan.create(config)
    e = Emitter()
    e.add(description(config, symbol+'_mechanism', kernels=False),
          '#include "chem_specialized.h"', '#include <float.h>', '#include <math.h>', '#include <string.h>', '#include <stdlib.h>')
    kernels = Emitter()
    chemistry(kernels, config, plan, chemistry_batch)
    algebra(kernels, plan, algebra_batch)
    e.lines += split_kernels(kernels.lines, function_lines, omit=() if flat else ('add_alpha',))
    if not flat:
        packed = Emitter()
        chemistry(packed, config, plan, chemistry_batch)
        algebra(packed, plan, algebra_batch, reciprocal=not division, lu=lu)
        e.lines += fixed_kernels(packed.lines, tile_size, function_lines, indirect=not direct_calls)
        e.lines += wrappers(len(plan.variables), len(config['reactions']), len(plan.flat), tile_size)
    rates(e, config)
    template = Path(__file__).with_name('ros3.c.in').read_text()
    if not flat:
        template = transform_solver(template, len(plan.variables), len(config['reactions']), len(plan.flat), tile_size)
    for key, value in [('MECHANISM', symbol+'_mechanism'), ('SPECIES', len(plan.variables)), ('NNZ', len(plan.flat))]:
        template = template.replace('@'+key+'@', str(value))
    e.add(template, f'const chem_specialized {symbol} = {{',
          f'    &{symbol}_mechanism, calculate_rates, forcing, jacobian, factor, linear, specialized_solve', '};')
    return '\n'.join(e.lines) + '\n'


def positive(value):
    n = int(value)
    if n < 1:
        raise argparse.ArgumentTypeError('must be positive')
    return n


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('configuration')
    parser.add_argument('output')
    parser.add_argument('--symbol', default='specialized')
    parser.add_argument('--chemistry-batch', type=positive, default=8)
    parser.add_argument('--algebra-batch', type=positive, default=8)
    parser.add_argument('--function-lines', type=positive, default=256)
    parser.add_argument('--tile-size', type=positive, default=4)
    parser.add_argument('--lu', choices=['mozart', 'doolittle'], default='doolittle')
    parser.add_argument('--division', action='store_true', help='Divide in triangular solves instead of storing reciprocal pivots')
    parser.add_argument('--direct-calls', action='store_true', help='Allow inlining between fixed kernel parts')
    parser.add_argument('--flat', action='store_true', help='Disable the private tiled layout')
    args = parser.parse_args()
    if not args.symbol.isascii() or not args.symbol.isidentifier():
        parser.error('symbol must be a C identifier')
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(emit(read(args.configuration), args.symbol, args.chemistry_batch, args.algebra_batch, args.function_lines, args.tile_size, args.flat, args.direct_calls, args.division, args.lu))
