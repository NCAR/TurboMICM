"""Generate a complete three-stage Rosenbrock solver from a CAMP configuration.

The reader and symbolic layout are shared with the generic library's generator.
All loops emitted here traverse cells; chemistry and sparse matrix operations
are resolved at generation time. Arithmetic order within a cell is preserved.
"""
import argparse
import re
from pathlib import Path

from configure import generate, literal, mechanism_layout, read_configuration


def number(value):
    return literal(value, True)


def custom_indices(config):
    """Stable parameter insertion order, matching chem_model_finalize."""
    labels = {}
    for reaction in config['reactions']:
        if reaction['type'] == 'USER_DEFINED':
            names = [reaction['MUSICA name']]
        elif reaction['type'] == 'SURFACE':
            names = [reaction['MUSICA name'] + suffix for suffix in
                     ('.effective radius [m]', '.particle number concentration [# m-3]')]
        else:
            names = []
        for label in names:
            if label not in labels:
                labels[label] = len(labels)
    return labels


def rate_body(reaction, species, labels, index):
    """Emit the library's rate formulas with literal configuration parameters.

    Keep exp/pow and operation order identical to rates.c. The compiler can
    eliminate constant powers and pressure factors without fast-math flags.
    Derived branched/surface constants use the loader's once-computed doubles.
    """
    r, kind = reaction, reaction['type']
    defaults = dict(A=1, B=0, C=0, D=300, E=0, X=1, Y=0, a0=0, n=0,
                    k0_A=1, k0_B=0, k0_C=0, k0_D=298,
                    kinf_A=1, kinf_B=0, kinf_C=0, kinf_D=298,
                    kint_A=1, kint_B=0, kint_C=0, kint_D=298, Fc=.6, N=1)

    def p(key):
        return number(r.get(key, defaults[key]))

    arrhenius = f'{p("A")} * exp({p("C")} / t) * pow(t / {p("D")}, {p("B")}) * (1.0 + {p("E")} * pressure)'
    body = []
    if kind == 'ARRHENIUS':
        expression = arrhenius
    elif kind in ('TROE', 'TERNARY_CHEMICAL_ACTIVATION', 'TERNARY_CHEMICAL_ACTIVATION_JPL19'):
        d0 = p('k0_D') if kind.endswith('JPL19') else '300.0'
        di = p('kinf_D') if kind.endswith('JPL19') else '300.0'
        body += [f'double k0 = {p("k0_A")} * exp({p("k0_C")} / t) * pow(t / {d0}, {p("k0_B")});',
                 f'double ki = {p("kinf_A")} * exp({p("kinf_C")} / t) * pow(t / {di}, {p("kinf_B")});',
                 'double ratio = k0 * density / ki;']
        scale = '1.0' if kind == 'TERNARY_CHEMICAL_ACTIVATION' else 'density'
        expression = f'k0 * {scale} / (1.0 + ratio) * pow({p("Fc")}, {p("N")} / ({p("N")} + pow(log10(ratio), 2.0)))'
        if kind.endswith('JPL19'):
            body += [f'double falloff = {expression};',
                     f'double kint = {p("kint_A")} * exp({p("kint_C")} / t) * pow(t / {p("kint_D")}, {p("kint_B")});']
            expression = 'falloff + kint * (1.0 - falloff / ki)'
    elif kind == 'BRANCHED':
        body += [f'double a = state->model->reactions[{index}].rate.k0 * density;',
                 f'double z = state->model->reactions[{index}].rate.z;',
                 'double b = 0.43 * pow(t / 298.0, -8.0);',
                 'double v = a / (1.0 + a / b) * pow(0.41, 1.0 / (1.0 + pow(log10(a / b), 2.0)));',
                 f'double pre = {p("X")} * exp(-({p("Y")}) / t);']
        expression = 'pre * (z / (z + v))' if r['branch'] == 0 else 'pre * (v / (v + z))'
    elif kind == 'TUNNELING':
        expression = f'{p("A")} * exp(-({p("B")}) / t + {p("C")} / (t * t * t))'
    elif kind == 'TAYLOR_SERIES':
        body += ['double polynomial = 0.0, power = 1.0;']
        for coefficient in r.get('coefficients', [1]):
            body += [f'polynomial += {number(coefficient)} * power;', 'power *= t;']
        expression = f'polynomial * {arrhenius}'
    elif kind == 'REVERSIBLE':
        expression = f'{p("A")} * exp({p("C")} / t) * {number(r.get("k_r", r.get("B", 0)))}'
    elif kind == 'USER_DEFINED':
        expression = f'u[{labels[r["MUSICA name"]]}*s+c] * {number(r.get("scaling factor", r.get("A", 1)))}'
    elif kind == 'SURFACE':
        label = r['MUSICA name']
        radius = labels[label + '.effective radius [m]']
        count = labels[label + '.particle number concentration [# m-3]']
        sp = species[r['gas-phase reactant']]
        body += [f'double radius = u[{radius}*s+c], count = u[{count}*s+c];',
                 f'double speed = sqrt(state->model->reactions[{index}].rate.mean_free_speed_factor * t);']
        expression = (f'4.0 * count * CHEM_PI * radius * radius / '
                      f'(radius / {number(sp["diffusion coefficient [m2 s-1]"])} + '
                      f'4.0 / (speed * {number(r.get("reaction probability", 1))}))')
    else:
        raise ValueError(f'unsupported rate type: {kind}')
    body += [f'double value = {expression};']
    for reactant in r['reactants']:
        if species[reactant].get('parameterization'):
            body += ['value *= density;']
    body += [f'k[{index}*s+c] = value;']
    return body


def emit_rates(config, name, chunk):
    species = {s['name']: s for s in config['species']}
    labels = custom_indices(config)
    blocks = [rate_body(r, species, labels, i) for i, r in enumerate(config['reactions'])]
    lines, calls = [], []
    for start in range(0, len(blocks), chunk):
        fn = f'{name}_rates_{start//chunk}'
        calls.append(f'  {fn}(state);')
        lines += [f'static void {fn}(ChemState *state) {{',
                  '  size_t s = state->stride;',
                  '  double *restrict k = state->rates;',
                  '  const double *restrict u = state->custom;',
                  '  (void)u;',
                  '  for (size_t c = 0; c < state->cells; ++c) {',
                  '    double t = state->conditions[c].temperature;',
                  '    double pressure = state->conditions[c].pressure;',
                  '    double density = state->conditions[c].air_density;',
                  '    (void)t; (void)pressure; (void)density;']
        for body in blocks[start:start+chunk]:
            lines += ['    {'] + ['      '+line for line in body] + ['    }']
        lines += ['  }', '}']
    return '\n'.join(lines + [f'void chem_rates_{name}(ChemState *state) {{'] + calls + ['}'])


def emit_linear(name, n, indices, chunk):
    """Unroll sparse Mozart elimination and forward/back substitution.

    A chunk contains successive pivots (or rows), fused inside one cell loop.
    This retains every multiply/subtract and division in reference order.
    Functions bound compiler work without imposing a cell vector length.
    """
    factor, forward, backward = [], [], []
    for i in range(n):
        lower = [j for j in range(i+1, n) if (j, i) in indices]
        upper = [j for j in range(i+1, n) if (i, j) in indices]
        block = []
        if lower:
            block += [f'double inverse = 1.0 / a[{indices[i,i]}*s+c];']
            block += [f'a[{indices[j,i]}*s+c] *= inverse;' for j in lower]
            block += [f'a[{indices[j,k]}*s+c] -= a[{indices[j,i]}*s+c] * a[{indices[i,k]}*s+c];'
                      for k in upper for j in lower]
        factor.append(block)
        forward.append([f'x[{i}*s+c] -= a[{indices[i,j]}*s+c] * x[{j}*s+c];'
                        for j in range(i) if (i, j) in indices])
    for i in reversed(range(n)):
        backward.append([f'x[{i}*s+c] -= a[{indices[i,j]}*s+c] * x[{j}*s+c];'
                         for j in range(i+1, n) if (i, j) in indices] +
                        [f'x[{i}*s+c] /= a[{indices[i,i]}*s+c];'])
    lines, wrappers = [], []
    for family, blocks in [('factor', factor), ('forward', forward), ('backward', backward)]:
        calls = []
        args = 'double *restrict a' if family == 'factor' else 'const double *restrict a, double *restrict x'
        actual = 'nc, s, a' if family == 'factor' else 'nc, s, a, x'
        for start in range(0, len(blocks), chunk):
            selected = [b for b in blocks[start:start+chunk] if b]
            if not selected:
                continue
            fn = f'{name}_{family}_{start//chunk}'
            calls.append(f'  {fn}({actual});')
            lines += [f'static void {fn}(size_t nc, size_t s, {args}) {{',
                      '  #pragma omp simd', '  for (size_t c = 0; c < nc; ++c) {']
            for body in selected:
                lines += ['    {'] + ['      '+line for line in body] + ['    }']
            lines += ['  }', '}']
        wrappers += [f'static void {name}_{family}(size_t nc, size_t s, {args}) {{',
                     '  (void)nc; (void)s; (void)a;']
        if family != 'factor':
            wrappers += ['  (void)x;']
        wrappers += calls + ['}']
    lines += wrappers
    lines += [f'void chem_factor_{name}(const ChemModel *m, size_t nc, size_t s, double *a) {{',
              f'  (void)m; {name}_factor(nc, s, a);', '}',
              f'void chem_linear_{name}(const ChemModel *m, size_t nc, size_t s, const double *a, double *x) {{',
              f'  (void)m; {name}_forward(nc, s, a, x); {name}_backward(nc, s, a, x);', '}']
    return '\n'.join(lines)


def specialize(config, name, kernel_chunk=32, linear_chunk=8, fuse_kernels=False,
               linear='unrolled', fuse_stages=True, specialize_rates=True):
    if linear_chunk < 1:
        raise ValueError('linear chunk must be positive')
    vmap, indices = mechanism_layout(config)
    n = len(vmap)
    if n == 0:
        raise ValueError('a solver needs at least one variable species')
    lines = [generate(config, name, kernel_chunk, fuse_kernels),
             '#include <math.h>', '#include <float.h>', '#include <string.h>']
    if specialize_rates:
        lines += [emit_rates(config, name, kernel_chunk)]
    else:
        lines += [f'void chem_rates_{name}(ChemState *s) {{ chem_calculate_rates(s); }}']
    if linear == 'unrolled':
        lines += [emit_linear(name, n, indices, linear_chunk)]
    else:
        lines += [f'void chem_factor_{name}(const ChemModel *m, size_t n, size_t s, double *a) {{ chem_lu_mozart(m,n,s,a); }}',
                  f'void chem_linear_{name}(const ChemModel *m, size_t n, size_t s, const double *a, double *x) {{ chem_linear_solve(m,n,s,a,x); }}']
    if fuse_stages:
        lines += [f'static void {name}_shift(ChemState *state, double alpha) {{',
                  '  size_t s = state->stride; double *restrict a = state->jacobian;',
                  '  #pragma omp simd', '  for (size_t c = 0; c < state->cells; ++c) {']
        lines += [f'    a[{indices[i,i]}*s+c] += alpha;' for i in range(n)]
        lines += ['  }', '}']
        template = Path(__file__).with_name('rosenbrock3.c.in').read_text()
        lines += [template.replace('@NAME@', name).replace('@N@', str(n)).replace('@NNZ@', str(len(indices)))]
    else:
        # Used for measured ablation of stage fusion; same generated chemistry/LU.
        lines += [f'ChemStatus chem_solve_{name}(double dt, ChemState *s, ChemStats *stats) {{',
                  '  ChemSolver solver = chem_solver(3);',
                  f'  solver.factor = chem_factor_{name}; solver.solve = chem_linear_{name};',
                  '  return chem_solve(&solver,dt,s,stats);', '}']
    return '\n'.join(lines) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config')
    parser.add_argument('output')
    parser.add_argument('--name', required=True)
    parser.add_argument('--kernel-chunk', type=int, default=32)
    parser.add_argument('--linear-chunk', type=int, default=8)
    parser.add_argument('--fuse-kernels', type=int, choices=(0, 1), default=0)
    parser.add_argument('--linear', choices=('unrolled', 'generic'), default='unrolled')
    parser.add_argument('--fuse-stages', type=int, choices=(0, 1), default=1)
    parser.add_argument('--specialize-rates', type=int, choices=(0, 1), default=1)
    args = parser.parse_args()
    if not re.fullmatch('[a-zA-Z_][a-zA-Z_0-9]*', args.name):
        parser.error('name must be a C identifier')
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(specialize(read_configuration(args.config), args.name, args.kernel_chunk,
                               args.linear_chunk, bool(args.fuse_kernels), args.linear,
                               bool(args.fuse_stages), bool(args.specialize_rates)))


if __name__ == '__main__':
    main()
