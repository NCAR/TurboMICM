"""Generate a complete three-stage Rosenbrock solver from a CAMP configuration.

The reader and symbolic layout are shared with the generic library's generator.
All loops emitted here traverse cells; chemistry and sparse matrix operations
are resolved at generation time. Arithmetic order within a cell is preserved.
"""
import argparse
import re
from pathlib import Path

from configure import generate, literal, mechanism_layout, read_configuration, write_dependencies


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
    return '\n'.join(lines + [f'void chem_rates_{name}(ChemState *state) {{', '  (void)state;'] + calls + ['}'])


def linear_plan(n, indices):
    """Tuples name an operation and its fixed array indices, in execution order."""
    factor, forward, backward = [], [], []
    for i in range(n):
        lower = [j for j in range(i+1, n) if (j, i) in indices]
        upper = [j for j in range(i+1, n) if (i, j) in indices]
        block = []
        if lower:
            block += [('inverse', indices[i,i])]
            block += [('scale', indices[j,i]) for j in lower]
            block += [('update', indices[j,k], indices[j,i], indices[i,k]) for k in upper for j in lower]
        factor.append(block)
        forward.append([('subtract', i, indices[i,j], j) for j in range(i) if (i, j) in indices])
    for i in reversed(range(n)):
        backward.append([('subtract', i, indices[i,j], j) for j in range(i+1, n) if (i, j) in indices] +
                        [('divide', i, indices[i,i])])
    return [('factor', factor), ('forward', forward), ('backward', backward)]


def linear_statement(op):
    kind, i, *rest = op
    if kind == 'inverse':
        return f'double inverse = 1.0 / a[{i}*s+c];'
    if kind == 'scale':
        return f'a[{i}*s+c] *= inverse;'
    if kind == 'update':
        j, k = rest
        return f'a[{i}*s+c] -= a[{j}*s+c] * a[{k}*s+c];'
    if kind == 'subtract':
        j, k = rest
        return f'x[{i}*s+c] -= a[{j}*s+c] * x[{k}*s+c];'
    j, = rest
    return f'x[{i}*s+c] /= a[{j}*s+c];'


def linear_operation_loop(op, pointers):
    """Small SIMD loop; restrict also describes distinct sparse element arrays."""
    kind, i, *rest = op
    statement = linear_statement(op).replace('double inverse =', 'inverse[c] =').replace('*= inverse;', '*= inverse[c];')
    declarations = []
    if pointers:
        if kind == 'inverse':
            declarations = [f'const double *restrict diagonal = a+{i}*s;']
            statement = 'inverse[c] = 1.0 / diagonal[c];'
        elif kind == 'scale':
            declarations = [f'double *restrict target = a+{i}*s;']
            statement = 'target[c] *= inverse[c];'
        elif kind == 'update':
            j, k = rest
            declarations = [f'double *restrict target = a+{i}*s;',
                            f'const double *restrict left = a+{j}*s, *restrict right = a+{k}*s;']
            statement = 'target[c] -= left[c]*right[c];'
        elif kind == 'subtract':
            j, k = rest
            declarations = [f'double *restrict target = x+{i}*s;',
                            f'const double *restrict left = a+{j}*s, *restrict right = x+{k}*s;']
            statement = 'target[c] -= left[c]*right[c];'
        else:
            j, = rest
            declarations = [f'double *restrict target = x+{i}*s;', f'const double *restrict diagonal = a+{j}*s;']
            statement = 'target[c] /= diagonal[c];'
    return ['  {'] + ['    '+line for line in declarations] + [
        '    #pragma omp simd', '    for (size_t c = 0; c < nc; ++c) {', '      '+statement, '    }', '  }']


def dot_factor_blocks(n, indices):
    """Left-looking exact LU: accumulate each element once, in pivot order.

    Each element sees the same subtractions as Mozart, followed by the same
    column reciprocal multiplication. Completed L/U entries are the operands.
    The sparse elimination tree is unchanged; only independent work moves.
    """
    blocks = []
    for i in range(n):
        block = []
        upper = [(i,k,False) for k in range(i,n) if (i,k) in indices]
        lower = [(j,i,True) for j in range(i+1,n) if (j,i) in indices]
        for row, col, scale in upper+lower:
            if scale and col == i and row == lower[0][0]:
                block += linear_operation_loop(('inverse',indices[i,i]), True)
            target = indices[row,col]
            body = [f'double value = target[c];']
            for p in range(i):
                if (row,p) in indices and (p,col) in indices:
                    body += [f'value -= a[{indices[row,p]}*s+c] * a[{indices[p,col]}*s+c];']
            body += ['target[c] = value * inverse[c];' if scale else 'target[c] = value;']
            if len(body) == 2 and not scale:
                continue  # unchanged U entry, including an isolated diagonal
            block += ['  {', f'    double *restrict target = a+{target}*s;',
                      '    #pragma omp simd', '    for (size_t c = 0; c < nc; ++c) {']
            block += ['      '+line for line in body] + ['    }', '  }']
        blocks.append(block)
    return blocks


def dot_solve_block(body):
    """Keep a triangular row's partial sum in a SIMD register."""
    i = body[0][1]
    statements = ['double value = target[c];']
    for op in body:
        if op[0] == 'subtract':
            _, _, j, k = op
            statements += [f'value -= a[{j}*s+c] * x[{k}*s+c];']
        else:
            statements += [f'value /= a[{op[2]}*s+c];']
    statements += ['target[c] = value;']
    return ['  {', f'    double *restrict target = x+{i}*s;',
            '    #pragma omp simd', '    for (size_t c = 0; c < nc; ++c) {'] + [
                '      '+line for line in statements] + ['    }','  }']


def emit_linear(name, n, indices, chunk, loop, cell_stride=0):
    """Unroll sparse Mozart elimination and forward/back substitution.

    A chunk contains successive pivots (or rows), fused inside one cell loop.
    This retains every multiply/subtract and division in reference order.
    Functions bound compiler work without imposing a cell vector length.
    """
    lines, wrappers = [], []
    for family, blocks in linear_plan(n, indices):
        if family == 'factor' and loop == 'dot':
            blocks = dot_factor_blocks(n,indices)
        calls = []
        args = 'double *restrict a, double *restrict inverse' if family == 'factor' else 'const double *restrict a, double *restrict x'
        actual = 'nc, s, a, inverse' if family == 'factor' else 'nc, s, a, x'
        for start in range(0, len(blocks), chunk):
            selected = [b for b in blocks[start:start+chunk] if b]
            if not selected:
                continue
            fn = f'{name}_{family}_{start//chunk}'
            calls.append(f'  {fn}({actual});')
            lines += [f'static void {fn}(size_t nc, size_t s, {args}) {{']
            if cell_stride:
                lines[-1] = lines[-1].replace('size_t nc, size_t s,', 'size_t unused_nc, size_t unused_s,')
                lines += [f'  const size_t nc = {cell_stride}, s = {cell_stride};', '  (void)unused_nc; (void)unused_s;']
            if family == 'factor':
                lines += ['  (void)inverse;']
            if loop == 'dot':
                for body in selected:
                    lines += body if family == 'factor' else dot_solve_block(body)
            elif loop in ('operations', 'pointers'):
                for body in selected:
                    for op in body:
                        lines += linear_operation_loop(op,loop == 'pointers')
            else:
                lines += ['  #pragma omp simd', '  for (size_t c = 0; c < nc; ++c) {']
                for body in selected:
                    lines += ['    {'] + ['      '+linear_statement(op) for op in body] + ['    }']
                lines += ['  }']
            lines += ['}']
        wrappers += [f'static void {name}_{family}(size_t nc, size_t s, {args}) {{',
                     '  (void)nc; (void)s; (void)a;']
        if family != 'factor':
            wrappers += ['  (void)x;']
        else:
            wrappers += ['  (void)inverse;']
        wrappers += calls + ['}']
    lines += wrappers
    lines += [f'void chem_factor_{name}(const ChemModel *m, size_t nc, size_t s, double *a) {{',
              '  double *inverse = malloc(nc*sizeof(double));',
              '  if (!inverse) { for (size_t c = 0; c < nc; ++c) a[m->diagonal[0]*s+c] = NAN; return; }',
              f'  {name}_factor(nc, s, a, inverse); free(inverse);', '}',
              f'void chem_linear_{name}(const ChemModel *m, size_t nc, size_t s, const double *a, double *x) {{',
              f'  (void)m; {name}_forward(nc, s, a, x); {name}_backward(nc, s, a, x);', '}']
    return '\n'.join(lines)


def specialize(config, name, kernel_chunk=32, linear_chunk=1, fuse_kernels=True,
               linear='unrolled', fuse_stages=True, specialize_rates=True, linear_loop='rows',
               vector_error=True, cell_stride=4, initialize=True, solver_tile=4):
    if linear_chunk < 1:
        raise ValueError('linear chunk must be positive')
    if cell_stride < 0:
        raise ValueError('cell stride must be nonnegative')
    if solver_tile < 0:
        raise ValueError('solver tile must be nonnegative')
    vmap, indices = mechanism_layout(config)
    n = len(vmap)
    if n == 0:
        raise ValueError('a solver needs at least one variable species')
    lines = ['/* Generated by tools/specialize.py; edit the generator or configuration. */']
    if cell_stride:
        lines += [generate(config,name+'_fixed',kernel_chunk,fuse_kernels,cell_stride,False,initialize=initialize)]
    lines += [generate(config, name, kernel_chunk, fuse_kernels,dispatch_stride=cell_stride,initialize=initialize),
             '#include <math.h>', '#include <float.h>', '#include <string.h>']
    if specialize_rates:
        lines += [emit_rates(config, name, kernel_chunk)]
    else:
        lines += [f'void chem_rates_{name}(ChemState *s) {{ chem_calculate_rates(s); }}']
    if linear == 'unrolled':
        if cell_stride:
            lines += [emit_linear(name+'_dynamic', n, indices, linear_chunk, linear_loop),
                      emit_linear(name+'_fixed', n, indices, linear_chunk, linear_loop,cell_stride),
                      f'void chem_factor_{name}(const ChemModel *m, size_t nc, size_t s, double *a) {{',
                      f'  if (nc == {cell_stride} && s == {cell_stride}) chem_factor_{name}_fixed(m,nc,s,a);',
                      f'  else chem_factor_{name}_dynamic(m,nc,s,a);', '}',
                      f'void chem_linear_{name}(const ChemModel *m, size_t nc, size_t s, const double *a, double *x) {{',
                      f'  if (nc == {cell_stride} && s == {cell_stride}) chem_linear_{name}_fixed(m,nc,s,a,x);',
                      f'  else chem_linear_{name}_dynamic(m,nc,s,a,x);', '}']
        else:
            lines += [emit_linear(name, n, indices, linear_chunk, linear_loop)]
    else:
        lines += [f'void chem_factor_{name}(const ChemModel *m, size_t n, size_t s, double *a) {{ chem_lu_mozart(m,n,s,a); }}',
                  f'void chem_linear_{name}(const ChemModel *m, size_t n, size_t s, const double *a, double *x) {{ chem_linear_solve(m,n,s,a,x); }}']
    if fuse_stages:
        lines += [f'static void {name}_shift(ChemState *state, double alpha) {{',
                  '  size_t s = state->stride; double *restrict a = state->jacobian;',
                  '  #pragma omp simd', '  for (size_t c = 0; c < state->cells; ++c) {']
        lines += [f'    a[{indices[i,i]}*s+c] += alpha;' for i in range(n)]
        lines += ['  }', '}']
        if solver_tile:
            template = Path(__file__).with_name('rosenbrock3_blocked.c.in').read_text()
            if linear == 'unrolled':
                if cell_stride:
                    factor = (f'if (nc == {cell_stride} && s == {cell_stride}) '
                              f'{name}_fixed_factor(nc,s,q->jacobian,q->error); '
                              f'else {name}_dynamic_factor(nc,s,q->jacobian,q->error);')
                else:
                    factor = f'{name}_factor(nc,s,q->jacobian,q->error);'
            else:
                factor = f'chem_factor_{name}(m,nc,s,q->jacobian);'
            template = template.replace('@FACTOR_CALL@',factor).replace('@TILE@',str(solver_tile)).replace('@R@',str(len(config['reactions'])))
        else:
            template = Path(__file__).with_name('rosenbrock3.c.in').read_text()
        lines += [template.replace('@NAME@', name).replace('@N@', str(n)).replace('@NNZ@', str(len(indices)))
                  .replace('@VECTOR_ERROR@', str(int(vector_error)))]
        lines[-1] = lines[-1].replace('@INITIALIZE@',str(int(initialize)))
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
    parser.add_argument('--linear-chunk', type=int, default=1)
    parser.add_argument('--fuse-kernels', type=int, choices=(0, 1), default=1)
    parser.add_argument('--linear', choices=('unrolled', 'generic'), default='unrolled')
    parser.add_argument('--linear-loop', choices=('operations', 'pointers', 'rows', 'dot'), default='rows')
    parser.add_argument('--fuse-stages', type=int, choices=(0, 1), default=1)
    parser.add_argument('--specialize-rates', type=int, choices=(0, 1), default=1)
    parser.add_argument('--vector-error', type=int, choices=(0, 1), default=1)
    parser.add_argument('--cell-stride', type=int, default=4)
    parser.add_argument('--initialize-kernels', type=int, choices=(0,1), default=1)
    parser.add_argument('--solver-tile', type=int, default=4)
    args = parser.parse_args()
    if not re.fullmatch('[a-zA-Z_][a-zA-Z_0-9]*', args.name):
        parser.error('name must be a C identifier')
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    config = read_configuration(args.config)
    path.write_text(specialize(config, args.name, args.kernel_chunk,
                               args.linear_chunk, bool(args.fuse_kernels), args.linear,
                               bool(args.fuse_stages), bool(args.specialize_rates), args.linear_loop,
                               bool(args.vector_error), args.cell_stride, bool(args.initialize_kernels), args.solver_tile))
    write_dependencies(path,config['sources'])


if __name__ == '__main__':
    main()
