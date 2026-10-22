"""Compile a validated CAMP mechanism into fixed-index C11 chemistry kernels.

The IR retains ordered reactant occurrences and products: neither summing
duplicate stoichiometry nor reassociating an accumulation is numerically free.
Symbolic elimination mirrors the generic library, including Markowitz ties.
Only the cell dimension remains a loop in the emitted arithmetic.
"""
from __future__ import annotations

import argparse
import copy
from pathlib import Path

from configure import emit_c, read_mechanism, literal, load_document

RATE_ORDER = "ARRHENIUS TROE TERNARY TERNARY_JPL19 BRANCHED TUNNELING TAYLOR REVERSIBLE USER SURFACE".split()


def real(x):
    """A double literal, including for integer-valued configuration fields."""
    s = literal(x)
    return s if "." in s or "e" in s else s + ".0"


class Plan:
    def __init__(self, mechanism):
        self.mechanism = mechanism
        species = [s["name"] for s in mechanism["species"] if not s.get("parameterization")]
        self.n = n = len(species)
        reactions = sorted(copy.deepcopy(mechanism["reactions"]), key=lambda r: RATE_ORDER.index(r["type"]))
        self.reactions = reactions
        original = {name: i for i, name in enumerate(species)}
        pattern = [[False] * n for _ in range(n)]
        custom = 0
        for r in reactions:
            r["parameters"] = len([s for s in r["reactants"] if s not in original])
            r["lhs"] = [original[s] for s in r["reactants"] if s in original]
            r["rhs"] = [(original[s], y) for s, y in r["products"] if s in original]
            r["custom"] = custom
            custom += 1 if r["type"] == "USER" else 2 if r["type"] == "SURFACE" else 0
            for col in r["lhs"]:
                for row in r["lhs"] + [s for s, _ in r["rhs"]]:
                    pattern[row][col] = True
        # MICM uses size_t arithmetic for empty rows/columns (also on 32-bit C).
        # All negative costs are -k; modulo ordering makes them larger than
        # nonnegative costs, and orders them by -k within that class.
        work = copy.deepcopy(pattern)
        alive = set(range(n))
        permutation = []
        while alive:
            def cost(v):
                value = (sum(work[v][j] for j in alive) - 1) * (sum(work[j][v] for j in alive) - 1)
                return value < 0, value, v
            pivot = min(alive, key=cost)
            permutation.append(pivot)
            alive.remove(pivot)
            for row in alive:
                if work[row][pivot]:
                    for col in alive:
                        work[row][col] |= work[pivot][col]
        inverse = {v: i for i, v in enumerate(permutation)}
        self.names = [species[i] for i in permutation]
        for r in reactions:
            r["lhs"] = [inverse[v] for v in r["lhs"]]
            r["rhs"] = [(inverse[v], y) for v, y in r["rhs"]]
        fill = [[pattern[i][j] for j in permutation] for i in permutation]
        for i in range(n):
            fill[i][i] = True
        for i in range(n):
            for col in range(i + 1, n):
                if fill[i][col]:
                    for row in range(i + 1, n):
                        fill[row][col] |= fill[row][i]
        self.index = {(r, c): i for i, (r, c) in enumerate(
            (r, c) for r in range(n) for c in range(n) if fill[r][c])}
        self.nnz = len(self.index)


class Emitter:
    def __init__(self, chunk, fusion, width):
        self.lines = []
        self.chunk = chunk
        self.fusion = fusion
        self.width = width

    def put(self, *lines):
        self.lines.extend(lines)

    def loop(self, statements, fixed=True):
        stride = str(self.width) if fixed and self.width else 'n'
        self.put('    #pragma omp simd', f'    for (size_t c=0;c<{stride};++c) {{')
        address_stride = f'((size_t){self.width})' if fixed and self.width else 'n'
        self.lines.extend('        ' + s.replace('*n', '*'+address_stride) for s in statements)
        self.put('    }')

    def blocks(self, name, signature, arguments, setup, operations):
        """Bound function size without prescribing the compiler's inlining."""
        parts = [operations[i:i+self.chunk] for i in range(0, len(operations), self.chunk)]
        for part, ops in enumerate(parts):
            self.put(f'static void {name}_{part}({signature})', '{', *setup)
            for start in range(0, len(ops), self.fusion):
                statements = []
                for op in ops[start:start+self.fusion]:
                    statements.extend(['{', *op, '}'])
                self.loop(statements)
            self.put('}')
        self.put(f'static void {name}({signature})', '{')
        if not parts:
            self.put(*['    (void)' + a.strip() + ';' for a in arguments.split(',')])
        for part in range(len(parts)):
            self.put(f'    {name}_{part}({arguments});')
        self.put('}')


def rate_statements(r, i):
    """Emit the same ordered double expressions as rates.c with fixed inputs."""
    p = {k: real(v) for k, v in r['params'].items()}
    kind = r['type']
    out = []

    def arr(prefix='', denom='300.0'):
        return f"{p[prefix+'A']} * exp({p[prefix+'C']} / t) * pow(t / {p.get(prefix+'D', denom)}, {p[prefix+'B']})"

    if kind in ('ARRHENIUS', 'TAYLOR'):
        expression = arr() + f" * (1.0 + {p['E']} * pressure)"
        if kind == 'TAYLOR':
            out += ['double poly=0.0, power=1.0;']
            for coefficient in r['coefficients']:
                out += [f'poly += {real(coefficient)} * power;', 'power *= t;']
            expression = 'poly * ' + expression
    elif kind in ('TROE', 'TERNARY', 'TERNARY_JPL19'):
        out += [f'double k0={arr("k0_")};', f'double ki={arr("kinf_")};']
        if kind == 'TERNARY_JPL19':
            out += [f'double km={arr("kint_")};', 'double k0m=k0*density, ratio=k0m/ki;']
            out += [f'double kf=k0m/(1.0+ratio)*pow({p["Fc"]}, {p["N"]}/({p["N"]}+pow(log10(ratio),2.0)));']
            expression = 'kf+km*(1.0-kf/ki)'
        else:
            out += ['double ratio=k0*density/ki;']
            scale = 'density' if kind == 'TROE' else '1.0'
            expression = f'k0*{scale}/(1.0+ratio)*pow({p["Fc"]},{p["N"]}/({p["N"]}+pow(log10(ratio),2.0)))'
    elif kind == 'BRANCHED':
        # Prepared once by the builder, using the platform's double libm.
        # Keeping these two derived constants avoids cross-libm folding drift.
        out += [f'const ChemBranched *bpar=&s->mechanism->processes[{i}].rate.data.branched;',
                'double a=bpar->k0*density, b=.43*pow(t/298.0,-8.0);',
                'double A=a/(1.0+a/b)*pow(.41,1.0/(1.0+pow(log10(a/b),2.0)));',
                f'double pre={p["X"]}*exp(-({p["Y"]})/t);']
        expression = 'pre*(bpar->z/(bpar->z+A))' if r['params']['branch'] == 0 else 'pre*(A/(A+bpar->z))'
    elif kind == 'TUNNELING':
        expression = f'{p["A"]}*exp(-({p["B"]})/t+{p["C"]}/(t*t*t))'
    elif kind == 'REVERSIBLE':
        expression = f'{p["A"]}*exp({p["C"]}/t)*{p["reverse_rate"]}'
    elif kind == 'USER':
        expression = f's->custom[{r["custom"]}*n+c]*{p["scaling"]}'
    elif kind == 'SURFACE':
        out += [f'double speed=sqrt(s->mechanism->processes[{i}].rate.data.surface.speed_factor*t);',
                f'double radius=s->custom[{r["custom"]}*n+c], count=s->custom[{r["custom"]+1}*n+c];']
        expression = f'4.0*count*CHEM_PI*radius*radius/(radius/{p["diffusion"]}+4.0/(speed*{p["probability"]}))'
    else:
        raise ValueError(f'No rate emitter for {kind}')
    out += [f'double rate={expression};']
    if r['parameters']:
        out += ['double mult=1.0;'] + ['mult *= density;'] * r['parameters'] + ['rate *= mult;']
    return out + [f's->rates[{i}*n+c]=rate;']


def emit_rates(e, plan, symbol, width, cache):
    ops = []
    for i, r in enumerate(plan.reactions):
        ops.append(['double t=s->conditions[c].temperature, pressure=s->conditions[c].pressure, density=s->conditions[c].air_density;',
                    '(void)t; (void)pressure; (void)density;'] + rate_statements(r, i))
    e.blocks('rates', 'ChemState *s', 's', ['    size_t n=s->cells; (void)n;'], ops)
    e.put(f'void {symbol}_update_rates(ChemState *s)', '{')
    if cache:
        custom_count = sum(1 if r['type']=='USER' else 2 if r['type']=='SURFACE' else 0 for r in plan.reactions)
        e.put('    size_t cells=s->cells;',
              f'    size_t custom_count={custom_count}*cells;',
              '    double *saved=s->cached_rate_inputs;',
              '    if (saved && s->rates_cached) {',
              '        int same=memcmp(saved+3*cells,s->custom,custom_count*sizeof(double))==0;',
              '        for (size_t c=0;same && c<cells;++c)',
              '            same=memcmp(saved+3*c,&s->conditions[c].temperature,sizeof(double))==0 && memcmp(saved+3*c+1,&s->conditions[c].pressure,sizeof(double))==0 && memcmp(saved+3*c+2,&s->conditions[c].air_density,sizeof(double))==0;',
              '        if (same) return;',
              '    }')
    if width:
        e.put(f'    if (s->cells!={width}) chem_update_rates(s); else rates(s);')
    else:
        e.put('    rates(s);')
    if cache:
        e.put('    if (!saved && cells<=SIZE_MAX/sizeof(double)/(3+'+str(custom_count)+'))',
              '        s->cached_rate_inputs=saved=malloc((3*cells+custom_count)*sizeof(double));',
              '    if (saved) {',
              '        for (size_t c=0;c<cells;++c) {',
              '            saved[3*c]=s->conditions[c].temperature; saved[3*c+1]=s->conditions[c].pressure; saved[3*c+2]=s->conditions[c].air_density;',
              '        }',
              '        memcpy(saved+3*cells,s->custom,custom_count*sizeof(double));',
              '        s->rates_cached=1;',
              '    }')
    e.put('}')


def emit_chemistry(e, plan, width, joint):
    forcing, jacobian = [], []
    for i, r in enumerate(plan.reactions):
        op = [f'double rate=k[{i}*n+c];'] + [f'rate *= y[{v}*n+c];' for v in r['lhs']]
        op += [f'f[{v}*n+c] -= rate;' for v in r['lhs']]
        op += [f'f[{v}*n+c] += {real(y)}*rate;' for v, y in r['rhs']]
        op += ['(void)rate;']
        forcing.append(op)
    for col in range(plan.n):
        for i, r in enumerate(plan.reactions):
            for independent in r['lhs']:
                if independent != col:
                    continue
                dependent = r['lhs'].copy()
                dependent.remove(col)  # first occurrence, exactly as the generic plan
                op = [f'double rate=k[{i}*n+c];'] + [f'rate *= y[{v}*n+c];' for v in dependent]
                op += [f'f[{plan.index[v,col]}*n+c] += rate;' for v in dependent + [col]]
                op += [f'f[{plan.index[v,col]}*n+c] -= {real(y)}*rate;' for v, y in r['rhs']]
                jacobian.append(op)
    for name, ops, count in [('forcing', forcing, plan.n), ('jacobian', jacobian, plan.nnz)]:
        e.blocks(name+'_add', 'size_t n, const double *restrict k, const double *restrict y, double *restrict f',
                 'n,k,y,f', ['    (void)y; (void)n;'], ops)
        e.put(f'static void {name}(const ChemState *s,const double *y,double *f)', '{')
        if width and name == 'forcing':
            e.put(f'    if (s->cells!={width}) {{ chem_forcing(s,y,f); return; }}')
        e.put(f'    memset(f,0,{count}*s->cells*sizeof(double));',
              f'    {name}_add(s->cells,s->rates,y,f);', '}')

    if joint:
        # Different Jacobian columns are disjoint. Traversing reactions first
        # preserves every individual entry's accumulation order, while exposing
        # common rate products to the C compiler in the forcing/Jacobian pass.
        operations = []
        for reaction, (r, force) in enumerate(zip(plan.reactions, forcing)):
            op = ['{', *force, '}']
            for col in r['lhs']:
                dependent = r['lhs'].copy()
                dependent.remove(col)
                derivative = [f'double rate=k[{reaction}*n+c];']
                derivative += [f'rate *= y[{v}*n+c];' for v in dependent]
                derivative += [f'j[{plan.index[v,col]}*n+c] += rate;' for v in dependent+[col]]
                derivative += [f'j[{plan.index[v,col]}*n+c] -= {real(y)}*rate;' for v, y in r['rhs']]
                op.extend(['{', *derivative, '}'])
            operations.append(op)
        e.blocks('initial_add', 'size_t n,const double *restrict k,const double *restrict y,double *restrict f,double *restrict j',
                 'n,k,y,f,j', ['    (void)n; (void)y; (void)j;'], operations)
        e.put('static void initial(ChemSolver *s,const ChemState *state)', '{')
        if width:
            e.put(f'    if (s->cells!={width}) {{ chem_forcing(state,state->values,s->initial); chem_solver_negative_jacobian(s,state,state->values); return; }}')
        e.put(f'    memset(s->initial,0,{plan.n}*s->cells*sizeof(double));',
              f'    memset(s->matrix,0,{plan.nnz}*s->cells*sizeof(double));',
              '    initial_add(s->cells,state->rates,state->values,s->initial,s->matrix);', '}')
        e.put('#define CHEM_INITIAL(s,state) initial(s,state)')


def emit_factor(e, plan, width, chunk, fused):
    index, n = plan.index, plan.n
    factor = []
    for pivot in range(n):
        diag = index[pivot,pivot]
        # Return cannot appear in an omp simd structured block. A reduction
        # checks each pivot below, immediately before its arithmetic is emitted.
        op = [f'double inverse=1.0/v[{diag}*n+c];']
        op += [f'v[{index[row,pivot]}*n+c] *= inverse;' for row in range(pivot+1,n) if (row,pivot) in index]
        for col in range(pivot+1,n):
            if (pivot,col) in index:
                for row in range(pivot+1,n):
                    if (row,pivot) in index:
                        op += [f'v[{index[row,col]}*n+c] -= v[{index[row,pivot]}*n+c]*v[{index[pivot,col]}*n+c];']
        factor.append(op)
    # One chunk per subset of pivots. The singular check uses a SIMD reduction,
    # then returns outside that loop. The normal path is all straight-line math.
    for part, start in enumerate(range(0, n, chunk)):
        e.put(f'static ChemStatus factor_{part}(size_t n,double *restrict v)', '{', '    (void)n;')
        for pivot in range(start,min(n,start+chunk)):
            diag = index[pivot,pivot]
            stride = width or 'n'
            e.put('    { int bad=0;', '    #pragma omp simd reduction(|:bad)',
                  f'    for (size_t c=0;c<{stride};++c) bad |= v[((size_t){diag})*{stride}+c]==0 || !isfinite(v[((size_t){diag})*{stride}+c]);',
                  '    if (bad) return CHEM_SINGULAR;', '    }')
            if fused:
                if len(factor[pivot]) > 1:
                    e.loop(factor[pivot])
            else:
                for statement in factor[pivot][1:]:
                    e.loop([factor[pivot][0], statement] if 'inverse' in statement else [statement])
        e.put('    return CHEM_OK;', '}')
    e.put('static ChemStatus factor(const ChemSparse *a,size_t n,double *v,void *ctx)', '{', '    (void)a; (void)ctx; (void)n; (void)v;')
    if width:
        e.put(f'    if (n!={width}) return chem_lu_factor(a,n,v,CHEM_LU_MOZART);')
    for part in range((n+chunk-1)//chunk):
        e.put(f'    {{ ChemStatus status=factor_{part}(n,v); if (status) return status; }}')
    e.put('    return CHEM_OK;', '}')


def emit_solve(e, plan, symbol, width):
    index, n = plan.index, plan.n
    ops = []
    for row in range(n):
        op = [f'double value=x[{row}*n+c];']
        op += [f'value -= v[{index[row,col]}*n+c]*x[{col}*n+c];' for col in range(row) if (row,col) in index]
        if len(op)>1:
            ops.append(op + [f'x[{row}*n+c]=value;'])
    for row in reversed(range(n)):
        op = [f'double value=x[{row}*n+c];']
        op += [f'value -= v[{index[row,col]}*n+c]*x[{col}*n+c];' for col in range(row+1,n) if (row,col) in index]
        ops.append(op + [f'x[{row}*n+c]=value/v[{index[row,row]}*n+c];'])
    e.blocks('triangular', 'size_t n,const double *restrict v,double *restrict x', 'n,v,x', ['    (void)n;'], ops)
    e.put('static void solve(const ChemSparse *a,size_t n,const double *v,double *x,void *ctx)',
          '{ (void)a; (void)ctx; ' + (f'if (n!={width}) {{ chem_lu_solve(a,n,v,x); return; }} ' if width else '') + 'triangular(n,v,x); }')
    e.put(f'void {symbol}_specialized_attach(ChemSolver *s)', '{',
          '    s->linear=(ChemLinearSolver){factor,solve,NULL};', '}')
    e.put('static void diagonal(size_t n,double alpha,double *v)', '{', '    (void)n; (void)alpha; (void)v;')
    for row in range(n):
        e.loop([f'v[{index[row,row]}*n+c] += alpha;'], fixed=False)
    e.put('}')


def emit_integrator(e, symbol, width, stages):
    if stages:
        e.put('#define CHEM_STAGE_COUNT 3', '#define CHEM_FUSED_ESTIMATE',
              '#define CHEM_ALGEBRAIC(s,v) 0', '#define chem_constraint_count(m) 0',
              '#define chem_initialize_constraints(s,state,p,stats) CHEM_OK')
    solve_name = symbol+'_solve3' if stages else symbol+'_solve'
    e.put(f'#define CHEM_SOLVE {solve_name}',
          '#define chem_forcing forcing',
          ('#define chem_solver_negative_jacobian(s,state,y) ' + (f'((s)->cells=={width} ? jacobian(state,y,(s)->matrix) : chem_solver_negative_jacobian(s,state,y))' if width else 'jacobian(state,y,(s)->matrix)')),
          '#define chem_add_mass_diagonal(m,pattern,cells,alpha,matrix) diagonal(cells,alpha,matrix)',
          '#include "rosenbrock_step.inc"')
    if stages:
        e.put(f'ChemStatus {symbol}_solve(ChemSolver *s,ChemState *state,double dt,const ChemRosenbrockParameters *p,ChemSolverStats *stats)',
              '{', '    if (!p) return CHEM_INVALID;',
              f'    if (p->stages==3) return {symbol}_solve3(s,state,dt,p,stats);',
              '    return chem_solve(s,state,dt,p,stats);', '}')


def generate(mechanism, symbol, chunk=32, fused=True, fusion=1, width=0, cache=True, stages=True, joint=False):
    plan = Plan(mechanism)
    e = Emitter(chunk, fusion, width)
    e.put(emit_c(mechanism, symbol), '#include "solver_internal.h"', '#include <math.h>', '#include <string.h>')
    # Expose a separate ordered constructor. The generic constructor remains
    # available for independent tests and for callers using the generic API.
    e.put(f'ChemStatus {symbol}_specialized_create(ChemMechanism **out)', '{',
          f'    ChemStatus status={symbol}_create(out);',
          '    if (!status) status=chem_mechanism_reorder(*out);', '    return status;', '}')
    emit_rates(e, plan, symbol, width, cache)
    emit_chemistry(e, plan, width, joint)
    emit_factor(e, plan, width, chunk, fused)
    emit_solve(e, plan, symbol, width)
    emit_integrator(e, symbol, width, stages)
    return '\n'.join(e.lines) + '\n'


def dependencies(configuration):
    """Record transitive CAMP includes for make, including newly added files."""
    seen = set()

    def visit(path):
        path = path.resolve()
        if path in seen:
            return
        seen.add(path)
        for name in load_document(path).get('camp-files', []):
            visit(path.parent / name)

    visit(configuration)
    return sorted(seen)


def make_path(path):
    import os
    text = os.path.relpath(path)
    return text.replace('\\', '\\\\').replace('$', '$$').replace('#', '\\#').replace(' ', '\\ ')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('configuration', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--symbol', default='configured')
    parser.add_argument('--depfile', type=Path, help='optional make dependency file')
    parser.add_argument('--chunk', type=int, default=32, help='operations/pivots per emitted function (positive)')
    parser.add_argument('--width', type=int, default=32, help='constant cell stride; 0 emits variable stride; other group sizes use generic kernels')
    parser.add_argument('--fusion', type=int, default=1, help='operations per SIMD cell loop (positive)')
    parser.add_argument('--joint-initial', action='store_true', help='combine initial forcing and Jacobian in one pass')
    parser.add_argument('--generic-stages', action='store_true', help='retain dynamic stages and separate estimate loops')
    parser.add_argument('--no-cache', action='store_true', help='recompute unchanged rate inputs')
    parser.add_argument('--split-lu', action='store_true', help='emit separate cell loops for each LU update')
    args = parser.parse_args()
    if args.chunk < 1 or args.fusion < 1 or args.width < 0:
        parser.error('--chunk and --fusion must be positive; --width must be nonnegative')
    result = generate(read_mechanism(args.configuration), args.symbol,
                      chunk=args.chunk, fused=not args.split_lu,
                      fusion=args.fusion, width=args.width, cache=not args.no_cache,
                      stages=not args.generic_stages, joint=args.joint_initial)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(result)
    if args.depfile:
        args.depfile.parent.mkdir(parents=True, exist_ok=True)
        args.depfile.write_text(make_path(args.output) + ': ' +
                               ' '.join(make_path(p) for p in dependencies(args.configuration)) + '\n')


if __name__ == '__main__':
    main()
