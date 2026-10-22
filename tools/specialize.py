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

from configure import emit_c, read_mechanism, literal

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
    def __init__(self, chunk, fusion):
        self.lines = []
        self.chunk = chunk
        self.fusion = fusion

    def put(self, *lines):
        self.lines.extend(lines)

    def loop(self, statements):
        self.put('    #pragma omp simd', '    for (size_t c=0;c<n;++c) {')
        self.lines.extend('        ' + s for s in statements)
        self.put('    }')

    def blocks(self, name, signature, arguments, setup, operations, returns=False):
        """Bound function size without prescribing the compiler's inlining."""
        parts = [operations[i:i+self.chunk] for i in range(0, len(operations), self.chunk)]
        for part, ops in enumerate(parts):
            self.put(f'static {"ChemStatus" if returns else "void"} {name}_{part}({signature})', '{', *setup)
            for start in range(0, len(ops), self.fusion):
                self.loop([line for op in ops[start:start+self.fusion] for line in ['{', *op, '}']])
            if returns:
                self.put('    return CHEM_OK;')
            self.put('}')
        self.put(f'static {"ChemStatus" if returns else "void"} {name}({signature})', '{')
        if not parts:
            # All signatures used here name these arguments.
            self.put(*['    (void)' + a.strip() + ';' for a in arguments.split(',')])
        for part in range(len(parts)):
            call = f'{name}_{part}({arguments})'
            self.put(f'    {{ ChemStatus status={call}; if (status) return status; }}' if returns else f'    {call};')
        if returns:
            self.put('    return CHEM_OK;')
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


def generate(mechanism, symbol, chunk=32, fused=True, fusion=1):
    plan = Plan(mechanism)
    e = Emitter(chunk, fusion)
    e.put(emit_c(mechanism, symbol), '#include "solver_internal.h"', '#include <math.h>', '#include <string.h>')
    # Expose a separate ordered constructor. The generic constructor remains
    # available for independent tests and for callers using the generic API.
    e.put(f'ChemStatus {symbol}_specialized_create(ChemMechanism **out)', '{',
          f'    ChemStatus status={symbol}_create(out);',
          '    if (!status) status=chem_mechanism_reorder(*out);', '    return status;', '}')
    ops = []
    for i, r in enumerate(plan.reactions):
        ops.append(['double t=s->conditions[c].temperature, pressure=s->conditions[c].pressure, density=s->conditions[c].air_density;',
                    '(void)t; (void)pressure; (void)density;'] + rate_statements(r, i))
    e.blocks('rates', 'ChemState *s', 's', ['    size_t n=s->cells;'], ops)
    e.put(f'void {symbol}_update_rates(ChemState *s) {{ rates(s); }}')

    forcing, jacobian = [], []
    for i, r in enumerate(plan.reactions):
        op = [f'double rate=k[{i}*n+c];'] + [f'rate *= y[{v}*n+c];' for v in r['lhs']]
        op += [f'f[{v}*n+c] -= rate;' for v in r['lhs']]
        op += [f'f[{v}*n+c] += {real(y)}*rate;' for v, y in r['rhs']]
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
                 'n,k,y,f', ['    (void)y;'], ops)
        e.put(f'static void {name}(const ChemState *s,const double *y,double *f)', '{',
              f'    memset(f,0,{count}*s->cells*sizeof(double));', f'    {name}_add(s->cells,s->rates,y,f);', '}')

    index, n = plan.index, plan.n
    factor = []
    for pivot in range(n):
        diag = index[pivot,pivot]
        # Error checks cannot appear in an omp simd structured block. A scalar
        # precheck is done in factor() below, each pivot immediately before use.
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
        e.put(f'static ChemStatus factor_{part}(size_t n,double *restrict v)', '{')
        for pivot in range(start,min(n,start+chunk)):
            diag = index[pivot,pivot]
            e.put('    { int bad=0;', '    #pragma omp simd reduction(|:bad)',
                  f'    for (size_t c=0;c<n;++c) bad |= v[{diag}*n+c]==0 || !isfinite(v[{diag}*n+c]);',
                  '    if (bad) return CHEM_SINGULAR;', '    }')
            if fused:
                e.loop(factor[pivot]) if len(factor[pivot]) > 1 else None
            else:
                for statement in factor[pivot][1:]:
                    e.loop([factor[pivot][0], statement] if 'inverse' in statement else [statement])
        e.put('    return CHEM_OK;', '}')
    e.put('static ChemStatus factor(const ChemSparse *a,size_t n,double *v,void *ctx)', '{', '    (void)a; (void)ctx; (void)n; (void)v;')
    for part in range((n+chunk-1)//chunk):
        e.put(f'    {{ ChemStatus status=factor_{part}(n,v); if (status) return status; }}')
    e.put('    return CHEM_OK;', '}')
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
    e.blocks('triangular', 'size_t n,const double *restrict v,double *restrict x', 'n,v,x', [], ops)
    e.put('static void solve(const ChemSparse *a,size_t n,const double *v,double *x,void *ctx)',
          '{ (void)a; (void)ctx; triangular(n,v,x); }')
    e.put(f'void {symbol}_specialized_attach(ChemSolver *s)', '{',
          '    s->linear=(ChemLinearSolver){factor,solve,NULL};', '}')
    e.put('static void diagonal(size_t n,double alpha,double *v)', '{', '    (void)n; (void)alpha; (void)v;')
    for row in range(n):
        e.loop([f'v[{index[row,row]}*n+c] += alpha;'])
    e.put('}')
    e.put(f'#define CHEM_SOLVE {symbol}_solve',
          '#define chem_forcing forcing',
          '#define chem_solver_negative_jacobian(s,state,y) jacobian(state,y,(s)->matrix)',
          '#define chem_add_mass_diagonal(m,pattern,cells,alpha,matrix) diagonal(cells,alpha,matrix)',
          '#include "rosenbrock_step.inc"')
    return '\n'.join(e.lines) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('configuration', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--symbol', default='configured')
    parser.add_argument('--chunk', type=int, default=32, help='operations/pivots per emitted function (positive)')
    parser.add_argument('--fusion', type=int, default=1, help='operations per SIMD cell loop (positive)')
    parser.add_argument('--split-lu', action='store_true', help='emit separate cell loops for each LU update')
    args = parser.parse_args()
    if args.chunk < 1 or args.fusion < 1:
        parser.error('--chunk and --fusion must be positive')
    result = generate(read_mechanism(args.configuration), args.symbol, args.chunk, not args.split_lu, args.fusion)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(result)


if __name__ == '__main__':
    main()
