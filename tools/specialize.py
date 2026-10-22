"""Compile a MICM configuration to fixed-index C11 chemistry and sparse solves.

Parsing/units live in configure.py. Plan builds the symbolic chemistry and LU;
Emitter turns ordered scalar operations into SIMD loops across independent cells.
No TS1 names, reaction counts, or measured input values occur in this compiler.
"""
import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path

from configure import generate, read_config


def literal(value):
    """Round-trip a finite binary64 value into a C floating-point constant."""
    value = float(value)
    if not math.isfinite(value):
        raise ValueError('non-finite generated constant')
    return repr(value)


def markowitz(pattern):
    """Stable diagonal ordering, with the generic library's tie breaking.

    Do not insert diagonal edges here: the library orders the chemistry pattern
    before adding the identity. A negative degree product sorts after finite
    costs, matching unsigned arithmetic without depending on size_t's width.
    """
    rows = [set(row) for row in pattern]
    alive = set(range(len(rows)))
    order = []
    while alive:
        def key(v):
            cost = (len(rows[v] & alive)-1) * (sum(v in rows[i] for i in alive)-1)
            return (cost if cost >= 0 else math.inf, v)
        pivot = min(alive, key=key)
        order.append(pivot)
        alive.remove(pivot)
        for i in alive:
            if pivot in rows[i]:
                rows[i].update(rows[pivot] & alive)
    return order


@dataclass
class Reaction:
    reactants: list
    products: list
    parameterized: int
    rate: dict


class Plan:
    """Symbolic model in exactly the generic constructor's state order."""
    def __init__(self, model):
        self.model = model
        active = [s['name'] for s in model.species if s.get('parameterization') != 'air density']
        if not active:
            raise ValueError('mechanism has no evolving species')
        old = {name: i for i, name in enumerate(active)}
        pattern = [set() for _ in active]
        for r in model.reactions:
            columns = [old[x] for x in r['reactants'] if x in old]
            rows = columns + [old[x] for x, _ in r['products'] if x in old]
            for i in rows:
                pattern[i].update(columns)
        self.order = markowitz(pattern)
        self.names = [active[i] for i in self.order]
        indices = {name: i for i, name in enumerate(self.names)}
        self.reactions = []
        for r in model.reactions:
            reactants = [indices[x] for x in r['reactants'] if x in indices]
            products = [(indices[x], y) for x, y in r['products'] if x in indices]
            self.reactions.append(Reaction(reactants, products, len(r['reactants'])-len(reactants), r))
        pattern = [set() for _ in active]
        for r in self.reactions:
            for i in r.reactants + [i for i, _ in r.products]:
                pattern[i].update(r.reactants)
        for i in range(len(active)):
            pattern[i].add(i)
        for pivot in range(len(active)):
            upper = {k for k in pattern[pivot] if k > pivot}
            for row in range(pivot+1, len(active)):
                if pivot in pattern[row]:
                    pattern[row].update(upper)
        self.rows = [sorted(row) for row in pattern]
        self.index = {(i, j): k for k, (i, j) in enumerate(
            (i, j) for i, row in enumerate(self.rows) for j in row)}

    def lu_operations(self, algorithm='mozart'):
        """Mozart elimination: inverse, lower scaling, outer products."""
        a = lambda i, j: f'a[{self.index[i,j]}*n+c]'
        if algorithm == 'doolittle':
            for i, row in enumerate(self.rows):
                for k in row:
                    if k >= i:
                        products = [f'value -= {a(i,j)}*{a(j,k)};' for j in range(i)
                                    if (i,j) in self.index and (j,k) in self.index]
                        if products:
                            yield [f'double value = {a(i,k)};'] + products + [f'{a(i,k)} = value;']
                yield [f'd[{i}*n+c] = 1.0/{a(i,i)};']
                for k in range(i+1,len(self.rows)):
                    if (k,i) in self.index:
                        products = [f'value -= {a(k,j)}*{a(j,i)};' for j in range(i)
                                    if (k,j) in self.index and (j,i) in self.index]
                        yield [f'double value = {a(k,i)};'] + products + [f'{a(k,i)} = value*d[{i}*n+c];']
            return
        for i, row in enumerate(self.rows):
            yield [f'd[{i}*n+c] = 1.0/{a(i,i)};']
            lower = [j for j in range(i+1, len(self.rows)) if (j, i) in self.index]
            for j in lower:
                yield [f'{a(j,i)} *= d[{i}*n+c];']
            for k in row:
                if k > i:
                    for j in lower:
                        yield [f'{a(j,k)} -= {a(j,i)}*{a(i,k)};']


class Emitter:
    """Bound compiler work by splitting large kernels into named C functions.

    Each unit is a reaction, derivative, scalar elimination, or triangular row.
    Fusing units across the cell loop preserves each cell's operation order.
    """
    def __init__(self, name, block, fixed, fusion):
        self.name, self.block, self.fixed, self.lines = name, block, fixed, []
        self.fusion = fusion
        self.kernels = []

    def kernel(self, kind, arguments, units, prelude=(), result='void', finish=()):
        units = list(units)
        calls, fixed_calls = [], []
        declarations = ', '.join(decl for decl, _ in arguments)
        values = ', '.join(value for _, value in arguments)
        for start in range(0, len(units), self.block):
            name = f'{self.name}_{kind}_{start//self.block}'
            for fixed in range(1+self.fixed):
                args = arguments[1:] if fixed else arguments
                decls = ', '.join(decl for decl, _ in args)
                if fixed:
                    decls = 'size_t count, '+decls
                suffix = '_fixed' if fixed else ''
                self.lines += [f'void {name}{suffix}({decls});']
                self.kernels += [f'void {name}{suffix}({decls}) {{']
                if fixed:
                    self.kernels += ['  const size_t n = CH_DEFAULT_GROUP_SIZE;']
                self.kernels += [f'  (void){value};' for _, value in args]
                limit = 'count' if fixed else 'n'
                chunk = units[start:start+self.block]
                for offset in range(0,len(chunk),self.fusion):
                    self.kernels += ['  #pragma omp simd', f'  for (size_t c = 0; c < {limit}; ++c) {{']
                    for unit in chunk[offset:offset+self.fusion]:
                        self.kernels += ['    {'] + ['      '+line for line in unit] + ['    }']
                    self.kernels += ['  }']
                self.kernels += ['}']
            calls.append(f'  {name}({values});')
            fixed_calls.append(f'  {name}_fixed(CH_KERNEL_WIDTH, '+', '.join(v+'+base' if '*' in d else v for d, v in arguments[1:])+');')
        self.lines += [f'{result} ch_{kind}_{self.name}({declarations}) {{']
        # Suppress unused parameters for valid mechanisms with no reactions.
        self.lines += [f'  (void){value};' for _, value in arguments]
        self.lines += ['  '+line for line in prelude]
        if self.fixed and units:
            self.lines += ['  if (n == CH_DEFAULT_GROUP_SIZE && n % CH_KERNEL_WIDTH == 0) {',
                           '    for (size_t base = 0; base < n; base += CH_KERNEL_WIDTH) {']
            self.lines += ['    '+s for s in fixed_calls] + ['    }', '  } else {']
            self.lines += ['  '+s for s in calls] + ['  }']
        else:
            self.lines += calls
        self.lines += ['  '+line for line in finish] + ['}']


def rate_expression(reaction, model):
    """The same evaluation order and rate formulas as src/rates.c, fixed p[]."""
    r = reaction.rate
    p = list(r['p'])
    kind = r['type']
    out = []
    if kind == 5:  # branched rate: prepare constants once, at generation time
        p.append(2e-22 * 6.02214076e23 * 1e-6 * math.exp(p[3]))
        a = p[5] * (2.45e19 / 6.02214076e23 * 1e6)
        b = .43 * (293.0/298.0)**-8.0
        v = a/(1+a/b) * .41**(1/(1+math.log10(a/b)**2))
        p.append(v*(1-p[2])/p[2])
    if kind == 6:  # surface: derive coefficients from the reactant metadata
        species = next(s for s in model.species if s['name'] == r['reactants'][0])
        p[0] = float(species['diffusion coefficient [m2 s-1]'])
        p[1] = 8.0*(1.380649e-23*6.02214076e23)/(math.pi*float(species['molecular weight [kg mol-1]']))
    p = list(map(literal, p))
    arr = lambda offset, ref: f'{p[offset]} * exp({p[offset+2]}/t) * pow(t/{ref}, {p[offset+1]})'
    if kind == 0:
        expression = arr(0, p[3])+f' * (1.0 + {p[4]}*pressure)'
    elif kind in (1, 2, 3):
        if kind == 3:
            out += [f'double k0 = {arr(0,p[3])};', f'double ki = {arr(4,p[7])};', f'double km = {arr(8,p[11])};']
            fc, exponent = p[12:14]
        else:
            out += [f'double k0 = {arr(0,"300.0")};', f'double ki = {arr(3,"300.0")};']
            fc, exponent = p[6:8]
        out += ['double ratio = k0 * air / ki;']
        expression = f'k0 * {"1.0" if kind == 2 else "air"} / (1.0 + ratio) * pow({fc}, {exponent}/({exponent}+pow(log10(ratio),2.0)))'
        if kind == 3:
            out += [f'double kf = {expression};']
            expression = 'kf + km * (1.0 - kf / ki)'
    elif kind == 4:
        expression = f'{p[0]} * exp(-({p[1]})/t + {p[2]}/(t*t*t))'
    elif kind == 5:
        out += [f'double a = {p[5]}*air, b = 0.43*pow(t/298.0,-8.0);',
                'double v = a/(1.0+a/b) * pow(0.41,1.0/(1.0+pow(log10(a/b),2.0)));',
                f'double pre = {p[0]}*exp(-({p[1]})/t);']
        expression = f'pre*({p[6]}/({p[6]}+v))' if r['p'][4] == 0 else f'pre*(v/(v+{p[6]}))'
    elif kind == 6:
        idx = r['custom_index']
        out += [f'double radius = custom[{idx}*n+c], number = custom[{idx+1}*n+c];']
        expression = f'4.0 * number * 3.14159265358979323846 * radius * radius / (radius/{p[0]} + 4.0/(sqrt({p[1]}*t)*{p[2]}))'
    elif kind == 7:
        expression = f'custom[{r["custom_index"]}*n+c] * {p[0]}'
    elif kind == 8:
        out += ['double poly = 0.0, power = 1.0;']
        for coef in p[5:]:
            out += [f'poly += {coef} * power;', 'power *= t;']
        expression = 'poly * '+arr(0,p[3])+f' * (1.0 + {p[4]}*pressure)'
    elif kind == 9:
        expression = f'{p[0]}*exp({p[1]}/t)*{p[2]}'
    else:
        raise ValueError(f'no emitter for rate type {kind}')
    out += [f'double value = {expression};']
    out += ['value *= air;'] * reaction.parameterized
    return out


def compile_model(model, name, block=32, fixed=1, lu='mozart', fusion=1):
    if block < 1 or fusion < 1:
        raise ValueError('kernel block and fusion must be positive')
    plan = Plan(model)
    emit = Emitter(name, block, fixed, fusion)
    nvars, nnz = len(plan.names), len(plan.index)
    preamble = ['/* Generated by tools/specialize.py. */', '#include "internal.h"',
                  '#ifndef CH_DEFAULT_GROUP_SIZE', '#define CH_DEFAULT_GROUP_SIZE 32', '#endif',
                  '#ifndef CH_KERNEL_WIDTH', '#define CH_KERNEL_WIDTH 4', '#endif',
                  '_Static_assert(CH_DEFAULT_GROUP_SIZE > 0 && CH_KERNEL_WIDTH > 0, "positive cell widths required");']
    emit.lines = preamble + [generate(model,name)]
    emit.kernels = preamble[:]
    emit.kernel('rates_kernel', [('size_t n','n'), ('const ch_conditions *restrict cond','cond'),
                ('const double *restrict custom','custom'), ('double *restrict rates','rates')],
        ([ 'double t = cond[c].temperature, air = cond[c].air_density, pressure = cond[c].pressure;',
           '(void)t; (void)air; (void)pressure; (void)custom;'] + rate_expression(r, model) +
         [f'rates[{i}*n+c] = value;'] for i,r in enumerate(plan.reactions)))
    emit.lines += [f'int ch_rates_{name}(ch_state *s) {{', '  if (!s) return CH_INVALID;',
                   f'  ch_rates_kernel_{name}(s->cells,s->conditions,s->custom,s->rates);', '  return CH_OK;', '}']
    args = [('size_t n','n'), ('const double *restrict rates','rates'),
            ('const double *restrict y','y'), ('double *restrict out','out')]
    units = []
    for ri,r in enumerate(plan.reactions):
        lines = [f'double rate = rates[{ri}*n+c];']
        lines += [f'rate *= y[{j}*n+c];' for j in r.reactants]
        lines += [f'out[{j}*n+c] -= rate;' for j in r.reactants]
        lines += [f'out[{j}*n+c] += {literal(v)}*rate;' for j,v in r.products]
        units.append(lines)
    emit.kernel('force_kernel',args,units,[f'memset(out,0,{nvars}*n*sizeof(double));'])
    units = []
    for col in range(nvars):
        for ri,r in enumerate(plan.reactions):
            for which,reactant in enumerate(r.reactants):
                if reactant != col:
                    continue
                # Retain multiplicity and ordered additions for repeated terms.
                first = r.reactants.index(col)
                others = r.reactants[:first]+r.reactants[first+1:]
                lines = [f'double rate = rates[{ri}*n+c];']
                lines += [f'rate *= y[{j}*n+c];' for j in others]
                lines += [f'out[{plan.index[j,col]}*n+c] += rate;' for j in others+[col]]
                lines += [f'out[{plan.index[j,col]}*n+c] -= {literal(v)}*rate;' for j,v in r.products]
                units.append(lines)
    emit.kernel('jacobian_kernel',args,units,[f'memset(out,0,{nnz}*n*sizeof(double));'])
    for kind in ('force','jacobian'):
        emit.lines += [f'void ch_{kind}_{name}(const ch_state *s, const double *y, double *out) {{',
                       f'  ch_{kind}_kernel_{name}(s->cells,s->rates,y,out);', '}']
    units = [[f'a[{plan.index[i,i]}*n+c] += alpha;'] for i in range(nvars)]
    emit.kernel('shift_kernel',[('size_t n','n'),('double *restrict a','a'),('double alpha','alpha')],units)
    emit.kernel('factor_kernel',[('size_t n','n'),('double *restrict a','a'),('double *restrict d','d')],plan.lu_operations(lu))
    emit.lines += [f'void ch_factor_{name}(ch_state *s, double alpha) {{',
                   f'  ch_shift_kernel_{name}(s->cells,s->jacobian,alpha);',
                   f'  ch_factor_kernel_{name}(s->cells,s->jacobian,s->inverse_diagonal);', '}']
    units = []
    for i,row in enumerate(plan.rows):
        lower = [j for j in row if j < i]
        if lower:
            units.append([f'double value = x[{i}*n+c];'] +
                         [f'value -= a[{plan.index[i,j]}*n+c]*x[{j}*n+c];' for j in lower] +
                         [f'x[{i}*n+c] = value;'])
    for i in reversed(range(nvars)):
        units.append([f'double value = x[{i}*n+c];'] +
                     [f'value -= a[{plan.index[i,j]}*n+c]*x[{j}*n+c];' for j in plan.rows[i] if j > i] +
                     [f'x[{i}*n+c] = value*d[{i}*n+c];'])
    emit.kernel('linear_kernel',[('size_t n','n'),('const double *restrict a','a'),
                ('const double *restrict d','d'),('double *restrict x','x')],units)
    emit.lines += [f'void ch_linear_{name}(ch_state *s, double *x) {{',
                   f'  ch_linear_kernel_{name}(s->cells,s->jacobian,s->inverse_diagonal,x);', '}',
                   f'#define CH_SPECIES {nvars}', f'#define CH_SOLVE ch_solve_{name}',
                   f'#define CH_FORCE ch_force_{name}', f'#define CH_JACOBIAN ch_jacobian_{name}',
                   f'#define CH_FACTOR ch_factor_{name}', f'#define CH_LINEAR ch_linear_{name}',
                   '#include "specialized_solver.h"']
    header = '\n'.join(['/* Generated by tools/specialize.py. */','#include "chem.h"',
        f'ch_model *ch_build_{name}(void);', f'int ch_init_{name}(ch_state *);',
        f'int ch_rates_{name}(ch_state *);',
        f'void ch_force_{name}(const ch_state *, const double *, double *);',
        f'void ch_jacobian_{name}(const ch_state *, const double *, double *);',
        f'void ch_factor_{name}(ch_state *, double);', f'void ch_linear_{name}(ch_state *, double *);',
        f'int ch_solve_{name}(ch_state *, double, const ch_rosenbrock_parameters *, ch_stats *);',''])
    metadata = dict(species=nvars, reactions=len(plan.reactions), nonzeros=nnz,
                    lu_operations=sum(1 for _ in plan.lu_operations(lu)), kernel_block=block, fixed_width=bool(fixed), lu=lu, fusion=fusion,
                    species_order=plan.names)
    return '\n'.join(emit.lines)+'\n', '\n'.join(emit.kernels)+'\n', header, metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--name',default='mechanism')
    parser.add_argument('--block',type=int,default=32)
    parser.add_argument('--fixed',type=int,choices=(0,1),default=1)
    parser.add_argument('--lu',choices=('mozart','doolittle'),default='mozart')
    parser.add_argument('--fusion',type=int,default=1)
    args = parser.parse_args()
    code, kernels, header, metadata = compile_model(read_config(args.config),args.name,args.block,args.fixed,args.lu,args.fusion)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(code)
    args.output.with_name(args.output.stem+'_kernels.c').write_text(kernels)
    args.output.with_suffix('.h').write_text(header)
    args.output.with_suffix('.json').write_text(json.dumps(metadata,indent=2)+'\n')


if __name__ == '__main__':
    main()
