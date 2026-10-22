"""Emit immutable C mechanism data and optional generic specialized kernels."""
import argparse
import json
from pathlib import Path
from configure import read, number, rate_initializer


def emit(config, symbol, kernels=True):
    species=config['species']; reactions=config['reactions']
    variables=[i for i,s in enumerate(species) if s.get('parameterization')!='air density']
    var={s:i for i,s in enumerate(variables)}
    n=len(var)
    # Reference Markowitz selection, including isolated-node unsigned costs.
    graph=set()
    for rx in reactions:
        rs=[var[x] for x in rx['reactants'] if x in var]
        ps=[var[x] for x,v in rx['products'] if x in var]
        graph.update((row,col) for col in rs for row in rs+ps)
    alive=set(range(n)); permutation=[]
    # Only the ordering of large unsigned costs matters; all supported hosts have
    # size_t at least 32 bits, far larger than mechanism degrees.
    mask=(1<<64)-1
    while alive:
        def cost(v):
            row=sum((v,j) in graph for j in alive)
            col=sum((i,v) in graph for i in alive)
            return ((row-1)*(col-1)) & mask
        pivot=min(alive,key=lambda v:(cost(v),v))
        permutation.append(pivot); alive.remove(pivot)
        ins=[i for i in alive if (i,pivot) in graph]
        cols=[j for j in alive if (pivot,j) in graph]
        graph.update((i,j) for i in ins for j in cols)
    variables=[variables[i] for i in permutation]
    var={s:i for i,s in enumerate(variables)}
    lines=['/* Generated from a mechanism configuration; do not edit. */', '#include "chem.h"', '']
    lines.append('static const chem_species species[]={')
    for sp in species:
        param=sp.get('parameterization')
        if param not in (None,'air density'): raise ValueError('Unknown parameterization: '+param)
        lines.append('    {.name='+json.dumps(sp['name'])+', .type=CHEM_SPECIES_'+('THIRD_BODY' if param else 'VARIABLE')+
                     ', .absolute_tolerance='+number(sp.get('absolute tolerance',1e-3))+
                     ', .initial='+number(sp.get('initial concentration [mol m-3]',0))+
                     ', .molecular_weight='+number(sp.get('molecular weight [kg mol-1]',0))+
                     ', .diffusion='+number(sp.get('diffusion coefficient [m2 s-1]',0))+'},')
    lines.append('};')
    for i,rx in enumerate(reactions):
        lines.append(f'static const size_t r{i}[]={{'+','.join(list(map(str,rx['reactants'])) or ['0'])+'};')
        prod=','.join('{'+str(s)+','+number(y)+'}' for s,y in rx['products']) or '{0,0}'
        lines.append(f'static const chem_product p{i}[]={{'+prod+'};')
    lines.append('static const chem_reaction reactions[]={')
    if not reactions: lines.append('    {0}')
    for i,rx in enumerate(reactions):
        lines.append(f'    {{{len(rx["reactants"])},{len(rx["products"])},r{i},p{i},'+rate_initializer(rx,species)+'},')
    lines.append('};')
    defaults=config.get('initial rate parameters',{})
    lines.append('static const chem_named_value defaults[]={')
    for key,value in defaults.items(): lines.append('    {'+json.dumps(key)+','+number(value)+'},')
    if not defaults: lines.append('    {0,0}')
    lines.append('};')
    if kernels:
        # Identical symbolic fill to the runtime planner, independent of LU algorithm.
        pattern={(i,i) for i in range(n)}
        jobs=[]
        for i,rx in enumerate(reactions):
            rs=[var[x] for x in rx['reactants'] if x in var]
            ps=[(var[x],v) for x,v in rx['products'] if x in var]
            for col in rs:
                for row in rs+[p[0] for p in ps]: pattern.add((row,col))
                deps=rs.copy(); deps.remove(col)
                jobs.append((col,i,deps,ps))
        for i in range(n):
            for k in range(i+1,n):
                if (i,k) in pattern:
                    for j in range(i+1,n):
                        if (j,i) in pattern: pattern.add((j,k))
        flat={pos:i for i,pos in enumerate(sorted(pattern))}
        def start(name):
            lines.extend([f'static void {name}(const chem_model *model, size_t cells, const double *restrict rates,',
                          '                      const double *restrict y, double *restrict out)', '{', '    (void)model;'])
        def loop(statements):
            lines.extend(['    #pragma omp simd', '    for (size_t c=0; c<cells; ++c) {'])
            lines.extend('        '+s for s in statements)
            lines.append('    }')
        start('forcing')
        for i,rx in enumerate(reactions):
            rs=[var[x] for x in rx['reactants'] if x in var]
            ps=[(var[x],v) for x,v in rx['products'] if x in var]
            statements=[f'double r=rates[{i}*cells+c];']
            statements += [f'r*=y[{s}*cells+c];' for s in rs]
            statements += [f'out[{s}*cells+c]-=r;' for s in rs]
            statements += [f'out[{s}*cells+c]+={number(v)}*r;' for s,v in ps]
            loop(statements)
        lines.append('}')
        start('jacobian')
        for col,i,rs,ps in sorted(jobs,key=lambda j:(j[0],j[1])):
            statements=[f'double r=rates[{i}*cells+c];']
            statements += [f'r*=y[{s}*cells+c];' for s in rs]
            statements += [f'out[{flat[s,col]}*cells+c]+=r;' for s in rs+[col]]
            statements += [f'out[{flat[s,col]}*cells+c]-={number(v)}*r;' for s,v in ps]
            loop(statements)
        lines.append('}')
    conditions=config.get('initial conditions',{})
    T=float(conditions.get('temperature',0)); P=float(conditions.get('pressure',0))
    density=P/(1.380649e-23*6.02214076e23*T) if T else 0
    lines.extend([f'const chem_mechanism {symbol}={{',
                  f'    .n_species={len(species)}, .n_reactions={len(reactions)}, .species=species, .reactions=reactions,',
                  f'    .n_defaults={len(defaults)}, .defaults=defaults,',
                  '    .conditions={'+','.join(map(number,[T,P,density]))+'},',
                  '    .forcing=forcing, .jacobian=jacobian' if kernels else '    .forcing=0, .jacobian=0', '};',''])
    return '\n'.join(lines)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('configuration'); p.add_argument('output'); p.add_argument('--symbol',default='mechanism')
    p.add_argument('--no-kernels',action='store_true')
    args=p.parse_args()
    if not args.symbol.isidentifier() or not args.symbol.isascii(): p.error('Invalid C symbol')
    Path(args.output).parent.mkdir(parents=True,exist_ok=True)
    Path(args.output).write_text(emit(read(args.configuration),args.symbol,not args.no_kernels))
