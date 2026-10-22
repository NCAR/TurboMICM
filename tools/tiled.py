"""Lower kernels and the solver to cell tiles without changing adaptive groups.

Public state stays species-major. Private workspace packs each tile's species,
reactions, and sparse elements together. Padding lanes have zero rates and are
excluded from the original species-major error reduction.
"""
import re
from split_kernels import split_kernels


def fixed_kernels(lines, width, function_lines, indirect=True):
    lines = split_kernels(lines, function_lines, indirect=indirect)
    output = []
    names = ('forcing', 'jacobian', 'add_alpha', 'factor', 'linear')
    for line in lines:
        line = re.sub(r'\b('+'|'.join(names)+r')(?=\b|_part|_call)', r'tile_\1', line)
        if line.startswith('static void ') or line.startswith('    const double *restrict rates,'):
            line = line.replace('size_t cells,', 'size_t unused_cells,')
        output.append(line)
        if line == '{':
            output.append(f'    const size_t cells = {width}; (void)unused_cells;')
    return output


def wrappers(species, reactions, nnz, width):
    lines = []
    for name, output_count in [('forcing', species), ('jacobian', nnz)]:
        lines += [f'''static void packed_{name}(const chem_model *m, size_t cells,
    const double *rates, const double *y, double *out)
{{
    for (size_t c=0; c<cells; c+={width})
        tile_{name}(m, {width}, rates+c*{reactions}, y+c*{species}, out+c*{output_count});
}}''']
    lines += [f'''static void packed_add_alpha(size_t cells, double *a, double alpha)
{{
    for (size_t c=0; c<cells; c+={width}) tile_add_alpha({width}, a+c*{nnz}, alpha);
}}
static void packed_factor(size_t cells, double *a, double *scratch)
{{
    for (size_t c=0; c<cells; c+={width}) tile_factor({width}, a+c*{nnz}, scratch);
}}
static void packed_linear(size_t cells, const double *a, double *x)
{{
    for (size_t c=0; c<cells; c+={width}) tile_linear({width}, a+c*{nnz}, x+c*{species});
}}''']
    return lines


def transform_solver(template, species, reactions, nnz, width):
    # Space for y, the six stage/work vectors, rates and the factorized Jacobian.
    count = 7*species + reactions + nnz
    setup = f'''
    chem_state *original = s;
    const size_t active = s->cells;
    if (active > SIZE_MAX-({width}-1)) {{ result.status=CHEM_SIZE; return result; }}
    const size_t padded = (active+{width}-1)/{width}*{width};
    if (padded > SIZE_MAX/sizeof(double)/{count}) {{ result.status=CHEM_SIZE; return result; }}
    const size_t required = padded*{count};
    if (s->specialized_capacity < required) {{
        double *grown = realloc(s->specialized_work, required*sizeof(double));
        if (!grown) {{ result.status=CHEM_MEMORY; return result; }}
        s->specialized_work = grown;
        s->specialized_capacity = required;
    }}
    chem_state view = *s;
    view.cells = padded;
    view.y = s->specialized_work;
    view.work = view.y + {species}*padded;
    view.rates = view.work + {6*species}*padded;
    view.jacobian = view.rates + {reactions}*padded;
    for (size_t i=0; i<{species}; ++i)
        for (size_t c=0; c<padded; ++c)
            view.y[(c/{width})*{species*width}+i*{width}+c%{width}] = c<active ? s->y[i*active+c] : 0.0;
    for (size_t i=0; i<{reactions}; ++i)
        for (size_t c=0; c<padded; ++c)
            view.rates[(c/{width})*{reactions*width}+i*{width}+c%{width}] = c<active ? s->rates[i*active+c] : 0.0;
    s = &view;
'''
    template = template.replace('    const size_t cells = s->cells,', setup+'    const size_t cells = s->cells,')
    for name in ('forcing', 'jacobian', 'add_alpha', 'factor', 'linear'):
        template = template.replace(name+'(', 'packed_'+name+'(')
    template = template.replace('c<cells; ++c) {\n                    size_t k = i*cells+c;',
                                f'c<active; ++c) {{\n                    size_t k = (c/{width})*{species*width}+i*{width}+c%{width};')
    template = template.replace('sqrt(sum/(double)n)', f'sqrt(sum/(double)({species}*active))')
    store = f'''    for (size_t i=0; i<{species}; ++i)
        for (size_t c=0; c<active; ++c)
            original->y[i*active+c] = y[(c/{width})*{species*width}+i*{width}+c%{width}];
'''
    template = template.replace('    result.stats.final_time = time;', store+'    result.stats.final_time = time;')
    return template
