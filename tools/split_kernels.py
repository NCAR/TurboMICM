"""Bound compiler work by splitting emitted loop sequences into C functions.

Only complete SIMD loops are split. Arguments are explicit and arithmetic order
is unchanged. This is an emitter pass over our own structured line sequences,
not a parser for arbitrary C.
"""


def split_kernels(lines, limit, indirect=False, omit=()):
    signatures = {
        'forcing': ('const chem_model *model, size_t cells, const double *restrict rates, const double *restrict y, double *restrict out', 'model, cells, rates, y, out'),
        'jacobian': ('const chem_model *model, size_t cells, const double *restrict rates, const double *restrict y, double *restrict out', 'model, cells, rates, y, out'),
        'add_alpha': ('size_t cells, double *restrict a, double alpha', 'cells, a, alpha'),
        'factor': ('size_t cells, double *restrict a, double *restrict scratch', 'cells, a, scratch'),
        'linear': ('size_t cells, const double *restrict a, double *restrict x', 'cells, a, x'),
    }
    output = []
    i = 0
    while i < len(lines):
        name = next((n for n in signatures if lines[i].startswith('static void '+n+'(')), None)
        if name is None:
            output.append(lines[i]); i += 1
            continue
        end = lines.index('}', i)
        if name in omit:
            i = end+1
            continue
        body = lines[lines.index('{', i)+1:end]
        declarations = []
        loops = []
        j = 0
        while j < len(body):
            if body[j] == '    #pragma omp simd':
                last = body.index('    }', j)
                loops.append(body[j:last+1]); j = last+1
            else:
                declarations.append(body[j]); j += 1
        groups = []
        current = []
        for loop in loops:
            if current and len(current)+len(loop)>limit:
                groups.append(current); current = []
            current += loop
        if current:
            groups.append(current)
        arguments, call = signatures[name]
        if len(groups) <= 1:
            output += lines[i:end+1]
        else:
            for index, group in enumerate(groups):
                output += [f'static void {name}_part{index}({arguments})', '{'] + declarations + group + ['}']
            if indirect:
                for index in range(len(groups)):
                    output += [f'static void (*volatile {name}_call{index})({arguments}) = {name}_part{index};']
            output += [f'static void {name}({arguments})', '{']
            target = 'call' if indirect else 'part'
            output += [f'    {name}_{target}{index}({call});' for index in range(len(groups))]
            output += ['}']
        i = end+1
    return output
