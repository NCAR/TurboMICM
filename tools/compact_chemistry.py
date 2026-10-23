"""Share generated mass-action kernels between reactions of the same shape.

Records remain in the exact source/derivative order. A shape fixes the number
of multiplies and signed updates, so the cell loop has no variable inner loops.
"""
from configure import number


def emit_kernel(e, name, jobs, output, count, fallback, negative_jacobian):
    # Each job is (rate index, ordered input columns, plain output columns,
    #               ordered (product column, yield) pairs).
    shapes = sorted({(len(inputs), len(plain), len(products))
                     for _, inputs, plain, products in jobs})
    shape_index = {shape: i for i, shape in enumerate(shapes)}
    max_inputs = max([len(inputs) for _, inputs, _, _ in jobs] + [1])
    max_outputs = max([len(plain) + len(products) for _, _, plain, products in jobs] + [1])
    max_products = max([len(products) for _, _, _, products in jobs] + [1])
    if jobs:
        e.add(f'static const struct {{ size_t rate, shape, inputs[{max_inputs}], outputs[{max_outputs}];',
              f'    double yields[{max_products}]; }} {name}_jobs[] = {{')
        for rate, inputs, plain, products in jobs:
            shape = shape_index[len(inputs), len(plain), len(products)]
            ints = lambda values: '{' + ','.join(map(str, values or [0])) + '}'
            yields = '{' + ','.join(number(v) for _, v in products or [(0, 0)]) + '}'
            e.add(f'    {{{rate},{shape},{ints(inputs)},{ints(plain + [i for i, _ in products])},{yields}}},')
        e.add('};')
    e.add(f'static void {name}(const ChemState *s, const double *restrict y, double *restrict {output})', '{')
    e.state_width(fallback)
    e.add('    const double *restrict rates = s->rates;', '    (void)y; (void)rates;',
          f'    memset({output}, 0, {count} * nc * sizeof(double));')
    if jobs:
        e.begin_tiles()
        e.add(f'    for (size_t j = 0; j < {len(jobs)}; ++j) {{',
              f'        switch ({name}_jobs[j].shape) {{')
        for shape, ordinal in shape_index.items():
            inputs, plain, products = shape
            # Resolve addresses outside the cell loop. C has no reference type;
            # fields retain their static table spelling instead of a copied job.
            job = f'{name}_jobs[j]'
            e.add(f'        case {ordinal}: {{',
                  f'            const double *restrict rate = rates + {job}.rate * nc;')
            for i in range(inputs):
                e.add(f'            const double *restrict input_{i} = y + {job}.inputs[{i}] * nc;')
            for i in range(plain + products):
                if output == 'a' and e.tile and e.tiled_matrix:
                    # Keep lane subtraction in the subscript so every pointer
                    # remains within its allocated matrix tile.
                    address = f'a + base * {e.nnz} + {job}.outputs[{i}] * (end - base)'
                else:
                    address = f'{output} + {job}.outputs[{i}] * nc'
                e.add(f'            double *output_{i} = {address};')
            for i in range(products):
                e.add(f'            const double yield_{i} = {job}.yields[{i}];')
            e.add('            #pragma omp simd', f'            for ({e.cell_range}) {{', '                double v = rate[c];')
            for i in range(inputs):
                e.add(f'                v *= input_{i}[c];')
            subscript = 'c - base' if output == 'a' and e.tile and e.tiled_matrix else 'c'
            for i in range(plain):
                e.add(f'                output_{i}[{subscript}] ' + ('+= v;' if negative_jacobian else '-= v;'))
            for i in range(products):
                e.add(f'                output_{plain + i}[{subscript}] ' + ('-=' if negative_jacobian else '+=') + f' yield_{i} * v;')
            e.add('            }', '            break;', '        }')
        e.add('        }', '    }')
        e.end_tiles()
    e.add('}')


def emit(e, mech, index):
    n = len(mech['species'])
    forcing = [(i, r['reactants'], r['reactants'], r['products']) for i, r in enumerate(mech['reactions'])]
    emit_kernel(e, 'sp_forcing', forcing, 'f', n, 'chem_forcing(s, y, f)', negative_jacobian=False)
    jacobian = []
    for col in range(n):
        for ri, r in enumerate(mech['reactions']):
            rr = r['reactants']
            if col not in rr:
                continue
            others = [i for j, i in enumerate(rr) if j != rr.index(col)]
            for _ in range(rr.count(col)):
                jacobian.append((ri, others, [index[i, col] for i in others + [col]],
                                 [(index[i, col], v) for i, v in r['products']]))
    emit_kernel(e, 'sp_jacobian', jacobian, 'a', len(index), 'chem_negative_jacobian(s, y, a)', negative_jacobian=True)
