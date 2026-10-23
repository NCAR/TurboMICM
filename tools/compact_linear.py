"""Emit a compact interpreter for a mechanism's fixed sparse operation plan.

This alternative shares short vector loops between structurally identical
updates. It trades immediate indices for a smaller instruction footprint.
"""


def emit(e, mech, index, reciprocal, batch=4, wide=True):
    n = len(mech['species'])
    batch = min(batch, max([sum((j, i) in index for j in range(i + 1, n)) for i in range(n)] + [1]))
    batch = max(batch, 1)
    rows = [sorted(j for i, j in index if i == row) for row in range(n)]
    diagonal = [index[i, i] for i in range(n)]
    row_start = [sum(len(row) for row in rows[:i]) for i in range(n + 1)]
    lower, lower_start, operations, operation_start = [], [0], [], [0]
    for i in range(n):
        lower_rows = [j for j in range(i + 1, n) if (j, i) in index]
        lower += [index[j, i] for j in lower_rows]
        lower_start.append(len(lower))
        for k in range(i + 1, n):
            if (i, k) not in index:
                continue
            # Adjacent lower updates share this upper entry. Batch them in one
            # cell loop while preserving the original column-then-row order.
            for begin in range(0, len(lower_rows), batch):
                group = lower_rows[begin:begin + batch]
                operations.append((index[i, k], [index[j, i] for j in group],
                                   [index[j, k] for j in group]))
        operation_start.append(len(operations))

    largest = max(n, len(index), len(lower), len(operations), batch)
    integer_type = 'size_t'
    if not wide:
        for bits in (8, 16, 32):
            if largest < 2 ** bits:
                integer_type = f'uint{bits}_t'
                break

    def array(name, values):
        e.add(f'static const {integer_type} sp_{name}[] = {{' + ','.join(map(str, values or [0])) + '};')

    array('diagonal', diagonal)
    array('row', row_start)
    array('column', [j for _, j in sorted(index)])
    array('lower', lower)
    array('lower_start', lower_start)
    array('operation_start', operation_start)
    e.add(f'static const struct {{ {integer_type} right, count, left[{batch}], destination[{batch}]; }} sp_operations[] = {{')
    for right, left, destination in operations or [(0, [], [])]:
        ints = lambda values: '{' + ','.join(map(str, values or [0])) + '}'
        e.add(f'    {{{right},{len(left)},{ints(left)},{ints(destination)}}},')
    e.add('};')

    def begin(signature, fallback, const=False, returns=False):
        e.add(signature, '{', '    (void)m;')
        e.matrix_width(fallback, returns=returns)
        e.begin_tiles()
        if not e.tile:
            e.add('    const size_t base = 0, end = nc;')
        e.add('    const size_t lanes = end - base;')
        stride = 'lanes' if e.tile and e.tiled_matrix else 'nc'
        offset = f'base * {len(index)}' if e.tile and e.tiled_matrix else 'base'
        e.add(f'    const size_t stride = {stride};',
              f'    {"const " if const else ""}double *restrict matrix = a + {offset};')

    def end(returns=False):
        e.end_tiles()
        if returns:
            e.add('    return CHEM_OK;')
        e.add('}')

    begin('static void sp_alpha(const ChemMechanism *m, size_t cells, double *restrict a, double alpha)',
          'chem_alpha_minus_jacobian(m, cells, a, alpha)')
    e.add(f'    for (size_t i = 0; i < {n}; ++i) {{',
          '        double *restrict pivot = matrix + sp_diagonal[i] * stride;',
          '        #pragma omp simd',
          '        for (size_t c = 0; c < lanes; ++c) pivot[c] += alpha;', '    }')
    end()
    begin('static ChemStatus sp_factor(const ChemMechanism *m, size_t cells, double *restrict a)',
          'chem_factor(m, cells, a)', returns=True)
    e.add(f'    for (size_t i = 0; i < {n}; ++i) {{',
          '        double *restrict pivot = matrix + sp_diagonal[i] * stride;',
          '        for (size_t c = 0; c < lanes; ++c)',
          '            if (pivot[c] == 0 || !isfinite(pivot[c])) return CHEM_SINGULAR;')
    if reciprocal:
        e.add('        #pragma omp simd', '        for (size_t c = 0; c < lanes; ++c) pivot[c] = 1.0 / pivot[c];')
    e.add('        for (size_t j = sp_lower_start[i]; j < sp_lower_start[i + 1]; ++j) {',
          '            double *restrict dst = matrix + sp_lower[j] * stride;',
          '            #pragma omp simd',
          '            for (size_t c = 0; c < lanes; ++c) dst[c] *= ' + ('pivot[c];' if reciprocal else '1.0 / pivot[c];'),
          '        }',
          '        for (size_t j = sp_operation_start[i]; j < sp_operation_start[i + 1]; ++j) {',
          '            const double *restrict right = matrix + sp_operations[j].right * stride;',
          '            (void)right;',
          '            switch (sp_operations[j].count) {')
    for count in sorted({len(left) for _, left, _ in operations}):
        e.add(f'            case {count}: {{')
        for k in range(count):
            e.add(f'                double *restrict dst_{k} = matrix + sp_operations[j].destination[{k}] * stride;',
                  f'                const double *restrict left_{k} = matrix + sp_operations[j].left[{k}] * stride;')
        e.add('                #pragma omp simd', '                for (size_t c = 0; c < lanes; ++c) {',
              '                    double v = right[c];')
        for k in range(count):
            e.add(f'                    dst_{k}[c] -= left_{k}[c] * v;')
        e.add('                }', '                break;', '            }')
    e.add('            }', '        }', '    }')
    end(returns=True)
    begin('static void sp_backsolve(const ChemMechanism *m, size_t cells, const double *restrict a, double *restrict x)',
          'chem_linear_solve(m, cells, a, x)', const=True)
    e.add(f'    for (size_t i = 0; i < {n}; ++i) {{',
          '        double *restrict xi = x + i * nc + base;',
          '        for (size_t j = sp_row[i]; j < sp_diagonal[i]; ++j) {',
          '            const double *restrict value = matrix + j * stride;',
          '            const double *restrict xj = x + sp_column[j] * nc + base;',
          '            #pragma omp simd',
          '            for (size_t c = 0; c < lanes; ++c) xi[c] -= value[c] * xj[c];',
          '        }', '    }',
          f'    for (size_t ii = {n}; ii > 0; --ii) {{',
          '        size_t i = ii - 1;', '        double *restrict xi = x + i * nc + base;',
          '        for (size_t j = sp_diagonal[i] + 1; j < sp_row[i + 1]; ++j) {',
          '            const double *restrict value = matrix + j * stride;',
          '            const double *restrict xj = x + sp_column[j] * nc + base;',
          '            #pragma omp simd',
          '            for (size_t c = 0; c < lanes; ++c) xi[c] -= value[c] * xj[c];',
          '        }', '        const double *restrict pivot = matrix + sp_diagonal[i] * stride;',
          '        #pragma omp simd',
          '        for (size_t c = 0; c < lanes; ++c) xi[c] ' + ('*= pivot[c];' if reciprocal else '/= pivot[c];'),
          '    }')
    end()
