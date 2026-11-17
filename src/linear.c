#include "internal.h"

#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static void *allocate(size_t n, size_t size)
{
    if (size && n > SIZE_MAX / size) return NULL;
    return calloc(n ? n : 1, size);
}

int chem_linear_add_alpha(ChemLinearSystem *s, double alpha)
{
    if (!s || !s->values || !s->lookup || !s->algebraic || !s->cells || !s->species_count ||
        s->species_count>SIZE_MAX/s->species_count || s->nonzeros>SIZE_MAX/s->cells || !isfinite(alpha)) return CHEM_INVALID;
    size_t n=s->species_count;
    for (size_t i=0;i<n;++i) if (!s->algebraic[i]) {
        int position=s->lookup[i*n+i];
        if (position<0 || (size_t)position>=s->nonzeros) return CHEM_INVALID;
        double *diag=s->values+(size_t)position*s->cells;
        #pragma omp simd
        for (size_t c=0;c<s->cells;++c) diag[c]+=alpha;
    }
    return CHEM_OK;
}

/* Symbolic elimination is independent of the concentrations. It is kept
 * separate from numeric factorization so another linear solver can replace it. */
int chem_make_symbolic(ChemState *s)
{
    const ChemMechanism *m = s->mechanism;
    size_t n = m->species_count;
    if (n > INT_MAX || n > SIZE_MAX / n || n*n > INT_MAX) return CHEM_INVALID;
    unsigned char *pattern = allocate(n*n, sizeof(*pattern));
    if (!pattern) return CHEM_MEMORY;
    for (size_t i = 0; i < n; ++i) pattern[i*n+i] = 1;
    for (size_t r = 0; r < m->reaction_count; ++r) {
        const ChemReaction *reaction = &m->reactions[r];
        if ((reaction->reactant_count && !m->reactants) || (reaction->product_count && !m->products)) {
            free(pattern); return CHEM_INVALID;
        }
        /* Source-only reactions still need their product indices validated. */
        for (size_t j = 0; j < reaction->product_count; ++j) {
            int row = m->products[reaction->product_start+j].species;
            if (row < 0 || (size_t)row >= n) { free(pattern); return CHEM_INVALID; }
        }
        for (size_t k = 0; k < reaction->reactant_count; ++k) {
            int col = m->reactants[reaction->reactant_start+k];
            if (col < 0 || (size_t)col >= n) { free(pattern); return CHEM_INVALID; }
            for (size_t j = 0; j < reaction->reactant_count; ++j) {
                int row = m->reactants[reaction->reactant_start+j];
                if (row < 0 || (size_t)row >= n) { free(pattern); return CHEM_INVALID; }
                if (!s->algebraic[row]) pattern[(size_t)row*n+(size_t)col] = 1;
            }
            for (size_t j = 0; j < reaction->product_count; ++j) {
                int row = m->products[reaction->product_start+j].species;
                if (row < 0 || (size_t)row >= n) { free(pattern); return CHEM_INVALID; }
                if (!s->algebraic[row]) pattern[(size_t)row*n+(size_t)col] = 1;
            }
        }
    }
    for (size_t k=0; k<s->constraint_count; ++k) {
        const ChemConstraint *constraint = &s->constraints[k];
        size_t row = (size_t)constraint->algebraic_species;
        for (size_t i=0; i<constraint->reactant_count; ++i)
            pattern[row*n+(size_t)constraint->reactants[i].species] = 1;
        for (size_t i=0; i<constraint->product_count; ++i)
            pattern[row*n+(size_t)constraint->products[i].species] = 1;
    }
    for (size_t i = 0; i < n; ++i)
        for (size_t k = i+1; k < n; ++k)
            if (pattern[i*n+k])
                for (size_t j = i+1; j < n; ++j)
                    if (pattern[j*n+i]) pattern[j*n+k] = 1;

    s->index = allocate(n*n, sizeof(*s->index));
    s->row_start = allocate(n+1, sizeof(*s->row_start));
    s->pivots = allocate(n, sizeof(*s->pivots));
    if (!s->index || !s->row_start || !s->pivots) { free(pattern); return CHEM_MEMORY; }
    for (size_t i = 0; i < n*n; ++i) s->index[i] = -1;
    for (size_t i = 0; i < n; ++i) {
        s->row_start[i] = s->nonzeros;
        for (size_t j = 0; j < n; ++j)
            if (pattern[i*n+j]) s->index[i*n+j] = (int)s->nonzeros++;
    }
    s->row_start[n] = s->nonzeros;
    s->columns = allocate(s->nonzeros, sizeof(*s->columns));
    if (!s->columns) { free(pattern); return CHEM_MEMORY; }
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            if (pattern[i*n+j]) s->columns[s->index[i*n+j]] = (int)j;

    size_t nl = 0, nu = 0, nz = 0;
    for (size_t i = 0; i < n; ++i) {
        size_t count_l = 0, count_u = 0;
        for (size_t j = i+1; j < n; ++j) count_l += pattern[j*n+i];
        for (size_t k = i+1; k < n; ++k) count_u += pattern[i*n+k];
        if (count_l && count_u > (SIZE_MAX - nz) / count_l) { free(pattern); return CHEM_MEMORY; }
        nz += count_l * count_u; nl += count_l; nu += count_u;
    }
    s->lower = allocate(nl, sizeof(*s->lower));
    s->upper = allocate(nu, sizeof(*s->upper));
    s->updates = allocate(nz, sizeof(*s->updates));
    if (!s->lower || !s->upper || !s->updates) { free(pattern); return CHEM_MEMORY; }
    nl = nu = nz = 0;
    for (size_t i = 0; i < n; ++i) {
        ChemPivot *pivot = &s->pivots[i];
        pivot->diagonal = s->index[i*n+i];
        pivot->first_lower = nl;
        for (size_t j = i+1; j < n; ++j)
            if (pattern[j*n+i]) { s->lower[nl++] = s->index[j*n+i]; ++pivot->lowers; }
        pivot->first_upper = nu;
        for (size_t k = i+1; k < n; ++k) if (pattern[i*n+k]) {
            ChemUpper *u = &s->upper[nu++];
            u->position = s->index[i*n+k];
            u->first_update = nz;
            for (size_t j = i+1; j < n; ++j) if (pattern[j*n+i]) {
                s->updates[nz++] = (ChemUpdate){s->index[j*n+k], s->index[j*n+i]};
                ++u->updates;
            }
            ++pivot->uppers;
        }
    }
    free(pattern);
    return CHEM_OK;
}

int chem_factor(ChemState *s, double alpha)
{
    size_t n = s->mechanism->species_count, cells = s->cells;
    double *a = s->matrix, *inverse = s->scratch;
    ChemLinearSystem view={n,cells,s->nonzeros,s->row_start,s->columns,s->index,s->algebraic,s->matrix};
    int shifted=chem_linear_add_alpha(&view,alpha);
    if (shifted) return shifted;
    for (size_t i = 0; i < n; ++i) {
        ChemPivot *p = &s->pivots[i];
        double *diag = a + (size_t)p->diagonal*cells;
        for (size_t c = 0; c < cells; ++c)
            if (diag[c] == 0.0 || !isfinite(diag[c])) return CHEM_SINGULAR;
        #pragma omp simd
        for (size_t c = 0; c < cells; ++c) inverse[c] = 1.0 / diag[c];
        for (size_t j = 0; j < p->lowers; ++j) {
            double *lower = a + (size_t)s->lower[p->first_lower+j]*cells;
            #pragma omp simd
            for (size_t c = 0; c < cells; ++c) lower[c] *= inverse[c];
        }
        for (size_t k = 0; k < p->uppers; ++k) {
            ChemUpper *u = &s->upper[p->first_upper+k];
            const double *upper = a + (size_t)u->position*cells;
            for (size_t j = 0; j < u->updates; ++j) {
                ChemUpdate *op = &s->updates[u->first_update+j];
                double *target = a + (size_t)op->target*cells;
                const double *lower = a + (size_t)op->lower*cells;
                #pragma omp simd
                for (size_t c = 0; c < cells; ++c) target[c] -= lower[c] * upper[c];
            }
        }
    }
    return CHEM_OK;
}

static ChemLinearSystem linear_view(const ChemState *s)
{
    return (ChemLinearSystem){s->mechanism->species_count,s->cells,s->nonzeros,
                              s->row_start,s->columns,s->index,s->algebraic,s->matrix};
}

static void substitute_view(const ChemLinearSystem *system, double *x)
{
    size_t n = system->species_count, cells = system->cells;
    const double *a = system->values;
    for (size_t i = 0; i < n; ++i) {
        double *xi = x + i*cells;
        for (size_t k = system->row_start[i]; k < system->row_start[i+1] && (size_t)system->columns[k] < i; ++k) {
            const double *lij = a + k*cells, *xj = x + (size_t)system->columns[k]*cells;
            #pragma omp simd
            for (size_t c = 0; c < cells; ++c) xi[c] -= lij[c]*xj[c];
        }
    }
    for (size_t i = n; i-- > 0;) {
        double *xi = x + i*cells;
        size_t diagonal = (size_t)system->lookup[i*n+i];
        for (size_t k = diagonal+1; k < system->row_start[i+1]; ++k) {
            const double *uij = a + k*cells, *xj = x + (size_t)system->columns[k]*cells;
            #pragma omp simd
            for (size_t c = 0; c < cells; ++c) xi[c] -= uij[c]*xj[c];
        }
        const double *diag = a + diagonal*cells;
        #pragma omp simd
        for (size_t c = 0; c < cells; ++c) xi[c] /= diag[c];
    }
}

void chem_substitute(const ChemState *s, double *x)
{
    ChemLinearSystem view=linear_view(s);
    substitute_view(&view,x);
}

int chem_linear_factor(ChemState *s, double alpha)
{
    if (!s->linear_solver.factor) return chem_factor(s,alpha);
    ChemLinearSystem view=linear_view(s);
    return s->linear_solver.factor(&view,alpha,s->linear_solver_context);
}

int chem_linear_solve(ChemState *s, double *values)
{
    if (!s->linear_solver.solve) { chem_substitute(s,values); return CHEM_OK; }
    ChemLinearSystem view=linear_view(s);
    return s->linear_solver.solve(&view,values,s->linear_solver_context);
}

static int doolittle_factor(ChemLinearSystem *s, double alpha, void *context)
{
    (void)context;
    size_t n=s->species_count, cells=s->cells;
    double *a=s->values;
    for (size_t i=0;i<n;++i) if (!s->algebraic[i]) {
        double *diag=a+(size_t)s->lookup[i*n+i]*cells;
        #pragma omp simd
        for (size_t cell=0;cell<cells;++cell) diag[cell]+=alpha;
    }
    for (size_t i=0;i<n;++i) {
        for (size_t k=i;k<n;++k) {
            int target=s->lookup[i*n+k];
            if (target<0) continue;
            double *aik=a+(size_t)target*cells;
            for (size_t j=0;j<i;++j) {
                int left=s->lookup[i*n+j], right=s->lookup[j*n+k];
                if (left<0 || right<0) continue;
                const double *aij=a+(size_t)left*cells, *ajk=a+(size_t)right*cells;
                #pragma omp simd
                for (size_t cell=0;cell<cells;++cell) aik[cell]-=aij[cell]*ajk[cell];
            }
        }
        const double *diag=a+(size_t)s->lookup[i*n+i]*cells;
        for (size_t cell=0;cell<cells;++cell)
            if (diag[cell]==0 || !isfinite(diag[cell])) return CHEM_SINGULAR;
        for (size_t k=i+1;k<n;++k) {
            int target=s->lookup[k*n+i];
            if (target<0) continue;
            double *aki=a+(size_t)target*cells;
            for (size_t j=0;j<i;++j) {
                int left=s->lookup[k*n+j], right=s->lookup[j*n+i];
                if (left<0 || right<0) continue;
                const double *akj=a+(size_t)left*cells, *aji=a+(size_t)right*cells;
                #pragma omp simd
                for (size_t cell=0;cell<cells;++cell) aki[cell]-=akj[cell]*aji[cell];
            }
            #pragma omp simd
            for (size_t cell=0;cell<cells;++cell) aki[cell]/=diag[cell];
        }
    }
    return CHEM_OK;
}

static int doolittle_solve(const ChemLinearSystem *s, double *rhs, void *context)
{ (void)context; substitute_view(s,rhs); return CHEM_OK; }

ChemLinearSolver chem_doolittle_linear_solver(void)
{ return (ChemLinearSolver){doolittle_factor,doolittle_solve}; }

/* The general view implementation exposes Mozart factorization to callers
 * with their own sparse matrices. Solver states use the precomputed traversal
 * above, which performs the same elimination without searching the pattern. */
static int mozart_factor(ChemLinearSystem *s, double alpha, void *context)
{
    (void)context;
    size_t n=s->species_count, cells=s->cells;
    double *a=s->values;
    for (size_t i=0;i<n;++i) if (!s->algebraic[i]) {
        int at=s->lookup[i*n+i]; if (at<0) return CHEM_INVALID;
        double *diag=a+(size_t)at*cells;
        #pragma omp simd
        for (size_t cell=0;cell<cells;++cell) diag[cell]+=alpha;
    }
    for (size_t i=0;i<n;++i) {
        int diagonal=s->lookup[i*n+i]; if (diagonal<0) return CHEM_INVALID;
        const double *diag=a+(size_t)diagonal*cells;
        for (size_t cell=0;cell<cells;++cell)
            if (diag[cell]==0 || !isfinite(diag[cell])) return CHEM_SINGULAR;
        for (size_t j=i+1;j<n;++j) {
            int lower=s->lookup[j*n+i]; if (lower<0) continue;
            double *aji=a+(size_t)lower*cells;
            #pragma omp simd
            for (size_t cell=0;cell<cells;++cell) aji[cell]*=1.0/diag[cell];
        }
        for (size_t k=i+1;k<n;++k) {
            int upper=s->lookup[i*n+k]; if (upper<0) continue;
            const double *aik=a+(size_t)upper*cells;
            for (size_t j=i+1;j<n;++j) {
                int lower=s->lookup[j*n+i]; if (lower<0) continue;
                int target=s->lookup[j*n+k]; if (target<0) return CHEM_INVALID;
                const double *aji=a+(size_t)lower*cells;
                double *ajk=a+(size_t)target*cells;
                #pragma omp simd
                for (size_t cell=0;cell<cells;++cell) ajk[cell]-=aji[cell]*aik[cell];
            }
        }
    }
    return CHEM_OK;
}

ChemLinearSolver chem_mozart_linear_solver(void)
{ return (ChemLinearSolver){mozart_factor,doolittle_solve}; }

int chem_linear_factor_copy(const ChemLinearSolver *solver, const ChemLinearSystem *input,
                            ChemLinearSystem *factors, double alpha, void *context)
{
    if (!solver || !solver->factor || !input || !factors || !input->values || !factors->values ||
        !input->cells || !input->row_start || !factors->row_start || !input->columns || !factors->columns ||
        !factors->lookup || !factors->algebraic ||
        input->values==factors->values || input->species_count!=factors->species_count ||
        input->cells!=factors->cells || input->nonzeros!=factors->nonzeros ||
        input->nonzeros > SIZE_MAX / input->cells ||
        input->nonzeros*input->cells > SIZE_MAX / sizeof(double)) return CHEM_INVALID;
    for (size_t i=0;i<=input->species_count;++i)
        if (input->row_start[i]!=factors->row_start[i]) return CHEM_INVALID;
    for (size_t i=0;i<input->nonzeros;++i)
        if (input->columns[i]!=factors->columns[i]) return CHEM_INVALID;
    memcpy(factors->values,input->values,input->nonzeros*input->cells*sizeof(double));
    return solver->factor(factors,alpha,context);
}
