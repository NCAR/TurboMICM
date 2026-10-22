/* Diagnostic only; all timed chemistry remains C, using the public interfaces. */
#include "chem.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
ChemModel *chem_load_ts1(int);
ChemModel *chem_load_ts1_specialized(int);
void chem_factor_ts1_specialized(const ChemModel *, size_t, size_t, double *);
void chem_linear_ts1_specialized(const ChemModel *, size_t, size_t, const double *, double *);
void chem_rates_ts1_specialized(ChemState *);
int main(int argc, char **argv) {
  size_t cells = argc > 1 ? (size_t)strtoul(argv[1],NULL,10) : 32;
  size_t repeats = argc > 2 ? (size_t)strtoul(argv[2],NULL,10) : 10000;
  for (int specialized = 0; specialized < 2; ++specialized) {
    ChemModel *m = specialized ? chem_load_ts1_specialized(1) : chem_load_ts1(1);
    ChemState *s = chem_state_create(m,cells);
    if (!m || !s) return 1;
    chem_calculate_rates(s); chem_jacobian(s,s->variables,s->jacobian); chem_alpha_minus_jacobian(s,1000.0);
    size_t nj = m->nnz*cells, ny = m->nvariables*cells;
    double *original = malloc(nj*sizeof(double));
    if (!original) return 1;
    memcpy(original,s->jacobian,nj*sizeof(double));
    ChemFactor factor = specialized ? chem_factor_ts1_specialized : chem_lu_mozart;
    ChemSolve solve = specialized ? chem_linear_ts1_specialized : chem_linear_solve;
    double start = omp_get_wtime();
    for (size_t i = 0; i < repeats; ++i) chem_forcing(s,s->variables,s->forcing);
    double forcing = omp_get_wtime()-start; start = omp_get_wtime();
    for (size_t i = 0; i < repeats; ++i) chem_jacobian(s,s->variables,s->jacobian);
    double jacobian = omp_get_wtime()-start; start = omp_get_wtime();
    for (size_t i = 0; i < repeats; ++i) { memcpy(s->jacobian,original,nj*sizeof(double)); factor(m,cells,cells,s->jacobian); }
    double lu = omp_get_wtime()-start; start = omp_get_wtime();
    for (size_t i = 0; i < repeats; ++i) { memcpy(s->stages[0],s->variables,ny*sizeof(double)); solve(m,cells,cells,s->jacobian,s->stages[0]); }
    double triangular = omp_get_wtime()-start; start = omp_get_wtime();
    for (size_t i = 0; i < repeats; ++i) { if (specialized) chem_rates_ts1_specialized(s); else chem_calculate_rates(s); }
    double rates = omp_get_wtime()-start;
    printf("mode=%s cells=%zu repeats=%zu forcing_ms=%.3f jacobian_ms=%.3f factor_ms=%.3f solve_ms=%.3f rates_ms=%.3f\n",
           specialized ? "specialized" : "tier1",cells,repeats,forcing*1000,jacobian*1000,lu*1000,triangular*1000,rates*1000);
    free(original); chem_state_free(s); chem_model_free(m);
  }
  return 0;
}
