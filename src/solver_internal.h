#ifndef CHEM_SOLVER_INTERNAL_H
#define CHEM_SOLVER_INTERNAL_H

#include "kinetics.h"

int chem_make_pattern(const chem_mechanism *m, chem_sparse_pattern *p);
void chem_free_pattern(chem_sparse_pattern *p);

#endif
