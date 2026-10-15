#ifndef CHEM_INTERNAL_H
#define CHEM_INTERNAL_H
#include "chem.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
char *chem_strdup(const char *s);
int chem_resize(void **p, size_t count, size_t size);
void chem_axpy(size_t n, double *restrict x, double factor, const double *restrict y);
#endif
