/* Errors, allocation helpers, small vectors, hashing, version. */
#include "internal.h"

#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>

void *tmc_malloc(size_t bytes)
{
  void *p = malloc(bytes ? bytes : 1);
  if (!p)
  {
    fprintf(stderr, "tmc: out of memory (%zu bytes)\n", bytes);
    abort();
  }
  return p;
}

void *tmc_calloc(size_t count, size_t bytes)
{
  void *p = calloc(count ? count : 1, bytes ? bytes : 1);
  if (!p)
  {
    fprintf(stderr, "tmc: out of memory (%zu x %zu bytes)\n", count, bytes);
    abort();
  }
  return p;
}

void *tmc_realloc(void *old, size_t bytes)
{
  void *p = realloc(old, bytes ? bytes : 1);
  if (!p)
  {
    fprintf(stderr, "tmc: out of memory (%zu bytes)\n", bytes);
    abort();
  }
  return p;
}

/* The solver's double arrays. With TMC_ALIGN (a power of two, at least 16; a
 * Derecho tuning, default 64) they come from C11's aligned_alloc at that
 * alignment, so a block's W doubles never straddle a cache line and vector
 * loads are aligned; TMC_ALIGN 0 keeps malloc. */
#ifndef TMC_ALIGN
#define TMC_ALIGN 64
#endif
double *tmc_alloc_doubles(size_t count)
{
#if TMC_ALIGN > 0
  size_t bytes = count * sizeof(double);
  bytes = (bytes + (size_t)TMC_ALIGN - 1) / (size_t)TMC_ALIGN * (size_t)TMC_ALIGN;
  void *p = aligned_alloc((size_t)TMC_ALIGN, bytes ? bytes : (size_t)TMC_ALIGN);
  if (!p)
  {
    fprintf(stderr, "tmc: out of memory (%zu bytes)\n", bytes);
    abort();
  }
  return p;
#else
  return tmc_malloc(count * sizeof(double));
#endif
}

void tmc_free_doubles(double *p)
{
  free(p);
}

char *tmc_strdup(const char *s)
{
  size_t n = strlen(s) + 1;
  char *p = tmc_malloc(n);
  memcpy(p, s, n);
  return p;
}

void tmc_strlcpy(char *dst, const char *src, size_t size)
{
  if (size == 0)
    return;
  size_t n = strlen(src);
  if (n >= size)
    n = size - 1;
  memcpy(dst, src, n);
  dst[n] = '\0';
}

void tmc_ivec_push(tmc_ivec *a, int x)
{
  if (a->n == a->cap)
  {
    a->cap = a->cap ? 2 * a->cap : 16;
    a->v = tmc_realloc(a->v, (size_t)a->cap * sizeof(int));
  }
  a->v[a->n++] = x;
}

void tmc_dvec_push(tmc_dvec *a, double x)
{
  if (a->n == a->cap)
  {
    a->cap = a->cap ? 2 * a->cap : 16;
    a->v = tmc_realloc(a->v, (size_t)a->cap * sizeof(double));
  }
  a->v[a->n++] = x;
}

void tmc_ivec_free(tmc_ivec *a)
{
  free(a->v);
  a->v = NULL;
  a->n = a->cap = 0;
}

void tmc_dvec_free(tmc_dvec *a)
{
  free(a->v);
  a->v = NULL;
  a->n = a->cap = 0;
}

unsigned long long tmc_fnv1a(unsigned long long h, const void *bytes, size_t n)
{
  const unsigned char *p = bytes;
  for (size_t i = 0; i < n; ++i)
  {
    h ^= p[i];
    h *= 1099511628211ULL;
  }
  return h;
}

tmc_error_category tmc_error_set(tmc_error *err, tmc_error_category category, int code, const char *format, ...)
{
  if (err)
  {
    err->category = category;
    err->code = code;
    va_list args;
    va_start(args, format);
    vsnprintf(err->message, sizeof err->message, format, args);
    va_end(args);
  }
  return category;
}

void tmc_error_clear(tmc_error *err)
{
  if (err)
  {
    err->category = TMC_OK;
    err->code = 0;
    err->message[0] = '\0';
  }
}

const char *tmc_error_category_name(tmc_error_category category)
{
  switch (category)
  {
    case TMC_OK: return "OK";
    case TMC_ERROR_CONFIGURATION: return "Configuration";
    case TMC_ERROR_PROCESS: return "Process";
    case TMC_ERROR_RATE_CONSTANT: return "Rate Constant";
    case TMC_ERROR_SPECIES: return "Species";
    case TMC_ERROR_MATRIX: return "Matrix";
    case TMC_ERROR_STATE: return "State";
    case TMC_ERROR_SOLVER: return "Solver";
    case TMC_ERROR_INTERNAL: return "Internal Error";
    case TMC_ERROR_CONSTRAINT: return "Constraint";
  }
  return "Unknown";
}

const char *tmc_solver_state_name(tmc_solver_state s)
{
  switch (s)
  {
    case TMC_NOT_YET_CALLED: return "Not Yet Called";
    case TMC_RUNNING: return "Running";
    case TMC_CONVERGED: return "Converged";
    case TMC_CONVERGENCE_EXCEEDED_MAX_STEPS: return "Convergence Exceeded Max Steps";
    case TMC_STEP_SIZE_TOO_SMALL: return "Step Size Too Small";
    case TMC_REPEATEDLY_SINGULAR_MATRIX: return "Repeatedly Singular Matrix";
    case TMC_NAN_DETECTED: return "NaNDetected";
    case TMC_INF_DETECTED: return "InfDetected";
    case TMC_ACCEPTING_UNCONVERGED_INTEGRATION: return "AcceptingUnconvergedIntegration";
    case TMC_CONSTRAINT_INITIALIZATION_FAILED: return "Constraint Initialization Failed";
  }
  return "Unknown";
}

/* The version tracks the MICM commit the port reproduces (v3.7.0-259). */
const char *tmc_version(void)
{
  return "3.7.0";
}
unsigned tmc_version_major(void)
{
  return 3;
}
unsigned tmc_version_minor(void)
{
  return 7;
}
unsigned tmc_version_patch(void)
{
  return 0;
}
