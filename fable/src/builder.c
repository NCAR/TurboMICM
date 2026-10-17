/* The solver builder: from a system and reactions to the compiled plans a
 * solver and its states use (MICM's SolverBuilder::Build). */
#include "internal.h"

#include <stdint.h>
#include <stdio.h>

struct tmc_builder
{
  tmc_method method;
  tmc_rosenbrock_parameters ros;
  tmc_backward_euler_parameters be;
  bool has_system;
  tmc_system system;
  int n_reactions;
  tmc_reaction *reactions;
  bool reorder_state;
  bool ignore_unused_species;
  int width;
  tmc_lu_algorithm lu_algorithm;
  const tmc_kernels *kernels;
  tmc_step_controller controller;
  void *controller_context;
};

static tmc_builder *create(void)
{
  tmc_builder *b = tmc_calloc(1, sizeof *b);
  b->reorder_state = true;
  b->width = TMC_WIDTH;
  b->lu_algorithm = TMC_LU_MOZART;
  b->controller = tmc_default_step_controller;
  return b;
}

tmc_builder *tmc_builder_create_rosenbrock(tmc_rosenbrock_parameters parameters)
{
  tmc_builder *b = create();
  b->method = TMC_METHOD_ROSENBROCK;
  b->ros = parameters;
  b->be = tmc_backward_euler_defaults();
  return b;
}

tmc_builder *tmc_builder_create_backward_euler(tmc_backward_euler_parameters parameters)
{
  tmc_builder *b = create();
  b->method = TMC_METHOD_BACKWARD_EULER;
  b->be = parameters;
  b->ros = tmc_rosenbrock_three_stage();
  return b;
}

void tmc_builder_free(tmc_builder *b)
{
  if (!b)
    return;
  if (b->has_system)
    tmc_system_free(&b->system);
  tmc_reactions_free(b->reactions, b->n_reactions);
  free(b);
}

void tmc_builder_set_system(tmc_builder *b, const tmc_system *system)
{
  if (b->has_system)
    tmc_system_free(&b->system);
  b->system = tmc_system_copy(system);
  b->has_system = true;
}

void tmc_builder_set_reactions(tmc_builder *b, const tmc_reaction *reactions, int n)
{
  tmc_reactions_free(b->reactions, b->n_reactions);
  b->reactions = tmc_malloc((size_t)(n > 0 ? n : 1) * sizeof *b->reactions);
  for (int i = 0; i < n; ++i)
    b->reactions[i] = tmc_reaction_copy(&reactions[i]);
  b->n_reactions = n;
}

void tmc_builder_set_reorder_state(tmc_builder *b, bool reorder)
{
  b->reorder_state = reorder;
}

void tmc_builder_set_ignore_unused_species(tmc_builder *b, bool ignore)
{
  b->ignore_unused_species = ignore;
}

void tmc_builder_set_width(tmc_builder *b, int L)
{
  b->width = L < 1 ? 1 : L;
}

void tmc_builder_set_lu_algorithm(tmc_builder *b, tmc_lu_algorithm algorithm)
{
  b->lu_algorithm = algorithm;
}

void tmc_builder_set_kernels(tmc_builder *b, const tmc_kernels *kernels)
{
  b->kernels = kernels;
}

void tmc_builder_set_step_controller(tmc_builder *b, tmc_step_controller controller, void *context)
{
  b->controller = controller ? controller : tmc_default_step_controller;
  b->controller_context = context;
}

/* ---- helpers ---- */

static int name_index(char (*names)[TMC_NAME_MAX], int n, const char *name)
{
  for (int i = 0; i < n; ++i)
    if (strcmp(names[i], name) == 0)
      return i;
  return -1;
}

static int compare_type(const void *a, const void *b)
{
  const tmc_reaction *ra = a, *rb = b;
  return (int)ra->rate_constant.type - (int)rb->rate_constant.type;
}

/* A stable sort by rate constant type (qsort is not stable). */
static void sort_reactions(tmc_reaction *reactions, int n)
{
  tmc_reaction *sorted = tmc_malloc((size_t)(n > 0 ? n : 1) * sizeof *sorted);
  int k = 0;
  for (int type = 0; type < TMC_RATE_TYPE_COUNT; ++type)
    for (int i = 0; i < n; ++i)
      if ((int)reactions[i].rate_constant.type == type)
        sorted[k++] = reactions[i];
  memcpy(reactions, sorted, (size_t)n * sizeof *sorted);
  free(sorted);
  (void)compare_type;
}

static void hash_bytes(unsigned long long *h, const void *p, size_t n)
{
  *h = tmc_fnv1a(*h, p, n);
}

static void hash_text(unsigned long long *h, const char *s)
{
  *h = tmc_fnv1a(*h, s, strlen(s) + 1);
}

static void hash_int(unsigned long long *h, int v)
{
  char buf[32];
  snprintf(buf, sizeof buf, "%d", v);
  hash_text(h, buf);
}

static void hash_double(unsigned long long *h, double v)
{
  /* the IEEE bit pattern, little endian, matches tools/mechanism.py */
  unsigned char bytes[8];
  memcpy(bytes, &v, 8);
  unsigned x = 1;
  if (*(unsigned char *)&x != 1)
    for (int i = 0; i < 4; ++i)
    {
      unsigned char t = bytes[i];
      bytes[i] = bytes[7 - i];
      bytes[7 - i] = t;
    }
  hash_bytes(h, bytes, 8);
}

static void free_solver_plans(tmc_solver *s)
{
  if (s->processes)
    for (int i = 0; i < s->n_reactions; ++i)
    {
      free(s->processes[i].reactants);
      free(s->processes[i].products);
      free(s->processes[i].yields);
    }
  free(s->processes);
  free(s->job_process);
  free(s->job_independent);
  free(s->job_n_dependent);
  free(s->job_n_products);
  free(s->job_dependent_ids);
  free(s->job_product_ids);
  free(s->job_yields);
  free(s->job_flat_ids);
  tmc_sparse_free(s->jacobian_pattern);
  tmc_lu_free(s->lu);
  free(s->variable_names);
  free(s->custom_labels);
  free(s->absolute_tolerance);
}

void tmc_solver_free(tmc_solver *s)
{
  if (!s)
    return;
  free_solver_plans(s);
  tmc_system_free(&s->system);
  tmc_reactions_free(s->reactions, s->n_reactions);
  free(s);
}

tmc_solver *tmc_builder_build(const tmc_builder *b, tmc_error *err)
{
  tmc_error_clear(err);
  if (!b->has_system)
  {
    tmc_error_set(err, TMC_ERROR_SOLVER, TMC_SOLVER_MISSING_CHEMICAL_SYSTEM, "Missing chemical system.");
    return NULL;
  }
  int n_vars = tmc_system_state_size(&b->system);
  if (n_vars == 0)
  {
    tmc_error_set(err, TMC_ERROR_SOLVER, TMC_SOLVER_MISSING_CHEMICAL_SPECIES, "Provided chemical system contains no species.");
    return NULL;
  }
  for (int i = 0; i < b->n_reactions; ++i)
    if (tmc_reaction_validate(&b->reactions[i], err) != TMC_OK)
      return NULL;

  tmc_solver *s = tmc_calloc(1, sizeof *s);
  s->method = b->method;
  s->ros = b->ros;
  s->be = b->be;
  s->controller = b->controller;
  s->controller_context = b->controller_context;
  s->width = b->width;
  s->reorder_state = b->reorder_state;
  s->system = tmc_system_copy(&b->system);
  s->n_reactions = b->n_reactions;
  s->reactions = tmc_malloc((size_t)(b->n_reactions > 0 ? b->n_reactions : 1) * sizeof *s->reactions);
  for (int i = 0; i < b->n_reactions; ++i)
    s->reactions[i] = tmc_reaction_copy(&b->reactions[i]);
  sort_reactions(s->reactions, s->n_reactions);

  /* Species map: declaration order, then the diagonal Markowitz reordering of
   * the Jacobian pattern built in that order. */
  s->n_vars = n_vars;
  s->variable_names = tmc_malloc((size_t)n_vars * sizeof *s->variable_names);
  char (*declared)[TMC_NAME_MAX] = tmc_malloc((size_t)n_vars * sizeof *declared);
  tmc_system_unique_names(&s->system, declared, n_vars);
  memcpy(s->variable_names, declared, (size_t)n_vars * sizeof *declared);

  /* Reactant and product ids in declared order (skipping parameterized species) */
  int n_rxn = s->n_reactions;
  tmc_ivec *rids = tmc_calloc((size_t)(n_rxn > 0 ? n_rxn : 1), sizeof *rids);
  tmc_ivec *pids = tmc_calloc((size_t)(n_rxn > 0 ? n_rxn : 1), sizeof *pids);
  for (int r = 0; r < n_rxn; ++r)
  {
    const tmc_reaction *rx = &s->reactions[r];
    for (int i = 0; i < rx->n_reactants; ++i)
    {
      if (rx->reactants[i].parameterized)
        continue;
      int id = name_index(declared, n_vars, rx->reactants[i].name);
      if (id < 0)
      {
        tmc_error_set(err, TMC_ERROR_PROCESS, TMC_PROCESS_REACTANT_DOES_NOT_EXIST, "%s", rx->reactants[i].name);
        goto fail;
      }
      tmc_ivec_push(&rids[r], id);
    }
    for (int i = 0; i < rx->n_products; ++i)
    {
      if (rx->products[i].species.parameterized)
        continue;
      int id = name_index(declared, n_vars, rx->products[i].species.name);
      if (id < 0)
      {
        tmc_error_set(err, TMC_ERROR_PROCESS, TMC_PROCESS_PRODUCT_DOES_NOT_EXIST, "%s", rx->products[i].species.name);
        goto fail;
      }
      tmc_ivec_push(&pids[r], id);
    }
  }

  /* Unused species check */
  if (!b->ignore_unused_species)
  {
    unsigned char *used = tmc_calloc((size_t)n_vars, 1);
    for (int r = 0; r < n_rxn; ++r)
    {
      for (int i = 0; i < rids[r].n; ++i)
        used[rids[r].v[i]] = 1;
      for (int i = 0; i < pids[r].n; ++i)
        used[pids[r].v[i]] = 1;
    }
    char msg[1024] = "Unused species in chemical system:";
    bool any = false;
    for (int i = 0; i < n_vars; ++i)
      if (!used[i])
      {
        any = true;
        size_t len = strlen(msg);
        snprintf(msg + len, sizeof msg - len, " '%s'", declared[i]);
      }
    free(used);
    if (any)
    {
      size_t len = strlen(msg);
      snprintf(msg + len, sizeof msg - len, ".");
      tmc_error_set(err, TMC_ERROR_SOLVER, TMC_SOLVER_UNUSED_SPECIES, "%s", msg);
      goto fail;
    }
  }

  int *perm = tmc_malloc((size_t)n_vars * sizeof(int)); /* perm[new] = declared index */
  for (int i = 0; i < n_vars; ++i)
    perm[i] = i;
  if (b->reorder_state)
  {
    unsigned char *pattern = tmc_calloc((size_t)n_vars * (size_t)n_vars, 1);
    for (int r = 0; r < n_rxn; ++r)
      for (int a = 0; a < rids[r].n; ++a)
      {
        int ind = rids[r].v[a];
        for (int d = 0; d < rids[r].n; ++d)
          pattern[(size_t)rids[r].v[d] * n_vars + ind] = 1;
        for (int d = 0; d < pids[r].n; ++d)
          pattern[(size_t)pids[r].v[d] * n_vars + ind] = 1;
      }
    tmc_markowitz_reorder(n_vars, pattern, perm);
    free(pattern);
  }
  int *new_of = tmc_malloc((size_t)n_vars * sizeof(int));
  for (int k = 0; k < n_vars; ++k)
  {
    new_of[perm[k]] = k;
    memcpy(s->variable_names[k], declared[perm[k]], TMC_NAME_MAX);
  }
  free(perm);

  /* Process plans in the new order */
  s->processes = tmc_calloc((size_t)(n_rxn > 0 ? n_rxn : 1), sizeof *s->processes);
  s->custom_labels = tmc_malloc((size_t)(2 * n_rxn + 1) * sizeof *s->custom_labels);
  s->n_custom = 0;
  for (int r = 0; r < n_rxn; ++r)
  {
    tmc_process_plan *p = &s->processes[r];
    const tmc_reaction *rx = &s->reactions[r];
    p->n_reactants = rids[r].n;
    p->reactants = tmc_malloc((size_t)(p->n_reactants > 0 ? p->n_reactants : 1) * sizeof(int));
    for (int i = 0; i < p->n_reactants; ++i)
      p->reactants[i] = new_of[rids[r].v[i]];
    p->n_products = pids[r].n;
    p->products = tmc_malloc((size_t)(p->n_products > 0 ? p->n_products : 1) * sizeof(int));
    p->yields = tmc_malloc((size_t)(p->n_products > 0 ? p->n_products : 1) * sizeof(double));
    int k = 0;
    for (int i = 0; i < rx->n_products; ++i)
    {
      if (rx->products[i].species.parameterized)
        continue;
      p->products[k] = new_of[pids[r].v[k]];
      p->yields[k] = rx->products[i].coefficient;
      ++k;
    }
    p->n_parameterized = 0;
    for (int i = 0; i < rx->n_reactants; ++i)
      if (rx->reactants[i].parameterized && p->n_parameterized < 4)
        p->parameterized[p->n_parameterized++] = rx->reactants[i];
    p->custom_index = -1;
    char labels[2][TMC_NAME_MAX];
    int n_labels = tmc_rate_constant_custom_parameter_labels(&rx->rate_constant, labels);
    if (n_labels > 0)
      p->custom_index = s->n_custom;
    for (int i = 0; i < n_labels; ++i)
    {
      if (name_index(s->custom_labels, s->n_custom, labels[i]) >= 0)
      {
        tmc_error_set(err, TMC_ERROR_SOLVER, TMC_SOLVER_DUPLICATE_PARAMETER, "Duplicate parameter labels detected:\n  - %s (from reaction)\n", labels[i]);
        free(new_of);
        goto fail;
      }
      memcpy(s->custom_labels[s->n_custom++], labels[i], TMC_NAME_MAX);
    }
  }

  /* Jacobian jobs: (independent variable, process) pairs, stably ordered by variable */
  {
    tmc_ivec job_ind = { 0 }, job_proc = { 0 };
    for (int r = 0; r < n_rxn; ++r)
      for (int i = 0; i < s->processes[r].n_reactants; ++i)
      {
        tmc_ivec_push(&job_ind, s->processes[r].reactants[i]);
        tmc_ivec_push(&job_proc, r);
      }
    int n_jobs = job_ind.n;
    int *order = tmc_malloc((size_t)(n_jobs > 0 ? n_jobs : 1) * sizeof(int));
    int k = 0;
    for (int v = 0; v < n_vars; ++v)
      for (int j = 0; j < n_jobs; ++j)
        if (job_ind.v[j] == v)
          order[k++] = j;
    s->n_jobs = n_jobs;
    s->job_process = tmc_malloc((size_t)(n_jobs > 0 ? n_jobs : 1) * sizeof(int));
    s->job_independent = tmc_malloc((size_t)(n_jobs > 0 ? n_jobs : 1) * sizeof(int));
    s->job_n_dependent = tmc_malloc((size_t)(n_jobs > 0 ? n_jobs : 1) * sizeof(int));
    s->job_n_products = tmc_malloc((size_t)(n_jobs > 0 ? n_jobs : 1) * sizeof(int));
    tmc_ivec dep = { 0 }, prod = { 0 };
    tmc_dvec yields = { 0 };
    for (int q = 0; q < n_jobs; ++q)
    {
      int j = order[q];
      int r = job_proc.v[j], ind = job_ind.v[j];
      const tmc_process_plan *p = &s->processes[r];
      s->job_process[q] = r;
      s->job_independent[q] = ind;
      bool found = false;
      int nd = 0;
      for (int i = 0; i < p->n_reactants; ++i)
      {
        if (p->reactants[i] == ind && !found)
        {
          found = true;
          continue;
        }
        tmc_ivec_push(&dep, p->reactants[i]);
        ++nd;
      }
      s->job_n_dependent[q] = nd;
      for (int i = 0; i < p->n_products; ++i)
      {
        tmc_ivec_push(&prod, p->products[i]);
        tmc_dvec_push(&yields, p->yields[i]);
      }
      s->job_n_products[q] = p->n_products;
    }
    s->job_dependent_ids = dep.v;
    s->job_product_ids = prod.v;
    s->job_yields = yields.v;
    s->n_job_dependent_total = dep.n;
    s->n_job_products_total = prod.n;
    free(order);
    tmc_ivec_free(&job_ind);
    tmc_ivec_free(&job_proc);
  }

  /* Jacobian pattern (with the diagonal) and the LU */
  {
    tmc_sparse_builder *sb = tmc_sparse_builder_create(n_vars);
    for (int r = 0; r < n_rxn; ++r)
    {
      const tmc_process_plan *p = &s->processes[r];
      for (int a = 0; a < p->n_reactants; ++a)
      {
        int ind = p->reactants[a];
        for (int d = 0; d < p->n_reactants; ++d)
          tmc_sparse_builder_with_element(sb, p->reactants[d], ind, NULL);
        for (int d = 0; d < p->n_products; ++d)
          tmc_sparse_builder_with_element(sb, p->products[d], ind, NULL);
      }
    }
    for (int i = 0; i < n_vars; ++i)
      tmc_sparse_builder_with_element(sb, i, i, NULL);
    s->jacobian_pattern = tmc_sparse_build(sb, TMC_SPARSE_CSR, s->width, true);
    tmc_sparse_builder_free(sb);
    s->lu = tmc_lu_create(s->jacobian_pattern, b->lu_algorithm, err);
    if (!s->lu)
    {
      free(new_of);
      goto fail;
    }
    s->lu_pattern = s->lu->pattern;
  }

  /* Flat LU element ids for the Jacobian jobs: dependents, the independent's
   * diagonal, then the products. */
  {
    s->n_job_flat_total = s->n_job_dependent_total + s->n_jobs + s->n_job_products_total;
    s->job_flat_ids = tmc_malloc((size_t)(s->n_job_flat_total > 0 ? s->n_job_flat_total : 1) * sizeof(int));
    int d = 0, pr = 0, f = 0;
    for (int q = 0; q < s->n_jobs; ++q)
    {
      int ind = s->job_independent[q];
      for (int i = 0; i < s->job_n_dependent[q]; ++i)
        s->job_flat_ids[f++] = tmc_sparse_element(s->lu_pattern, s->job_dependent_ids[d++], ind);
      s->job_flat_ids[f++] = tmc_sparse_element(s->lu_pattern, ind, ind);
      for (int i = 0; i < s->job_n_products[q]; ++i)
        s->job_flat_ids[f++] = tmc_sparse_element(s->lu_pattern, s->job_product_ids[pr++], ind);
    }
  }

  /* Absolute tolerances: 1e-3 unless the species sets a property */
  s->absolute_tolerance = tmc_malloc((size_t)n_vars * sizeof(double));
  for (int k = 0; k < n_vars; ++k)
  {
    s->absolute_tolerance[k] = 1.0e-3;
    const tmc_phase_species *ps = tmc_phase_find(&s->system.gas_phase, s->variable_names[k]);
    double v;
    if (ps && tmc_species_get_double(&ps->species, TMC_PROPERTY_ABSOLUTE_TOLERANCE, &v, NULL) == TMC_OK)
      s->absolute_tolerance[k] = v;
  }

  /* Fingerprint of what generated kernels depend on */
  {
    unsigned long long h = TMC_FNV_OFFSET;
    hash_int(&h, n_vars);
    for (int k = 0; k < n_vars; ++k)
      hash_text(&h, s->variable_names[k]);
    hash_int(&h, n_rxn);
    for (int r = 0; r < n_rxn; ++r)
    {
      const tmc_process_plan *p = &s->processes[r];
      hash_int(&h, (int)s->reactions[r].rate_constant.type);
      hash_int(&h, p->n_reactants);
      for (int i = 0; i < p->n_reactants; ++i)
        hash_int(&h, p->reactants[i]);
      hash_int(&h, p->n_products);
      for (int i = 0; i < p->n_products; ++i)
      {
        hash_int(&h, p->products[i]);
        hash_double(&h, p->yields[i]);
      }
    }
    hash_int(&h, s->lu_pattern->nnz);
    for (int e = 0; e < s->lu_pattern->nnz; ++e)
    {
      hash_int(&h, s->lu_pattern->row_of[e]);
      hash_int(&h, s->lu_pattern->col_of[e]);
    }
    s->fingerprint = h;
  }

  if (b->kernels)
  {
    if (b->kernels->width != s->width || b->kernels->fingerprint != s->fingerprint)
    {
      tmc_error_set(err, TMC_ERROR_SOLVER, TMC_SOLVER_KERNEL_MISMATCH, "Generated kernels (width %d, fingerprint %llx) do not match the solver (width %d, fingerprint %llx)", b->kernels->width, b->kernels->fingerprint, s->width, s->fingerprint);
      free(new_of);
      goto fail;
    }
    s->kernels = b->kernels;
  }

  free(new_of);
  for (int r = 0; r < n_rxn; ++r)
  {
    tmc_ivec_free(&rids[r]);
    tmc_ivec_free(&pids[r]);
  }
  free(rids);
  free(pids);
  free(declared);
  return s;

fail:
  for (int r = 0; r < n_rxn; ++r)
  {
    tmc_ivec_free(&rids[r]);
    tmc_ivec_free(&pids[r]);
  }
  free(rids);
  free(pids);
  free(declared);
  tmc_solver_free(s);
  return NULL;
}

/* ---- solver accessors ---- */

void tmc_solver_set_rosenbrock_parameters(tmc_solver *s, const tmc_rosenbrock_parameters *p)
{
  s->ros = *p;
}

const tmc_rosenbrock_parameters *tmc_solver_rosenbrock_parameters(const tmc_solver *s)
{
  return &s->ros;
}

void tmc_solver_set_backward_euler_parameters(tmc_solver *s, const tmc_backward_euler_parameters *p)
{
  s->be = *p;
}

tmc_method tmc_solver_method(const tmc_solver *s)
{
  return s->method;
}

bool tmc_solver_uses_kernels(const tmc_solver *s)
{
  return s->kernels != NULL;
}

int tmc_solver_num_variables(const tmc_solver *s)
{
  return s->n_vars;
}

int tmc_solver_num_reactions(const tmc_solver *s)
{
  return s->n_reactions;
}

int tmc_solver_num_custom_parameters(const tmc_solver *s)
{
  return s->n_custom;
}

const char *tmc_solver_variable_name(const tmc_solver *s, int index)
{
  return s->variable_names[index];
}

const char *tmc_solver_custom_parameter_label(const tmc_solver *s, int index)
{
  return s->custom_labels[index];
}

int tmc_solver_width(const tmc_solver *s)
{
  return s->width;
}

const tmc_system *tmc_solver_system(const tmc_solver *s)
{
  return &s->system;
}

const tmc_reaction *tmc_solver_reactions(const tmc_solver *s, int *n)
{
  if (n)
    *n = s->n_reactions;
  return s->reactions;
}

const tmc_sparse *tmc_solver_lu_pattern(const tmc_solver *s)
{
  return s->lu_pattern;
}

const tmc_lu *tmc_solver_lu(const tmc_solver *s)
{
  return s->lu;
}

unsigned long long tmc_solver_fingerprint(const tmc_solver *s)
{
  return s->fingerprint;
}

tmc_result tmc_solver_solve(tmc_solver *s, double time_step, tmc_state *state)
{
  tmc_result r;
  if (s->method == TMC_METHOD_BACKWARD_EULER)
    r = tmc_backward_euler_solve(s, time_step, state);
  else
    r = tmc_rosenbrock_solve(s, time_step, state);
  tmc_state_clamp(state);
  return r;
}
