/* The solver builder: from a system, reactions, constraints, and external
 * models to the compiled plans a solver and its states use (MICM's
 * SolverBuilder::Build, step for step). */
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
  int n_constraints;
  tmc_constraint *constraints;
  int n_external;
  tmc_external_model *externals;
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
  for (int i = 0; i < b->n_constraints; ++i)
    tmc_constraint_free(&b->constraints[i]);
  free(b->constraints);
  free(b->externals);
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

void tmc_builder_add_constraint(tmc_builder *b, const tmc_constraint *constraint)
{
  b->constraints = tmc_realloc(b->constraints, (size_t)(b->n_constraints + 1) * sizeof *b->constraints);
  b->constraints[b->n_constraints++] = tmc_constraint_copy(constraint);
}

void tmc_builder_add_external_model(tmc_builder *b, const tmc_external_model *model)
{
  b->externals = tmc_realloc(b->externals, (size_t)(b->n_external + 1) * sizeof *b->externals);
  b->externals[b->n_external++] = *model;
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

typedef char tmc_name[TMC_NAME_MAX];

static int name_index(tmc_name *names, int n, const char *name)
{
  for (int i = 0; i < n; ++i)
    if (strcmp(names[i], name) == 0)
      return i;
  return -1;
}

/* A growable list of names. */
typedef struct names_list
{
  tmc_name *v;
  int n, cap;
} names_list;

static void names_push(names_list *l, const char *name)
{
  if (l->n == l->cap)
  {
    l->cap = l->cap ? 2 * l->cap : 16;
    l->v = tmc_realloc(l->v, (size_t)l->cap * sizeof *l->v);
  }
  tmc_strlcpy(l->v[l->n++], name, TMC_NAME_MAX);
}

static bool names_contains(const names_list *l, const char *name)
{
  return name_index(l->v, l->n, name) >= 0;
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
}

/* Dense 0/1 pattern with set semantics. */
typedef struct pattern
{
  int n;
  unsigned char *p;
} pattern;

static void pattern_add(pattern *pt, int row, int col)
{
  pt->p[(size_t)row * pt->n + col] = 1;
}

static bool pattern_has(const pattern *pt, int row, int col)
{
  return pt->p[(size_t)row * pt->n + col] != 0;
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
  *h = tmc_fnv1a(*h, bytes, 8);
}

void tmc_solver_free(tmc_solver *s)
{
  if (!s)
    return;
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
  for (int i = 0; i < s->n_constraints; ++i)
  {
    tmc_constraint_free(&s->constraints[i]);
    if (s->constraint_plans)
    {
      free(s->constraint_plans[i].dependencies);
      free(s->constraint_plans[i].flat_ids);
    }
  }
  free(s->constraints);
  free(s->constraint_plans);
  free(s->algebraic_ids);
  free(s->is_algebraic);
  free(s->mass_diagonal);
  free(s->externals);
  free(s->external_constraint_active);
  free(s->algebraic_flat_ids);
  tmc_system_free(&s->system);
  tmc_reactions_free(s->reactions, s->n_reactions);
  free(s);
}

/* Everything the build allocates before it can fail. */
typedef struct build_scratch
{
  names_list declared;      /* merged state variable names in declaration order */
  tmc_ivec *rids, *pids;    /* per sorted reaction: reactant / product declared ids */
  int n_rxn;
  int *new_of;              /* declared id -> solver id */
  names_list params;        /* custom parameter labels */
  names_list constraint_params; /* per constraint parameter, in order */
  tmc_ivec alg_ids;         /* algebraic variable ids, in insertion order */
  tmc_ivec *ext_var_ids;    /* per external model: its algebraic ids */
  pattern jac;              /* Jacobian pattern in solver ids */
} build_scratch;

static void scratch_free(build_scratch *sc)
{
  free(sc->declared.v);
  if (sc->rids)
    for (int r = 0; r < sc->n_rxn; ++r)
    {
      tmc_ivec_free(&sc->rids[r]);
      tmc_ivec_free(&sc->pids[r]);
    }
  free(sc->rids);
  free(sc->pids);
  free(sc->new_of);
  free(sc->params.v);
  free(sc->constraint_params.v);
  tmc_ivec_free(&sc->alg_ids);
  free(sc->ext_var_ids);
  free(sc->jac.p);
}

tmc_solver *tmc_builder_build(const tmc_builder *b, tmc_error *err)
{
  tmc_error_clear(err);
  if (!b->has_system)
  {
    tmc_error_set(err, TMC_ERROR_SOLVER, TMC_SOLVER_MISSING_CHEMICAL_SYSTEM, "Missing chemical system.");
    return NULL;
  }
  for (int i = 0; i < b->n_reactions; ++i)
  {
    const tmc_reaction *rx = &b->reactions[i];
    if (tmc_reaction_validate(rx, err) != TMC_OK)
      return NULL;
    if (rx->rate_constant.type == TMC_RATE_SURFACE)
    {
      const tmc_surface_parameters *sp = &rx->rate_constant.p.surface;
      if (!sp->species.has_diffusion_coefficient)
      {
        tmc_error_set(err, TMC_ERROR_SPECIES, TMC_SPECIES_PROPERTY_NOT_FOUND, "Diffusion coefficient for species '%s' is not defined", sp->species.species.name);
        return NULL;
      }
      double mw;
      if (tmc_species_get_double(&sp->species.species, TMC_PROPERTY_MOLECULAR_WEIGHT, &mw, err) != TMC_OK)
        return NULL;
    }
  }

  build_scratch sc;
  memset(&sc, 0, sizeof sc);
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
  s->n_constraints = b->n_constraints;
  s->constraints = tmc_calloc((size_t)(b->n_constraints > 0 ? b->n_constraints : 1), sizeof *s->constraints);
  for (int i = 0; i < b->n_constraints; ++i)
    s->constraints[i] = tmc_constraint_copy(&b->constraints[i]);
  s->n_external = b->n_external;
  s->externals = tmc_calloc((size_t)(b->n_external > 0 ? b->n_external : 1), sizeof *s->externals);
  memcpy(s->externals, b->externals, (size_t)b->n_external * sizeof *s->externals);
  s->external_constraint_active = tmc_calloc((size_t)(b->n_external > 0 ? b->n_external : 1), 1);

  /* 1. The state variable names: the gas phase's, then each external model's
   *    (in std::set order), all in declaration order. */
  {
    int n_gas = tmc_system_state_size(&s->system);
    tmc_name *gas = tmc_malloc((size_t)(n_gas > 0 ? n_gas : 1) * sizeof *gas);
    tmc_system_unique_names(&s->system, gas, n_gas);
    for (int i = 0; i < n_gas; ++i)
      names_push(&sc.declared, gas[i]);
    free(gas);
    for (int m = 0; m < s->n_external; ++m)
    {
      tmc_name *names;
      int n = tmc_external_names(&s->externals[m], s->externals[m].state_variable_names, &names);
      for (int i = 0; i < n; ++i)
        names_push(&sc.declared, names[i]);
      free(names);
    }
  }
  const int n_vars = sc.declared.n;
  if (n_vars == 0)
  {
    tmc_error_set(err, TMC_ERROR_SOLVER, TMC_SOLVER_MISSING_CHEMICAL_SPECIES, "Provided chemical system contains no species.");
    goto fail;
  }
  s->n_vars = n_vars;

  /* 2. Reactant and product ids in declared order (parameterized species skipped) */
  sc.n_rxn = s->n_reactions;
  sc.rids = tmc_calloc((size_t)(sc.n_rxn > 0 ? sc.n_rxn : 1), sizeof *sc.rids);
  sc.pids = tmc_calloc((size_t)(sc.n_rxn > 0 ? sc.n_rxn : 1), sizeof *sc.pids);
  for (int r = 0; r < sc.n_rxn; ++r)
  {
    const tmc_reaction *rx = &s->reactions[r];
    for (int i = 0; i < rx->n_reactants; ++i)
    {
      if (rx->reactants[i].parameterized)
        continue;
      int id = name_index(sc.declared.v, n_vars, rx->reactants[i].name);
      if (id < 0)
      {
        tmc_error_set(err, TMC_ERROR_PROCESS, TMC_PROCESS_REACTANT_DOES_NOT_EXIST, "%s", rx->reactants[i].name);
        goto fail;
      }
      tmc_ivec_push(&sc.rids[r], id);
    }
    for (int i = 0; i < rx->n_products; ++i)
    {
      if (rx->products[i].species.parameterized)
        continue;
      int id = name_index(sc.declared.v, n_vars, rx->products[i].species.name);
      if (id < 0)
      {
        tmc_error_set(err, TMC_ERROR_PROCESS, TMC_PROCESS_PRODUCT_DOES_NOT_EXIST, "%s", rx->products[i].species.name);
        goto fail;
      }
      tmc_ivec_push(&sc.pids[r], id);
    }
  }

  /* 3. Unused species: every state variable must appear in a reaction or be
   *    used by an external model. */
  if (!b->ignore_unused_species)
  {
    unsigned char *used = tmc_calloc((size_t)n_vars, 1);
    for (int r = 0; r < sc.n_rxn; ++r)
    {
      for (int i = 0; i < sc.rids[r].n; ++i)
        used[sc.rids[r].v[i]] = 1;
      for (int i = 0; i < sc.pids[r].n; ++i)
        used[sc.pids[r].v[i]] = 1;
    }
    for (int m = 0; m < s->n_external; ++m)
    {
      const tmc_external_model *em = &s->externals[m];
      int (*callbacks[3])(void *, char (*)[TMC_NAME_MAX], int) = { em->species_used, em->constraint_species_dependencies, em->constraint_algebraic_variable_names };
      for (int k = 0; k < 3; ++k)
      {
        tmc_name *names;
        int n = tmc_external_names(em, callbacks[k], &names);
        for (int i = 0; i < n; ++i)
        {
          int id = name_index(sc.declared.v, n_vars, names[i]);
          if (id >= 0)
            used[id] = 1;
        }
        free(names);
      }
    }
    /* a constraint's algebraic species and dependencies count as used, as in MICM */
    for (int k = 0; k < s->n_constraints; ++k)
    {
      const tmc_constraint *c = &s->constraints[k];
      int id = name_index(sc.declared.v, n_vars, c->algebraic_species.name);
      if (id >= 0)
        used[id] = 1;
      int nd = tmc_constraint_dependencies(c, NULL, 0);
      tmc_name *deps = tmc_malloc((size_t)(nd > 0 ? nd : 1) * sizeof *deps);
      tmc_constraint_dependencies(c, deps, nd);
      for (int i = 0; i < nd; ++i)
      {
        id = name_index(sc.declared.v, n_vars, deps[i]);
        if (id >= 0)
          used[id] = 1;
      }
      free(deps);
    }
    char msg[1024] = "Unused species in chemical system:";
    bool any = false;
    for (int i = 0; i < n_vars; ++i)
      if (!used[i])
      {
        any = true;
        size_t len = strlen(msg);
        snprintf(msg + len, sizeof msg - len, " '%s'", sc.declared.v[i]);
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

  /* 4. The species order: declaration order, then the Markowitz reordering of
   *    the kinetic Jacobian pattern (and the external process patterns). */
  {
    int *perm = tmc_malloc((size_t)n_vars * sizeof(int));
    for (int i = 0; i < n_vars; ++i)
      perm[i] = i;
    if (b->reorder_state)
    {
      unsigned char *pat = tmc_calloc((size_t)n_vars * (size_t)n_vars, 1);
      for (int r = 0; r < sc.n_rxn; ++r)
        for (int a = 0; a < sc.rids[r].n; ++a)
        {
          int ind = sc.rids[r].v[a];
          for (int d = 0; d < sc.rids[r].n; ++d)
            pat[(size_t)sc.rids[r].v[d] * n_vars + ind] = 1;
          for (int d = 0; d < sc.pids[r].n; ++d)
            pat[(size_t)sc.pids[r].v[d] * n_vars + ind] = 1;
        }
      if (s->n_external)
      {
        tmc_name_map declared_map = { n_vars, sc.declared.v };
        for (int m = 0; m < s->n_external; ++m)
        {
          const tmc_external_model *em = &s->externals[m];
          if (!em->nonzero_jacobian_elements)
            continue;
          int n = em->nonzero_jacobian_elements(em->context, &declared_map, NULL, 0);
          if (n > 0)
          {
            tmc_index_pair *e = tmc_malloc((size_t)n * sizeof *e);
            n = em->nonzero_jacobian_elements(em->context, &declared_map, e, n);
            for (int i = 0; i < n; ++i)
              pat[(size_t)e[i].row * n_vars + e[i].col] = 1;
            free(e);
          }
        }
      }
      tmc_markowitz_reorder(n_vars, pat, perm);
      free(pat);
    }
    sc.new_of = tmc_malloc((size_t)n_vars * sizeof(int));
    s->variable_names = tmc_malloc((size_t)n_vars * sizeof *s->variable_names);
    for (int k = 0; k < n_vars; ++k)
    {
      sc.new_of[perm[k]] = k;
      memcpy(s->variable_names[k], sc.declared.v[perm[k]], TMC_NAME_MAX);
    }
    free(perm);
  }
  s->variable_map.n = n_vars;
  s->variable_map.names = s->variable_names;

  /* 5. Process plans in solver ids, and the custom parameter labels */
  s->processes = tmc_calloc((size_t)(sc.n_rxn > 0 ? sc.n_rxn : 1), sizeof *s->processes);
  for (int r = 0; r < sc.n_rxn; ++r)
  {
    tmc_process_plan *p = &s->processes[r];
    const tmc_reaction *rx = &s->reactions[r];
    p->n_reactants = sc.rids[r].n;
    p->reactants = tmc_malloc((size_t)(p->n_reactants > 0 ? p->n_reactants : 1) * sizeof(int));
    for (int i = 0; i < p->n_reactants; ++i)
      p->reactants[i] = sc.new_of[sc.rids[r].v[i]];
    p->n_products = sc.pids[r].n;
    p->products = tmc_malloc((size_t)(p->n_products > 0 ? p->n_products : 1) * sizeof(int));
    p->yields = tmc_malloc((size_t)(p->n_products > 0 ? p->n_products : 1) * sizeof(double));
    int k = 0;
    for (int i = 0; i < rx->n_products; ++i)
    {
      if (rx->products[i].species.parameterized)
        continue;
      p->products[k] = sc.new_of[sc.pids[r].v[k]];
      p->yields[k] = rx->products[i].coefficient;
      ++k;
    }
    p->n_parameterized = 0;
    for (int i = 0; i < rx->n_reactants; ++i)
      if (rx->reactants[i].parameterized && p->n_parameterized < 4)
        p->parameterized[p->n_parameterized++] = rx->reactants[i];
    p->custom_index = -1;
    tmc_name labels[2];
    int n_labels = tmc_rate_constant_custom_parameter_labels(&rx->rate_constant, labels);
    if (n_labels > 0)
      p->custom_index = sc.params.n;
    for (int i = 0; i < n_labels; ++i)
    {
      if (names_contains(&sc.params, labels[i]))
      {
        tmc_error_set(err, TMC_ERROR_SOLVER, TMC_SOLVER_DUPLICATE_PARAMETER, "Duplicate parameter labels detected:\n  - %s (from reaction)\n", labels[i]);
        goto fail;
      }
      names_push(&sc.params, labels[i]);
    }
  }
  for (int m = 0; m < s->n_external; ++m)
  {
    tmc_name *names;
    int n = tmc_external_names(&s->externals[m], s->externals[m].state_parameter_names, &names);
    for (int i = 0; i < n; ++i)
    {
      if (names_contains(&sc.params, names[i]))
      {
        tmc_error_set(err, TMC_ERROR_SOLVER, TMC_SOLVER_DUPLICATE_PARAMETER, "Duplicate parameter labels detected:\n  - %s (from external_model)\n", names[i]);
        free(names);
        goto fail;
      }
      names_push(&sc.params, names[i]);
    }
    free(names);
  }

  /* 6. Jacobian jobs: (independent variable, process) pairs, stably ordered by variable */
  {
    tmc_ivec job_ind = { 0 }, job_proc = { 0 };
    for (int r = 0; r < sc.n_rxn; ++r)
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

  /* 7. The kinetic Jacobian pattern in solver ids, plus the external process patterns */
  sc.jac.n = n_vars;
  sc.jac.p = tmc_calloc((size_t)n_vars * (size_t)n_vars, 1);
  for (int r = 0; r < sc.n_rxn; ++r)
  {
    const tmc_process_plan *p = &s->processes[r];
    for (int a = 0; a < p->n_reactants; ++a)
    {
      int ind = p->reactants[a];
      for (int d = 0; d < p->n_reactants; ++d)
        pattern_add(&sc.jac, p->reactants[d], ind);
      for (int d = 0; d < p->n_products; ++d)
        pattern_add(&sc.jac, p->products[d], ind);
    }
  }
  pattern ext_processes = { n_vars, tmc_calloc((size_t)n_vars * (size_t)n_vars, 1) };
  for (int m = 0; m < s->n_external; ++m)
  {
    const tmc_external_model *em = &s->externals[m];
    if (!em->nonzero_jacobian_elements)
      continue;
    s->has_external_processes = true;
    int n = em->nonzero_jacobian_elements(em->context, &s->variable_map, NULL, 0);
    if (n > 0)
    {
      tmc_index_pair *e = tmc_malloc((size_t)n * sizeof *e);
      n = em->nonzero_jacobian_elements(em->context, &s->variable_map, e, n);
      for (int i = 0; i < n; ++i)
      {
        pattern_add(&sc.jac, e[i].row, e[i].col);
        pattern_add(&ext_processes, e[i].row, e[i].col);
      }
      free(e);
    }
  }
  for (int m = 0; m < s->n_external; ++m)
    if (s->externals[m].add_forcing_terms)
      s->has_external_processes = true;

  /* 8. Constraints: algebraic variables, their parameters, and their pattern */
  s->mass_diagonal = tmc_malloc((size_t)n_vars * sizeof(double));
  s->is_algebraic = tmc_calloc((size_t)n_vars, 1);
  for (int i = 0; i < n_vars; ++i)
    s->mass_diagonal[i] = 1.0;
  if (s->n_constraints)
  {
    s->constraint_plans = tmc_calloc((size_t)s->n_constraints, sizeof *s->constraint_plans);
    /* parameter labels, made unique as MICM's SetUniqueParameterNames does */
    names_list used = { 0 };
    for (int k = 0; k < s->n_constraints; ++k)
    {
      tmc_name labels[1];
      int n = tmc_constraint_parameter_names(&s->constraints[k], labels, 1);
      for (int i = 0; i < n; ++i)
      {
        tmc_name label;
        memcpy(label, labels[i], TMC_NAME_MAX);
        int count = 0;
        while (names_contains(&used, label))
          snprintf(label, TMC_NAME_MAX, "%s_%d", labels[i], ++count);
        names_push(&used, label);
        names_push(&sc.constraint_params, label);
      }
    }
    free(used.v);
    for (int i = 0; i < sc.constraint_params.n; ++i)
    {
      if (names_contains(&sc.params, sc.constraint_params.v[i]))
      {
        tmc_error_set(err, TMC_ERROR_SOLVER, TMC_SOLVER_DUPLICATE_PARAMETER, "Duplicate parameter name: %s", sc.constraint_params.v[i]);
        free(ext_processes.p);
        goto fail;
      }
      names_push(&sc.params, sc.constraint_params.v[i]);
    }
    int param_k = 0;
    for (int k = 0; k < s->n_constraints; ++k)
    {
      const tmc_constraint *c = &s->constraints[k];
      tmc_constraint_plan *plan = &s->constraint_plans[k];
      plan->row = name_index(s->variable_names, n_vars, c->algebraic_species.name);
      if (plan->row < 0)
      {
        tmc_error_set(err, TMC_ERROR_CONSTRAINT, TMC_CONSTRAINT_UNKNOWN_SPECIES, "Constraint '%s' targets unknown algebraic species '%s'", c->name, c->algebraic_species.name);
        free(ext_processes.p);
        goto fail;
      }
      if (s->is_algebraic[plan->row])
      {
        tmc_error_set(err, TMC_ERROR_CONSTRAINT, TMC_CONSTRAINT_INVALID_STOICHIOMETRY, "Multiple constraints map to the same algebraic species row '%s'", c->algebraic_species.name);
        free(ext_processes.p);
        goto fail;
      }
      s->is_algebraic[plan->row] = 1;
      tmc_ivec_push(&sc.alg_ids, plan->row);
      plan->n_dependencies = tmc_constraint_dependencies(c, NULL, 0);
      plan->dependencies = tmc_malloc((size_t)(plan->n_dependencies > 0 ? plan->n_dependencies : 1) * sizeof(int));
      plan->flat_ids = tmc_malloc((size_t)(plan->n_dependencies > 0 ? plan->n_dependencies : 1) * sizeof(int));
      tmc_name *deps = tmc_malloc((size_t)(plan->n_dependencies > 0 ? plan->n_dependencies : 1) * sizeof *deps);
      tmc_constraint_dependencies(c, deps, plan->n_dependencies);
      for (int i = 0; i < plan->n_dependencies; ++i)
      {
        plan->dependencies[i] = name_index(s->variable_names, n_vars, deps[i]);
        if (plan->dependencies[i] < 0)
        {
          tmc_error_set(err, TMC_ERROR_CONSTRAINT, TMC_CONSTRAINT_UNKNOWN_SPECIES, "Constraint '%s' depends on unknown species '%s'", c->name, deps[i]);
          free(deps);
          free(ext_processes.p);
          goto fail;
        }
      }
      free(deps);
      plan->parameter_index = -1;
      if (c->type == TMC_CONSTRAINT_EQUILIBRIUM)
        plan->parameter_index = name_index(sc.params.v, sc.params.n, sc.constraint_params.v[param_k++]);
    }
    /* algebraic rows take no kinetic terms; the constraints' own pattern replaces them */
    for (int i = 0; i < n_vars; ++i)
      if (s->is_algebraic[i])
        memset(sc.jac.p + (size_t)i * n_vars, 0, (size_t)n_vars);
    for (int k = 0; k < s->n_constraints; ++k)
    {
      const tmc_constraint_plan *plan = &s->constraint_plans[k];
      pattern_add(&sc.jac, plan->row, plan->row);
      for (int i = 0; i < plan->n_dependencies; ++i)
        pattern_add(&sc.jac, plan->row, plan->dependencies[i]);
    }
  }
  /* external constraints */
  {
    names_list seen_params = { 0 }, seen_init = { 0 };
    for (int m = 0; m < s->n_external; ++m)
    {
      const tmc_external_model *em = &s->externals[m];
      tmc_name *alg;
      int n_alg = tmc_external_names(em, em->constraint_algebraic_variable_names, &alg);
      if (n_alg == 0)
      {
        free(alg);
        continue;
      }
      s->external_constraint_active[m] = 1;
      for (int i = 0; i < n_alg; ++i)
      {
        int id = name_index(s->variable_names, n_vars, alg[i]);
        if (id < 0)
        {
          tmc_error_set(err, TMC_ERROR_CONSTRAINT, TMC_CONSTRAINT_UNKNOWN_SPECIES, "External model constraint targets unknown algebraic species '%s'", alg[i]);
          free(alg);
          free(seen_params.v);
          free(seen_init.v);
          free(ext_processes.p);
          goto fail;
        }
        if (s->is_algebraic[id])
        {
          tmc_error_set(err, TMC_ERROR_CONSTRAINT, TMC_CONSTRAINT_DUPLICATE_ALGEBRAIC_SPECIES, "Multiple constraints map to the same algebraic species row '%s'", alg[i]);
          free(alg);
          free(seen_params.v);
          free(seen_init.v);
          free(ext_processes.p);
          goto fail;
        }
        s->is_algebraic[id] = 1;
        tmc_ivec_push(&sc.alg_ids, id);
      }
      free(alg);
      tmc_name *names;
      int n = tmc_external_names(em, em->constraint_state_parameter_names, &names);
      for (int i = 0; i < n; ++i)
      {
        if (names_contains(&seen_params, names[i]))
        {
          tmc_error_set(err, TMC_ERROR_CONSTRAINT, TMC_CONSTRAINT_DUPLICATE_PARAMETER, "Duplicate external constraint parameter name across models: %s", names[i]);
          free(names);
          free(seen_params.v);
          free(seen_init.v);
          free(ext_processes.p);
          goto fail;
        }
        if (names_contains(&sc.params, names[i]))
        {
          tmc_error_set(err, TMC_ERROR_SOLVER, TMC_SOLVER_DUPLICATE_PARAMETER, "Duplicate parameter name: %s", names[i]);
          free(names);
          free(seen_params.v);
          free(seen_init.v);
          free(ext_processes.p);
          goto fail;
        }
        names_push(&seen_params, names[i]);
        names_push(&sc.params, names[i]);
      }
      free(names);
      n = tmc_external_names(em, em->initialize_constraint_parameter_names, &names);
      for (int i = 0; i < n; ++i)
      {
        if (names_contains(&seen_init, names[i]))
        {
          tmc_error_set(err, TMC_ERROR_CONSTRAINT, TMC_CONSTRAINT_DUPLICATE_PARAMETER, "Duplicate external initialize constraint parameter name across models: %s", names[i]);
          free(names);
          free(seen_params.v);
          free(seen_init.v);
          free(ext_processes.p);
          goto fail;
        }
        if (names_contains(&sc.params, names[i]))
        {
          tmc_error_set(err, TMC_ERROR_SOLVER, TMC_SOLVER_DUPLICATE_PARAMETER, "Duplicate parameter name: %s", names[i]);
          free(names);
          free(seen_params.v);
          free(seen_init.v);
          free(ext_processes.p);
          goto fail;
        }
        names_push(&seen_init, names[i]);
        names_push(&sc.params, names[i]);
      }
      free(names);
      if (em->nonzero_constraint_jacobian_elements)
      {
        int ne = em->nonzero_constraint_jacobian_elements(em->context, &s->variable_map, NULL, 0);
        if (ne > 0)
        {
          tmc_index_pair *e = tmc_malloc((size_t)ne * sizeof *e);
          ne = em->nonzero_constraint_jacobian_elements(em->context, &s->variable_map, e, ne);
          for (int i = 0; i < ne; ++i)
            pattern_add(&sc.jac, e[i].row, e[i].col);
          free(e);
        }
      }
      for (int i = 0; i < n_vars; ++i)
        if (s->is_algebraic[i])
          pattern_add(&sc.jac, i, i);
    }
    free(seen_params.v);
    free(seen_init.v);
  }
  /* algebraic ids ascending (MICM keeps them in a std::set) */
  s->n_algebraic = sc.alg_ids.n;
  s->algebraic_ids = tmc_malloc((size_t)(s->n_algebraic > 0 ? s->n_algebraic : 1) * sizeof(int));
  {
    int k = 0;
    for (int i = 0; i < n_vars; ++i)
      if (s->is_algebraic[i])
      {
        s->algebraic_ids[k++] = i;
        s->mass_diagonal[i] = 0.0;
      }
  }
  /* external process elements in algebraic rows stay in the pattern and are cleared after evaluation */
  tmc_ivec alg_flat_pairs = { 0 };
  if (s->n_algebraic)
    for (int i = 0; i < n_vars; ++i)
      if (s->is_algebraic[i])
        for (int j = 0; j < n_vars; ++j)
          if (pattern_has(&ext_processes, i, j))
          {
            pattern_add(&sc.jac, i, j);
            tmc_ivec_push(&alg_flat_pairs, i);
            tmc_ivec_push(&alg_flat_pairs, j);
          }
  free(ext_processes.p);

  /* 9. The Jacobian pattern (with the diagonal) and the LU */
  {
    tmc_sparse_builder *sb = tmc_sparse_builder_create(n_vars);
    for (int i = 0; i < n_vars; ++i)
      for (int j = 0; j < n_vars; ++j)
        if (pattern_has(&sc.jac, i, j))
          tmc_sparse_builder_with_element(sb, i, j, NULL);
    for (int i = 0; i < n_vars; ++i)
      tmc_sparse_builder_with_element(sb, i, i, NULL);
    s->jacobian_pattern = tmc_sparse_build(sb, TMC_SPARSE_CSR, s->width, true);
    tmc_sparse_builder_free(sb);
    s->lu = tmc_lu_create(s->jacobian_pattern, b->lu_algorithm, err);
    if (!s->lu)
    {
      tmc_ivec_free(&alg_flat_pairs);
      goto fail;
    }
    s->lu_pattern = s->lu->pattern;
  }
  s->n_algebraic_flat = alg_flat_pairs.n / 2;
  s->algebraic_flat_ids = tmc_malloc((size_t)(s->n_algebraic_flat > 0 ? s->n_algebraic_flat : 1) * sizeof(int));
  for (int i = 0; i < s->n_algebraic_flat; ++i)
    s->algebraic_flat_ids[i] = tmc_sparse_element(s->lu_pattern, alg_flat_pairs.v[2 * i], alg_flat_pairs.v[2 * i + 1]);
  tmc_ivec_free(&alg_flat_pairs);

  /* 10. The custom parameter labels, and the flat LU ids of the jobs and constraints */
  s->n_custom = sc.params.n;
  s->custom_labels = tmc_malloc((size_t)(s->n_custom > 0 ? s->n_custom : 1) * sizeof *s->custom_labels);
  for (int i = 0; i < s->n_custom; ++i)
    memcpy(s->custom_labels[i], sc.params.v[i], TMC_NAME_MAX);
  s->parameter_map.n = s->n_custom;
  s->parameter_map.names = s->custom_labels;
  {
    s->n_job_flat_total = s->n_job_dependent_total + s->n_jobs + s->n_job_products_total;
    s->job_flat_ids = tmc_malloc((size_t)(s->n_job_flat_total > 0 ? s->n_job_flat_total : 1) * sizeof(int));
    int d = 0, pr = 0, f = 0;
    for (int q = 0; q < s->n_jobs; ++q)
    {
      int ind = s->job_independent[q];
      for (int i = 0; i < s->job_n_dependent[q]; ++i, ++d)
        s->job_flat_ids[f++] = s->is_algebraic[s->job_dependent_ids[d]] ? -1 : tmc_sparse_element(s->lu_pattern, s->job_dependent_ids[d], ind);
      s->job_flat_ids[f++] = s->is_algebraic[ind] ? -1 : tmc_sparse_element(s->lu_pattern, ind, ind);
      for (int i = 0; i < s->job_n_products[q]; ++i, ++pr)
        s->job_flat_ids[f++] = s->is_algebraic[s->job_product_ids[pr]] ? -1 : tmc_sparse_element(s->lu_pattern, s->job_product_ids[pr], ind);
    }
    for (int k = 0; k < s->n_constraints; ++k)
    {
      tmc_constraint_plan *plan = &s->constraint_plans[k];
      for (int i = 0; i < plan->n_dependencies; ++i)
        plan->flat_ids[i] = tmc_sparse_element(s->lu_pattern, plan->row, plan->dependencies[i]);
    }
  }

  /* 11. External models bind their indices */
  for (int m = 0; m < s->n_external; ++m)
  {
    const tmc_external_model *em = &s->externals[m];
    if (em->finalize_process_setup)
      em->finalize_process_setup(em->context, &s->parameter_map, &s->variable_map, s->lu_pattern);
    if (em->finalize_constraint_setup && s->external_constraint_active[m])
      em->finalize_constraint_setup(em->context, &s->parameter_map, &s->variable_map, s->lu_pattern);
  }

  /* 12. Absolute tolerances: 1e-3 unless the species sets a property */
  s->absolute_tolerance = tmc_malloc((size_t)n_vars * sizeof(double));
  for (int k = 0; k < n_vars; ++k)
  {
    s->absolute_tolerance[k] = 1.0e-3;
    const tmc_phase_species *ps = tmc_phase_find(&s->system.gas_phase, s->variable_names[k]);
    double v;
    if (ps && tmc_species_get_double(&ps->species, TMC_PROPERTY_ABSOLUTE_TOLERANCE, &v, NULL) == TMC_OK)
      s->absolute_tolerance[k] = v;
  }

  /* 13. Fingerprint of what generated kernels depend on */
  {
    unsigned long long h = TMC_FNV_OFFSET;
    hash_int(&h, n_vars);
    for (int k = 0; k < n_vars; ++k)
      hash_text(&h, s->variable_names[k]);
    hash_int(&h, sc.n_rxn);
    for (int r = 0; r < sc.n_rxn; ++r)
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
    if (s->n_constraints || s->n_external || s->n_algebraic)
    {
      tmc_error_set(err, TMC_ERROR_SOLVER, TMC_SOLVER_KERNEL_MISMATCH, "Generated kernels cover kinetic mechanisms only; this solver has constraints or external models");
      goto fail;
    }
    if (b->kernels->width != s->width || b->kernels->fingerprint != s->fingerprint)
    {
      tmc_error_set(err, TMC_ERROR_SOLVER, TMC_SOLVER_KERNEL_MISMATCH, "Generated kernels (width %d, fingerprint %llx) do not match the solver (width %d, fingerprint %llx)", b->kernels->width, b->kernels->fingerprint, s->width, s->fingerprint);
      goto fail;
    }
    s->kernels = b->kernels;
  }

  scratch_free(&sc);
  return s;

fail:
  scratch_free(&sc);
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

int tmc_solver_num_constraints(const tmc_solver *s)
{
  return s->n_algebraic;
}

const double *tmc_solver_mass_matrix_diagonal(const tmc_solver *s)
{
  return s->mass_diagonal;
}

int tmc_solver_algebraic_variables(const tmc_solver *s, const int **ids)
{
  if (ids)
    *ids = s->algebraic_ids;
  return s->n_algebraic;
}

const tmc_name_map *tmc_solver_variable_map(const tmc_solver *s)
{
  return &s->variable_map;
}

const tmc_name_map *tmc_solver_parameter_map(const tmc_solver *s)
{
  return &s->parameter_map;
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
  if (s->n_external)
    tmc_external_initialize_constraint_parameters(s, state->variables, state->custom_rate_parameters);
  tmc_result r;
  if (s->method == TMC_METHOD_BACKWARD_EULER)
    r = tmc_backward_euler_solve(s, time_step, state);
  else
    r = tmc_rosenbrock_solve(s, time_step, state);
  tmc_state_clamp(state);
  return r;
}
