/* SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"
#include <string.h>

static char *copy_string(const char *a, const char *suffix)
{
    size_t n = strlen(a), m = strlen(suffix);
    if (n > SIZE_MAX - m - 1) return NULL;
    char *s = malloc(n + m + 1);
    if (s) { memcpy(s, a, n); memcpy(s + n, suffix, m + 1); }
    return s;
}
void chem_mechanism_free(ChemMechanism *m)
{
    if (!m) return;
    if (m->species) for (size_t i = 0; i < m->species_count; ++i) free((void *)m->species[i].name);
    if (m->processes) for (size_t i = 0; i < m->reaction_count; ++i) {
        free(m->processes[i].reactants); free(m->processes[i].parameters);
        free(m->processes[i].products); free((void *)m->processes[i].rate.label);
    }
    if (m->custom_names) for (size_t i = 0; i < m->custom_count; ++i) free(m->custom_names[i]);
    chem_constraints_free(m->constraints,m->constraint_count); free(m->algebraic);
    free(m->species); free(m->variable_species); free(m->processes); free(m->custom_names); free(m);
}
ChemStatus chem_mechanism_create(size_t ns, const ChemSpecies *species, size_t nr,
                                 const ChemReaction *reactions, ChemMechanism **out)
{
    if (!out || (ns && !species) || (nr && !reactions)) return CHEM_INVALID;
    *out = NULL;
    ChemMechanism *m = calloc(1, sizeof(*m));
    if (!m) return CHEM_MEMORY;
    ChemStatus status = CHEM_MEMORY;
    size_t *map = chem_alloc(ns, sizeof(*map));
    m->species_count = ns; m->reaction_count = nr;
    m->species = chem_alloc(ns, sizeof(*m->species));
    m->variable_species = chem_alloc(ns, sizeof(*m->variable_species));
    m->processes = chem_alloc(nr, sizeof(*m->processes));
    if (!map || !m->species || !m->variable_species || !m->processes) goto fail;
    for (size_t i = 0; i < ns; ++i) {
        status = CHEM_INVALID;
        if (!species[i].name || !*species[i].name) goto fail;
        for (size_t j = 0; j < i; ++j) if (!strcmp(species[i].name, species[j].name)) goto fail;
        m->species[i] = species[i];
        m->species[i].name = copy_string(species[i].name, "");
        status = CHEM_MEMORY;
        if (!m->species[i].name) goto fail;
        map[i] = SIZE_MAX;
        if (!species[i].parameter) {
            map[i] = m->variable_count;
            m->variable_species[m->variable_count++] = i;
        }
    }
    /* Counting sort is stable, preserving MICM's reaction accumulation order. */
    for (size_t i = 0; i < nr; ++i) {
        status = CHEM_INVALID;
        unsigned type = (unsigned)reactions[i].rate.type;
        if (type >= CHEM_RATE_TYPE_COUNT) goto fail;
        ++m->offsets[type + 1];
    }
    for (size_t t = 1; t <= CHEM_RATE_TYPE_COUNT; ++t) m->offsets[t] += m->offsets[t-1];
    size_t next[CHEM_RATE_TYPE_COUNT];
    memcpy(next, m->offsets, sizeof(next));
    for (size_t i = 0; i < nr; ++i) {
        const ChemReaction *r = reactions + i;
        ChemProcess *p = m->processes + next[r->rate.type]++;
        status = CHEM_INVALID;
        if ((r->reactant_count && !r->reactants) || (r->product_count && !r->products)) goto fail;
        p->rate = r->rate;
        p->rate.label = copy_string(r->rate.label ? r->rate.label : "", "");
        p->reactants = chem_alloc(r->reactant_count, sizeof(*p->reactants));
        p->parameters = chem_alloc(r->reactant_count, sizeof(*p->parameters));
        p->products = chem_alloc(r->product_count, sizeof(*p->products));
        status = CHEM_MEMORY;
        if (!p->rate.label || !p->reactants || !p->parameters || !p->products) goto fail;
        status = chem_rate_prepare(&p->rate);
        if (status != CHEM_OK) goto fail;
        for (size_t j = 0; j < r->reactant_count; ++j) {
            size_t s = r->reactants[j];
            status = CHEM_NOT_FOUND;
            if (s >= ns) goto fail;
            if (map[s] == SIZE_MAX) p->parameters[p->parameter_count++] = s;
            else p->reactants[p->reactant_count++] = map[s];
        }
        for (size_t j = 0; j < r->product_count; ++j) {
            size_t s = r->products[j].species;
            status = CHEM_NOT_FOUND;
            if (s >= ns) goto fail;
            if (map[s] != SIZE_MAX) p->products[p->product_count++] = (ChemProduct){map[s], r->products[j].yield};
        }
    }
    for (size_t i = 0; i < nr; ++i) {
        ChemRate *r = &m->processes[i].rate;
        r->custom_index = m->custom_count;
        m->custom_count += r->type == CHEM_USER ? 1 : r->type == CHEM_SURFACE ? 2 : 0;
    }
    m->custom_names = chem_alloc(m->custom_count, sizeof(*m->custom_names));
    status = CHEM_MEMORY;
    if (!m->custom_names) goto fail;
    for (size_t i = 0; i < nr; ++i) {
        ChemRate *r = &m->processes[i].rate;
        size_t j = r->custom_index;
        if (r->type == CHEM_USER) {
            m->custom_names[j] = copy_string(r->label, "");
            if (!m->custom_names[j]) goto fail;
        } else if (r->type == CHEM_SURFACE) {
            m->custom_names[j] = copy_string(r->label, ".effective radius [m]");
            m->custom_names[j+1] = copy_string(r->label, ".particle number concentration [# m-3]");
            if (!m->custom_names[j] || !m->custom_names[j+1]) goto fail;
        }
    }
    m->reaction_custom_count=m->custom_count;
    m->algebraic=chem_alloc(m->variable_count,1);
    if (!m->algebraic) { status=CHEM_MEMORY; goto fail; }
    free(map); *out = m; return CHEM_OK;
fail:
    free(map); chem_mechanism_free(m); return status;
}
size_t chem_species_count(const ChemMechanism *m) { return m->variable_count; }
size_t chem_reaction_count(const ChemMechanism *m) { return m->reaction_count; }
size_t chem_custom_count(const ChemMechanism *m) { return m->custom_count; }
const char *chem_species_name(const ChemMechanism *m, size_t v)
{ return v < m->variable_count ? m->species[m->variable_species[v]].name : NULL; }
const char *chem_custom_name(const ChemMechanism *m, size_t p)
{ return p < m->custom_count ? m->custom_names[p] : NULL; }
size_t chem_rate_offset(const ChemMechanism *m, ChemRateType t)
{ return (unsigned)t <= CHEM_RATE_TYPE_COUNT ? m->offsets[t] : m->reaction_count; }
ChemStatus chem_species_index(const ChemMechanism *m, const char *name, size_t *out)
{
    if (!m || !name || !out) return CHEM_INVALID;
    for (size_t i = 0; i < m->variable_count; ++i)
        if (!strcmp(chem_species_name(m, i), name)) { *out = i; return CHEM_OK; }
    return CHEM_NOT_FOUND;
}
ChemStatus chem_custom_index(const ChemMechanism *m, const char *name, size_t *out)
{
    if (!m || !name || !out) return CHEM_INVALID;
    for (size_t i = 0; i < m->custom_count; ++i)
        if (!strcmp(m->custom_names[i], name)) { *out = i; return CHEM_OK; }
    return CHEM_NOT_FOUND;
}
ChemStatus chem_set_rate_callback(ChemMechanism *m, const char *name, ChemRateFunction fn, void *ctx)
{
    if (!m || !name || !fn) return CHEM_INVALID;
    for (size_t i = m->offsets[CHEM_CALLBACK]; i < m->reaction_count; ++i) {
        ChemRate *r = &m->processes[i].rate;
        if (!strcmp(r->label, name)) {
            r->data.callback.function = fn; r->data.callback.context = ctx; return CHEM_OK;
        }
    }
    return CHEM_NOT_FOUND;
}

#include "chem_solver.h"
unsigned char *chem_mechanism_pattern(const ChemMechanism *m)
{
    size_t n=m->variable_count;
    if (n && n>SIZE_MAX/n) return NULL;
    unsigned char *pattern=chem_alloc(n*n,1);
    if (!pattern) return NULL;
    for (size_t r=0;r<m->reaction_count;++r) {
        const ChemProcess *p=&m->processes[r];
        for (size_t j=0;j<p->reactant_count;++j) {
            size_t col=p->reactants[j];
            for (size_t k=0;k<p->reactant_count;++k) pattern[p->reactants[k]*n+col]=1;
            for (size_t k=0;k<p->product_count;++k) pattern[p->products[k].species*n+col]=1;
        }
    }
    for (size_t row=0;row<n;++row) if (m->algebraic[row]) {
        memset(pattern+row*n,0,n);
        for (size_t col=0;col<n;++col) pattern[row*n+col]=(unsigned char)chem_constraint_has_element(m,row,col);
    }
    return pattern;
}
ChemStatus chem_mechanism_reorder(ChemMechanism *m)
{
    if (!m) return CHEM_INVALID;
    size_t n=m->variable_count;
    unsigned char *pattern=chem_mechanism_pattern(m);
    size_t *perm=chem_alloc(n,sizeof(size_t)), *inverse=chem_alloc(n,sizeof(size_t)), *species=chem_alloc(n,sizeof(size_t));
    ChemStatus status=CHEM_MEMORY;
    if (!pattern || !perm || !inverse || !species) goto finish;
    status=chem_markowitz(n,pattern,perm);
    if (status) goto finish;
    for (size_t i=0;i<n;++i) { inverse[perm[i]]=i; species[i]=m->variable_species[perm[i]]; }
    memcpy(m->variable_species,species,n*sizeof(size_t));
    for (size_t r=0;r<m->reaction_count;++r) {
        ChemProcess *p=&m->processes[r];
        for (size_t j=0;j<p->reactant_count;++j) p->reactants[j]=inverse[p->reactants[j]];
        for (size_t j=0;j<p->product_count;++j) p->products[j].species=inverse[p->products[j].species];
    }
    memset(m->algebraic,0,n);
    for (size_t i=0;i<m->constraint_count;++i) {
        ChemBoundConstraint *c=&m->constraints[i];
        c->row=inverse[c->row]; m->algebraic[c->row]=1;
        for (size_t j=0;j<c->definition.reactant_count;++j) c->reactants[j]=inverse[c->reactants[j]];
        for (size_t j=0;j<c->definition.product_count;++j) c->products[j]=inverse[c->products[j]];
    }
finish:
    free(pattern); free(perm); free(inverse); free(species); return status;
}

int chem_reactions_use_species(const ChemMechanism *m,const char *name)
{
    for (size_t i=0;i<m->reaction_count;++i) {
        const ChemProcess *p=&m->processes[i];
        for (size_t j=0;j<p->reactant_count;++j)
            if (!strcmp(chem_species_name(m,p->reactants[j]),name)) return 1;
        for (size_t j=0;j<p->product_count;++j)
            if (!strcmp(chem_species_name(m,p->products[j].species),name)) return 1;
        for (size_t j=0;j<p->parameter_count;++j)
            if (!strcmp(m->species[p->parameters[j]].name,name)) return 1;
    }
    return 0;
}
int chem_kinetic_has_element(const ChemMechanism *m,size_t row,size_t column)
{
    for (size_t i=0;i<m->reaction_count;++i) {
        const ChemProcess *p=&m->processes[i];int independent=0,dependent=0;
        for (size_t j=0;j<p->reactant_count;++j) {
            independent|=p->reactants[j]==column;dependent|=p->reactants[j]==row;
        }
        for (size_t j=0;j<p->product_count;++j) dependent|=p->products[j].species==row;
        if (independent && dependent) return 1;
    }
    return 0;
}

const ChemRate *chem_rate_definition(const ChemMechanism *m,size_t reaction)
{ return reaction<m->reaction_count?&m->processes[reaction].rate:NULL; }
size_t chem_parameterized_reactant_count(const ChemMechanism *m,size_t reaction)
{ return reaction<m->reaction_count?m->processes[reaction].parameter_count:0; }
double chem_reaction_multiplier(const ChemMechanism *m,size_t reaction,const ChemConditions *conditions)
{
    if (reaction>=m->reaction_count) return 1;
    const ChemProcess *p=&m->processes[reaction];double multiplier=1;
    for (size_t i=0;i<p->parameter_count;++i) {
        const ChemSpecies *species=&m->species[p->parameters[i]];
        multiplier*=species->parameter(conditions,species->parameter_context);
    }
    return multiplier;
}
