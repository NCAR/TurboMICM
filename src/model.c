#include "chem.h"
#include "external_internal.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char *join(const char *label, const char *suffix)
{
    size_t n=strlen(label), k=strlen(suffix);
    char *s=malloc(n+k+1);
    if (s) { memcpy(s,label,n); memcpy(s+n,suffix,k+1); }
    return s;
}

const char *chem_status_string(chem_status s)
{
    static const char *names[]={"Converged","Invalid argument","Allocation failure","Unknown name","Size mismatch",
        "Singular matrix","Convergence exceeded max steps","Step size too small","NaN detected","Inf detected",
        "Constraint Initialization Failed","Accepting unconverged integration"};
    return (unsigned)s < sizeof(names)/sizeof(*names) ? names[s] : "Unknown status";
}

size_t chem_species_index(const chem_model *m, const char *name)
{
    if (!name) return CHEM_NONE;
    for (size_t i=0; i<m->n_species; ++i)
        if (!strcmp(m->mechanism->species[m->var_to_spec[i]].name,name)) return i;
    return CHEM_NONE;
}
size_t chem_custom_index(const chem_model *m, const char *name)
{
    if (!name) return CHEM_NONE;
    for (size_t i=0; i<m->n_custom; ++i) if (!strcmp(m->custom_names[i],name)) return i;
    return CHEM_NONE;
}

void chem_model_destroy(chem_model *m)
{
    for (size_t i=0; i<m->n_custom; ++i) free(m->custom_names ? m->custom_names[i] : NULL);
    for (size_t i=0; i<m->n_jobs; ++i) {
        free(m->jobs[i].reactants); free(m->jobs[i].reactant_flat);
        free(m->jobs[i].products); free(m->jobs[i].product_flat); free(m->jobs[i].yields);
    }
    free(m->jobs); free(m->custom_names); free(m->custom_offset); free(m->rates);
    free(m->constraint_rows); free(m->constraint_offsets);
    free(m->spec_to_var); free(m->var_to_spec); free(m->algebraic); free(m->jacobian_pattern);
    chem_lu_destroy(&m->lu);
    chem_external_destroy(m);
    memset(m,0,sizeof(*m));
}

chem_status chem_model_create(chem_model *m, const chem_mechanism *s, chem_lu_algorithm algorithm)
{
    memset(m,0,sizeof(*m));
    if (!s || !s->n_species || !s->species || (s->n_reactions && !s->reactions) ||
        (s->n_constraints && !s->constraints) || (s->n_parameters && !s->parameters)) return CHEM_INVALID;
    chem_status status=CHEM_MEMORY;
    unsigned char *pattern=NULL;
    status=chem_external_prepare(m,&s);
    if (status!=CHEM_OK) goto fail;
    status=CHEM_MEMORY;
    m->mechanism=s; m->n_reactions=s->n_reactions;
    m->spec_to_var=malloc(s->n_species*sizeof(size_t));
    m->var_to_spec=malloc(s->n_species*sizeof(size_t));
    m->custom_offset=calloc(s->n_reactions+1,sizeof(size_t));
    m->constraint_rows=calloc(s->n_constraints+1,sizeof(size_t));
    m->constraint_offsets=calloc(s->n_constraints+1,sizeof(size_t));
    m->rates=calloc(s->n_reactions+1,sizeof(chem_rate));
    m->custom_names=calloc(2*s->n_reactions+s->n_constraints+s->n_parameters+1,sizeof(char *));
    if (!m->spec_to_var || !m->var_to_spec || !m->custom_offset || !m->rates || !m->custom_names || !m->constraint_rows || !m->constraint_offsets) goto fail;
    for (size_t i=0; i<s->n_species; ++i) {
        if (!s->species[i].name) { status=CHEM_INVALID; goto fail; }
        for (size_t j=0; j<i; ++j) if (!strcmp(s->species[i].name,s->species[j].name)) { status=CHEM_INVALID; goto fail; }
        if (s->species[i].type == CHEM_SPECIES_VARIABLE) {
            m->spec_to_var[i]=m->n_species; m->var_to_spec[m->n_species++]=i;
        } else m->spec_to_var[i]=CHEM_NONE;
    }
    size_t n=m->n_species, njobs=0;
    if (!n || n>SIZE_MAX/n) { status=CHEM_INVALID; goto fail; }
    pattern=calloc(n*n,1); m->algebraic=calloc(n,1);
    if (!pattern || !m->algebraic) goto fail;
    for (size_t r=0; r<s->n_reactions; ++r) {
        const chem_reaction *rx=&s->reactions[r];
        m->rates[r]=rx->rate;
        status=chem_rate_prepare(&m->rates[r]);
        if (status!=CHEM_OK) goto fail;
        m->custom_offset[r]=m->n_custom;
        size_t nc=rx->rate.type==CHEM_SURFACE ? 2 : rx->rate.type==CHEM_USER ? 1 : 0;
        if (nc && !rx->rate.label) { status=CHEM_INVALID; goto fail; }
        for (size_t j=0; j<nc; ++j) {
            const char *suffix=nc==1 ? "" : j==0 ? ".effective radius [m]" : ".particle number concentration [# m-3]";
            char *label=join(rx->rate.label,suffix);
            if (!label) { status=CHEM_MEMORY; goto fail; }
            m->custom_names[m->n_custom++]=label;
        }
        for (size_t j=0; j<rx->n_products; ++j)
            if (rx->products[j].species>=s->n_species) { status=CHEM_INVALID; goto fail; }
        for (size_t j=0; j<rx->n_reactants; ++j)
            if (rx->reactants[j]>=s->n_species) { status=CHEM_INVALID; goto fail; }
        for (size_t j=0; j<rx->n_reactants; ++j) {
            size_t col=m->spec_to_var[rx->reactants[j]];
            if (col==CHEM_NONE) continue;
            ++njobs;
            for (size_t k=0; k<rx->n_reactants; ++k) {
                size_t row=m->spec_to_var[rx->reactants[k]];
                if (row!=CHEM_NONE) pattern[row*n+col]=1;
            }
            for (size_t k=0; k<rx->n_products; ++k) {
                size_t row=m->spec_to_var[rx->products[k].species];
                if (row!=CHEM_NONE) pattern[row*n+col]=1;
            }
        }
    }
    m->custom_offset[s->n_reactions]=m->n_custom;
    for (size_t i=0;i<s->n_parameters;++i) {
        if (!s->parameters[i].name || chem_custom_index(m,s->parameters[i].name)!=CHEM_NONE) { status=CHEM_INVALID; goto fail; }
        char *name=join(s->parameters[i].name,""); if (!name) { status=CHEM_MEMORY; goto fail; }
        m->custom_names[m->n_custom++]=name;
    }
    chem_external_pattern(m,pattern);
    if (!s->preserve_order) {
        size_t *perm=malloc(n*sizeof(size_t)), *original=malloc(n*sizeof(size_t));
        unsigned char *reordered=malloc(n*n);
        if (!perm || !original || !reordered) {
            free(perm); free(original); free(reordered); status=CHEM_MEMORY; goto fail;
        }
        status=chem_markowitz(n,pattern,perm);
        if (status!=CHEM_OK) { free(perm); free(original); free(reordered); goto fail; }
        memcpy(original,m->var_to_spec,n*sizeof(size_t));
        for (size_t i=0; i<n; ++i) {
            m->var_to_spec[i]=original[perm[i]];
            m->spec_to_var[m->var_to_spec[i]]=i;
            for (size_t j=0; j<n; ++j) reordered[i*n+j]=pattern[perm[i]*n+perm[j]];
        }
        memcpy(pattern,reordered,n*n); free(perm); free(original); free(reordered);
    }
    for (size_t k=0; k<s->n_constraints; ++k) {
        const chem_constraint *constraint=&s->constraints[k];
        status=chem_constraint_validate(constraint,s->n_species);
        if (status!=CHEM_OK) goto fail;
        size_t row=m->spec_to_var[constraint->algebraic_species];
        if (row==CHEM_NONE || m->algebraic[row]) { status=CHEM_INVALID; goto fail; }
        m->algebraic[row]=1; m->constraint_rows[k]=row;
        memset(pattern+row*n,0,n);
        size_t count=constraint->n_reactants+constraint->n_products;
        if (constraint->type==CHEM_CUSTOM_CONSTRAINT) count=constraint->n_dependencies;
        for (size_t j=0; j<count; ++j) {
            size_t sp=constraint->type==CHEM_CUSTOM_CONSTRAINT ? constraint->dependencies[j] :
                j<constraint->n_reactants ? constraint->reactants[j].species : constraint->products[j-constraint->n_reactants].species;
            size_t col=m->spec_to_var[sp];
            if (col==CHEM_NONE) { status=CHEM_INVALID; goto fail; }
            pattern[row*n+col]=1;
        }
        m->constraint_offsets[k]=CHEM_NONE;
        if (constraint->type==CHEM_EQUILIBRIUM_CONSTRAINT) {
            m->constraint_offsets[k]=m->n_custom;
            char *label=join(constraint->name, "");
            if (!label) { status=CHEM_MEMORY; goto fail; }
            m->custom_names[m->n_custom++]=label;
        }
    }
    /* External callbacks finalize every declared process element, even when a
     * constraint masks that row during evaluation. Keep those storage slots. */
    chem_external_pattern(m,pattern);
    for (size_t i=0; i<n; ++i) pattern[i*n+i]=1;
    status=chem_lu_create(&m->lu,n,pattern,algorithm);
    if (status!=CHEM_OK) goto fail;
    status=CHEM_MEMORY;
    m->jobs=calloc(njobs+1,sizeof(chem_jac_job));
    if (!m->jobs) goto fail;
    /* Stable independent-variable/reaction ordering, including duplicated reactants. */
    for (size_t col=0; col<n; ++col) for (size_t r=0; r<s->n_reactions; ++r) {
        const chem_reaction *rx=&s->reactions[r];
        for (size_t occurrence=0; occurrence<rx->n_reactants; ++occurrence) {
            if (m->spec_to_var[rx->reactants[occurrence]]!=col) continue;
            chem_jac_job *job=&m->jobs[m->n_jobs++];
            job->reaction=r; job->independent=col; job->diagonal=m->lu.map[col*n+col];
            job->reactants=malloc((rx->n_reactants+1)*sizeof(size_t));
            job->reactant_flat=malloc((rx->n_reactants+1)*sizeof(size_t));
            job->products=malloc((rx->n_products+1)*sizeof(size_t));
            job->product_flat=malloc((rx->n_products+1)*sizeof(size_t));
            job->yields=malloc((rx->n_products+1)*sizeof(double));
            if (!job->reactants || !job->reactant_flat || !job->products || !job->product_flat || !job->yields) goto fail;
            int skipped=0;
            for (size_t j=0; j<rx->n_reactants; ++j) {
                size_t row=m->spec_to_var[rx->reactants[j]];
                if (row==CHEM_NONE) continue;
                if (row==col && !skipped) { skipped=1; continue; }
                size_t k=job->n_reactants++;
                job->reactants[k]=row; job->reactant_flat[k]=m->lu.map[row*n+col];
            }
            for (size_t j=0; j<rx->n_products; ++j) {
                size_t row=m->spec_to_var[rx->products[j].species];
                if (row==CHEM_NONE) continue;
                size_t k=job->n_products++;
                job->products[k]=row; job->product_flat[k]=m->lu.map[row*n+col]; job->yields[k]=rx->products[j].yield;
            }
        }
    }
    m->forcing=s->forcing && !s->n_constraints ? s->forcing : chem_forcing;
    m->jacobian=s->jacobian && !s->n_constraints ? s->jacobian : chem_jacobian;
    m->jacobian_pattern=pattern; pattern=NULL;
    status=chem_external_finalize(m);
    if (status!=CHEM_OK) goto fail;
    return CHEM_OK;
fail:
    free(pattern); chem_model_destroy(m); return status;
}

void chem_forcing(const chem_model *m, size_t cells, const double *restrict rates,
                  const double *restrict y, double *restrict f)
{
    for (size_t r=0; r<m->n_reactions; ++r) {
        const chem_reaction *rx=&m->mechanism->reactions[r];
        #pragma omp simd
        for (size_t c=0; c<cells; ++c) {
            double rate=rates[r*cells+c];
            for (size_t j=0; j<rx->n_reactants; ++j) {
                size_t v=m->spec_to_var[rx->reactants[j]];
                if (v!=CHEM_NONE) rate*=y[v*cells+c];
            }
            for (size_t j=0; j<rx->n_reactants; ++j) {
                size_t v=m->spec_to_var[rx->reactants[j]];
                if (v!=CHEM_NONE && !m->algebraic[v]) f[v*cells+c]-=rate;
            }
            for (size_t j=0; j<rx->n_products; ++j) {
                size_t v=m->spec_to_var[rx->products[j].species];
                if (v!=CHEM_NONE && !m->algebraic[v]) f[v*cells+c]+=rx->products[j].yield*rate;
            }
        }
    }
}

void chem_jacobian(const chem_model *m, size_t cells, const double *restrict rates,
                   const double *restrict y, double *restrict jac)
{
    for (size_t j=0; j<m->n_jobs; ++j) {
        const chem_jac_job *job=&m->jobs[j];
        #pragma omp simd
        for (size_t c=0; c<cells; ++c) {
            double rate=rates[job->reaction*cells+c];
            for (size_t k=0; k<job->n_reactants; ++k) rate*=y[job->reactants[k]*cells+c];
            for (size_t k=0; k<job->n_reactants; ++k)
                if (!m->algebraic[job->reactants[k]]) jac[job->reactant_flat[k]*cells+c]+=rate;
            if (!m->algebraic[job->independent]) jac[job->diagonal*cells+c]+=rate;
            for (size_t k=0; k<job->n_products; ++k) {
                size_t flat=job->product_flat[k];
                /* Sparse flat ids are compiled before any constraint rows are substituted. */
                if (!m->algebraic[job->products[k]]) jac[flat*cells+c]-=job->yields[k]*rate;
            }
        }
    }
}

chem_status chem_model_set_lambda(chem_model *model, const char *name, chem_rate_function function, void *context)
{
    if (!name || !function) return CHEM_INVALID;
    for (size_t i=0; i<model->n_reactions; ++i) {
        chem_rate *r=&model->rates[i];
        if (r->type==CHEM_LAMBDA && r->label && !strcmp(r->label,name)) {
            r->p.lambda.function=function; r->p.lambda.context=context; return CHEM_OK;
        }
    }
    return CHEM_NOT_FOUND;
}

chem_status chem_model_create_checked(chem_model *m, const chem_mechanism *s, chem_lu_algorithm algorithm, chem_error *error)
{
    chem_error local={CHEM_OK,{0}};
    if (!error) error=&local;
    memset(m,0,sizeof(*m)); memset(error,0,sizeof(*error));
    if (s && s->reactions) for (size_t i=0; i<s->n_reactions; ++i) {
        const chem_reaction *rx=&s->reactions[i];
        if (rx->rate.type!=CHEM_SURFACE) continue;
        const char *name="unknown";
        if (rx->n_reactants && rx->reactants[0]<s->n_species) name=s->species[rx->reactants[0]].name;
        if (rx->rate.p.surface.diffusion<=0) {
            error->code=CHEM_NOT_FOUND;
            snprintf(error->message,sizeof(error->message),"Diffusion coefficient for species '%s' is not defined",name);
            return error->code;
        }
        if (rx->rate.p.surface.molecular_weight<=0) {
            error->code=CHEM_NOT_FOUND;
            snprintf(error->message,sizeof(error->message),"Species: '%s' Property: 'molecular weight [kg mol-1]'",name);
            return error->code;
        }
    }
    error->code=chem_model_create(m,s,algorithm);
    if (error->code!=CHEM_OK) snprintf(error->message,sizeof(error->message),"%s",chem_status_string(error->code));
    return error->code;
}
