#include "model.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

void chem_model_system_destroy(ChemModelSystem *s)
{
    if (!s) return;
    chem_state_destroy(&s->state); free(s->components); free(s->constraint_rows);
    chem_sparse_matrix_destroy(&s->structure); chem_sparse_plan_destroy(s->plan);
    *s=(ChemModelSystem){0};
}
ChemModelStatus chem_model_system_init(ChemModelSystem *destination,const ChemStateConfig *base,
    const ChemOdeComponent *base_process,const unsigned char *base_pattern,
    const ChemConstraintSet *built_in,size_t model_count,const ChemExternalModel *models,
    size_t width,ChemSparseOrdering ordering)
{
    if (!destination || !base || (model_count && !models) ||
        (ordering!=CHEM_SPARSE_CSR && ordering!=CHEM_SPARSE_CSC)) return CHEM_MODEL_INVALID;
    size_t n=base->species_count, nc=built_in ? built_in->count : 0;
    for (size_t i=0;i<model_count;++i) {
        if ((models[i].jacobian_count && !models[i].jacobian_entries) ||
            (models[i].constraint_count && !models[i].constraints) ||
            models[i].state.species_count>SIZE_MAX-n || models[i].constraint_count>SIZE_MAX-nc)
            return CHEM_MODEL_INVALID;
        n+=models[i].state.species_count; nc+=models[i].constraint_count;
    }
    if (!n || n>SIZE_MAX/n || n*n>SIZE_MAX/sizeof(ChemSparseEntry) ||
        model_count==SIZE_MAX || model_count+1>SIZE_MAX/sizeof(ChemOdeComponent) ||
        model_count>SIZE_MAX/sizeof(ChemStateFragment) || nc>SIZE_MAX/sizeof(ChemConstraint) ||
        (built_in && (built_in->variables!=n || (built_in->count && !built_in->constraints)))) return CHEM_MODEL_INVALID;
    ChemModelSystem copy={0}; ChemModelStatus status=CHEM_MODEL_OK;
    ChemStateFragment *fragments=model_count ? malloc(model_count*sizeof(*fragments)) : NULL;
    double *defaults=malloc(n*sizeof(*defaults));
    copy.components=calloc(model_count+1,sizeof(*copy.components));
    copy.constraint_rows=nc ? malloc(nc*sizeof(*copy.constraint_rows)) : NULL;
    unsigned char *pattern=calloc(n*n,1);
    if ((model_count && !fragments) || !defaults || !copy.components || (nc && !copy.constraint_rows) || !pattern) {
        status=CHEM_MODEL_NO_MEMORY; goto done;
    }
    for (size_t i=0;i<n;++i) defaults[i]=.001;
    size_t component_count=0,ci=0;
    if (base_process) copy.components[component_count++]=*base_process;
    if (built_in && built_in->count) {
        memcpy(copy.constraint_rows,built_in->constraints,built_in->count*sizeof(*copy.constraint_rows)); ci=built_in->count;
    }
    for (size_t i=0;i<model_count;++i) {
        fragments[i]=models[i].state;
        if (!fragments[i].absolute_tolerances) fragments[i].absolute_tolerances=defaults;
        copy.components[component_count++]=models[i].process;
        if (models[i].constraint_count) {
            memcpy(copy.constraint_rows+ci,models[i].constraints,models[i].constraint_count*sizeof(*copy.constraint_rows));
            ci+=models[i].constraint_count;
        }
    }
    copy.constraints=(ChemConstraintSet){n,nc,copy.constraint_rows};
    if (chem_constraint_validate(&copy.constraints)!=CHEM_CONSTRAINT_OK) { status=CHEM_MODEL_INVALID; goto done; }
    ChemStateConfig config=*base;
    if (!config.absolute_tolerances) config.absolute_tolerances=defaults;
    ChemStateStatus state_status=nc ? chem_state_init_constrained(&copy.state,&config,model_count,fragments,&copy.constraints) :
                                       chem_state_init_composed(&copy.state,&config,model_count,fragments);
    if (state_status!=CHEM_STATE_OK) { status=CHEM_MODEL_STATE_ERROR; goto done; }
    if (base_pattern) memcpy(pattern,base_pattern,n*n);
    for (size_t i=0;i<model_count;++i) for (size_t j=0;j<models[i].jacobian_count;++j) {
        ChemSparseEntry e=models[i].jacobian_entries[j];
        if (e.row>=n || e.column>=n) { status=CHEM_MODEL_INVALID; goto done; }
        pattern[e.row*n+e.column]=1;
    }
    if (nc) chem_constraint_sparsity(&copy.constraints,pattern);
    for (size_t i=0;i<n;++i) pattern[i*n+i]=1;
    size_t count=0; for (size_t i=0;i<n*n;++i) count+=pattern[i]!=0;
    ChemSparseEntry *entries=count ? malloc(count*sizeof(*entries)) : NULL;
    if (count && !entries) { status=CHEM_MODEL_NO_MEMORY; goto done; }
    size_t pos=0;
    for (size_t i=0;i<n;++i) for (size_t j=0;j<n;++j)
        if (pattern[i*n+j]) entries[pos++]=(ChemSparseEntry){i,j};
    ChemMatrixStatus matrix_status=chem_sparse_matrix_init(&copy.structure,n,base->cells,width,ordering,entries,count,0);
    free(entries);
    if (matrix_status!=CHEM_MATRIX_OK) { status=CHEM_MODEL_NO_MEMORY; goto done; }
    copy.plan=chem_sparse_plan_from_pattern(n,pattern);
    if (!copy.plan) { status=CHEM_MODEL_NO_MEMORY; goto done; }
    copy.ode=(ChemCompositeOde){n,component_count,copy.components};
    chem_model_system_destroy(destination); *destination=copy; copy=(ChemModelSystem){0};
    for (size_t i=0;i<model_count;++i) if (models[i].finalize)
        models[i].finalize(models[i].process.context,&destination->state,&destination->structure);
done:
    free(pattern); free(fragments); free(defaults); chem_model_system_destroy(&copy); return status;
}
