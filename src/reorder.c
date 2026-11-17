#include "reorder.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

void chem_permutation_destroy(ChemPermutation *p)
{
    if (!p) return;
    free(p->old_to_new); free(p->new_to_old); *p=(ChemPermutation){0};
}
ChemReorderStatus chem_permutation_init(ChemPermutation *p,size_t n,const unsigned char *pattern)
{
    if (!p || !n || !pattern) return CHEM_REORDER_INVALID;
    if (n>SIZE_MAX/n || n>SIZE_MAX/sizeof(size_t)) return CHEM_REORDER_NO_MEMORY;
    unsigned char *edges=malloc(n*n),*active=malloc(n);
    ChemPermutation copy={.size=n,.old_to_new=malloc(n*sizeof(size_t)),.new_to_old=malloc(n*sizeof(size_t))};
    if (!edges || !active || !copy.old_to_new || !copy.new_to_old) {
        free(edges); free(active); chem_permutation_destroy(&copy); return CHEM_REORDER_NO_MEMORY;
    }
    memcpy(edges,pattern,n*n); memset(active,1,n);
    for (size_t i=0;i<n;++i) edges[i*n+i]=1;
    for (size_t index=0;index<n;++index) {
        size_t best=SIZE_MAX,best_cost=SIZE_MAX;
        for (size_t i=0;i<n;++i) if (active[i]) {
            size_t row=0,col=0;
            for (size_t j=0;j<n;++j) if (active[j]) { row+=edges[i*n+j]!=0; col+=edges[j*n+i]!=0; }
            size_t cost=(row-1)*(col-1);
            if (best==SIZE_MAX || cost<best_cost) { best=i; best_cost=cost; }
        }
        copy.new_to_old[index]=best; copy.old_to_new[best]=index; active[best]=0;
        for (size_t row=0;row<n;++row) if (active[row] && edges[row*n+best])
            for (size_t col=0;col<n;++col) if (active[col] && edges[best*n+col]) edges[row*n+col]=1;
    }
    free(edges); free(active); chem_permutation_destroy(p); *p=copy; return CHEM_REORDER_OK;
}
void chem_reordered_constraints_destroy(ChemReorderedConstraints *c)
{
    if (!c) return;
    free(c->rows); free(c->terms); *c=(ChemReorderedConstraints){0};
}
ChemReorderStatus chem_reordered_constraints_init(ChemReorderedConstraints *out,
    const ChemConstraintSet *source,const ChemPermutation *p)
{
    if (!out || !source || !p || !p->old_to_new || !p->new_to_old || source->variables!=p->size ||
        chem_constraint_validate(source)!=CHEM_CONSTRAINT_OK) return CHEM_REORDER_INVALID;
    size_t total=0;
    for (size_t i=0;i<source->count;++i) {
        const ChemConstraint *c=source->constraints+i;
        if (c->type==CHEM_CONSTRAINT_CUSTOM) return CHEM_REORDER_INVALID;
        if (c->reactant_count>SIZE_MAX-total) return CHEM_REORDER_NO_MEMORY;
        total+=c->reactant_count;
        if (c->product_count>SIZE_MAX-total) return CHEM_REORDER_NO_MEMORY;
        total+=c->product_count;
    }
    if (source->count>SIZE_MAX/sizeof(ChemConstraint) || total>SIZE_MAX/sizeof(ChemConstraintTerm)) return CHEM_REORDER_NO_MEMORY;
    ChemReorderedConstraints copy={0};
    copy.rows=source->count ? malloc(source->count*sizeof(*copy.rows)) : NULL;
    copy.terms=total ? malloc(total*sizeof(*copy.terms)) : NULL;
    if ((source->count && !copy.rows) || (total && !copy.terms)) {
        chem_reordered_constraints_destroy(&copy); return CHEM_REORDER_NO_MEMORY;
    }
    size_t offset=0;
    for (size_t i=0;i<source->count;++i) {
        const ChemConstraint *c=source->constraints+i; copy.rows[i]=*c;
        copy.rows[i].row=p->old_to_new[c->row];
        copy.rows[i].reactants=c->reactant_count ? copy.terms+offset : NULL;
        for (size_t j=0;j<c->reactant_count;++j) {
            copy.terms[offset]=c->reactants[j]; copy.terms[offset++].species=p->old_to_new[c->reactants[j].species];
        }
        copy.rows[i].products=c->product_count ? copy.terms+offset : NULL;
        for (size_t j=0;j<c->product_count;++j) {
            copy.terms[offset]=c->products[j]; copy.terms[offset++].species=p->old_to_new[c->products[j].species];
        }
    }
    copy.set=(ChemConstraintSet){source->variables,source->count,copy.rows};
    chem_reordered_constraints_destroy(out); *out=copy; return CHEM_REORDER_OK;
}
ChemReorderStatus chem_permutation_for_system(ChemPermutation *p,const ChemMechanism *m,const ChemConstraintSet *constraints)
{
    if (!p || !m || !m->species_count || (constraints &&
        (constraints->variables!=m->species_count || chem_constraint_validate(constraints)!=CHEM_CONSTRAINT_OK)))
        return CHEM_REORDER_INVALID;
    size_t n=m->species_count;
    if (n>SIZE_MAX/n) return CHEM_REORDER_NO_MEMORY;
    unsigned char *pattern=calloc(n*n,1);
    if (!pattern) return CHEM_REORDER_NO_MEMORY;
    if (!chem_jacobian_pattern(m,pattern)) { free(pattern); return CHEM_REORDER_INVALID; }
    if (constraints) chem_constraint_sparsity(constraints,pattern);
    ChemReorderStatus status=chem_permutation_init(p,n,pattern);
    free(pattern); return status;
}
void chem_reordered_mechanism_destroy(ChemReorderedMechanism *m)
{
    if (!m) return;
    free(m->reactions); free(m->reactants); free(m->products); *m=(ChemReorderedMechanism){0};
}
ChemReorderStatus chem_reordered_mechanism_init(ChemReorderedMechanism *out,const ChemMechanism *source,
    const ChemPermutation *p)
{
    if (!out || !source || !p || !p->old_to_new || !p->new_to_old || source->species_count!=p->size ||
        (source->reaction_count && !source->reactions)) return CHEM_REORDER_INVALID;
    size_t nr=0,np=0;
    for (size_t i=0;i<source->reaction_count;++i) {
        const ChemReaction *r=source->reactions+i;
        if (r->reactant_count>SIZE_MAX-r->reactant_begin || r->product_count>SIZE_MAX-r->product_begin) return CHEM_REORDER_INVALID;
        size_t re=r->reactant_begin+r->reactant_count,pe=r->product_begin+r->product_count;
        if (re>nr) nr=re;
        if (pe>np) np=pe;
    }
    if ((nr && !source->reactants) || (np && !source->products)) return CHEM_REORDER_INVALID;
    if (nr>SIZE_MAX/sizeof(size_t) || np>SIZE_MAX/sizeof(ChemProduct) ||
        source->reaction_count>SIZE_MAX/sizeof(ChemReaction)) return CHEM_REORDER_NO_MEMORY;
    for (size_t i=0;i<source->reaction_count;++i) {
        const ChemReaction *r=source->reactions+i;
        for (size_t j=0;j<r->reactant_count;++j) if (source->reactants[r->reactant_begin+j]>=p->size) return CHEM_REORDER_INVALID;
        for (size_t j=0;j<r->product_count;++j) if (source->products[r->product_begin+j].species>=p->size) return CHEM_REORDER_INVALID;
    }
    ChemReorderedMechanism copy={0};
    copy.reactions=source->reaction_count ? malloc(source->reaction_count*sizeof(*copy.reactions)) : NULL;
    copy.reactants=nr ? calloc(nr,sizeof(*copy.reactants)) : NULL;
    copy.products=np ? calloc(np,sizeof(*copy.products)) : NULL;
    if ((source->reaction_count && !copy.reactions) || (nr && !copy.reactants) || (np && !copy.products)) {
        chem_reordered_mechanism_destroy(&copy); return CHEM_REORDER_NO_MEMORY;
    }
    if (source->reaction_count) memcpy(copy.reactions,source->reactions,source->reaction_count*sizeof(*copy.reactions));
    for (size_t i=0;i<source->reaction_count;++i) {
        const ChemReaction *r=source->reactions+i;
        for (size_t j=0;j<r->reactant_count;++j) {
            size_t index=r->reactant_begin+j; copy.reactants[index]=p->old_to_new[source->reactants[index]];
        }
        for (size_t j=0;j<r->product_count;++j) {
            size_t index=r->product_begin+j; copy.products[index]=source->products[index];
            copy.products[index].species=p->old_to_new[source->products[index].species];
        }
    }
    copy.mechanism=(ChemMechanism){p->size,source->reaction_count,copy.reactions,copy.reactants,copy.products};
    chem_reordered_mechanism_destroy(out); *out=copy; return CHEM_REORDER_OK;
}
