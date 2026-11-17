#include "chem.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

/* An independent implementation of the reference aerosol transfer and
 * constraint models, using only the public extension interface. */
typedef struct { int source,product,solvent,ss,ps,sv,pv,finalized; double rate; } Transfer;
typedef struct { size_t count; int species[4]; double coefficients[4],constant; int parameter_term,parameter; } Expression;
typedef struct {
    ChemSpecies species[4]; ChemMechanism mechanism; ChemJacobianEntry pattern[4];
    ChemExternalConstraint external_constraint; ChemExternalModel model;
    ChemConstraint builtin; ChemTerm reactant,product; int dependencies[2];
    Transfer transfer; Expression expression; ChemState *state;
} Aerosol;

static int transfer_finalize(const ChemMechanism *m,const ChemLinearSystem *j,void *context)
{
    (void)m; Transfer *t=context; size_t n=j->species_count;
    t->ss=j->lookup[(size_t)t->source*n+(size_t)t->source];
    t->ps=j->lookup[(size_t)t->product*n+(size_t)t->source];
    if (t->solvent>=0) {
        t->sv=j->lookup[(size_t)t->source*n+(size_t)t->solvent];
        t->pv=j->lookup[(size_t)t->product*n+(size_t)t->solvent];
        if (t->sv<0 || t->pv<0) return CHEM_INVALID;
    }
    ++t->finalized; return t->ss<0 || t->ps<0 ? CHEM_INVALID : CHEM_OK;
}
static void transfer_forcing(const ChemState *s,const double *y,double *out,void *context)
{
    Transfer *t=context; size_t cells=chem_cell_count(s);
    for (size_t c=0;c<cells;++c) {
        double rate=t->rate*y[(size_t)t->source*cells+c];
        if (t->solvent>=0) rate*=y[(size_t)t->solvent*cells+c];
        out[(size_t)t->source*cells+c]-=rate; out[(size_t)t->product*cells+c]+=rate;
    }
}
static void transfer_jacobian(const ChemState *s,const double *y,ChemLinearSystem *j,void *context)
{
    (void)s; Transfer *t=context;
    for (size_t c=0;c<j->cells;++c) {
        double ds=t->rate*(t->solvent<0 ? 1 : y[(size_t)t->solvent*j->cells+c]);
        j->values[(size_t)t->ss*j->cells+c]+=ds; j->values[(size_t)t->ps*j->cells+c]-=ds;
        if (t->solvent>=0) {
            double dv=t->rate*y[(size_t)t->source*j->cells+c];
            j->values[(size_t)t->sv*j->cells+c]+=dv; j->values[(size_t)t->pv*j->cells+c]-=dv;
        }
    }
}
static double evaluate_expression(const ChemState *s,size_t cell,const double *y,size_t stride,double *gradient,void *context)
{
    Expression *e=context; double result=e->constant;
    for (size_t i=0;i<e->count;++i) {
        double coefficient=(int)i==e->parameter_term ? chem_parameter_values_const(s,(size_t)e->parameter)[cell] : e->coefficients[i];
        result+=coefficient*y[(size_t)e->species[i]*stride];
        if (gradient) gradient[e->species[i]]+=coefficient;
    }
    return result;
}
static int aerosol_init(Aerosol *a,size_t cells,int external,int builtin,int solvent,double rate)
{
    memset(a,0,sizeof(*a));
    a->species[0]=(ChemSpecies){"A_GAS",0,1e-3,0,0};
    a->species[1]=(ChemSpecies){"AEROSOL.A_AQ",0,1e-3,0,0};
    a->species[2]=(ChemSpecies){solvent ? "S" : "B",0,1e-3,0,0};
    a->mechanism=(ChemMechanism){builtin || solvent ? 3 : 2,0,0,a->species,NULL,NULL,NULL,NULL,NULL,298,101325};
    a->transfer=(Transfer){.source=0,.product=1,.solvent=solvent ? 2 : -1,.rate=rate};
    a->pattern[0]=(ChemJacobianEntry){0,0}; a->pattern[1]=(ChemJacobianEntry){1,0};
    a->pattern[2]=(ChemJacobianEntry){0,2}; a->pattern[3]=(ChemJacobianEntry){1,2};
    a->dependencies[0]=0; a->dependencies[1]=1;
    a->expression=(Expression){.count=2,.species={0,1},.coefficients={1,1},.constant=-1,.parameter_term=-1};
    a->external_constraint=(ChemExternalConstraint){"mass",1,a->dependencies,2,evaluate_expression,&a->expression};
    a->model=(ChemExternalModel){.jacobian_pattern=a->pattern,.jacobian_count=solvent ? 4 : 2,
        .constraints=external ? &a->external_constraint : NULL,.constraint_count=(size_t)external,
        .finalize=transfer_finalize,.add_forcing=transfer_forcing,.subtract_jacobian=transfer_jacobian,.context=&a->transfer};
    if (builtin) {
        a->reactant=(ChemTerm){2,1}; a->product=(ChemTerm){1,1};
        a->builtin=(ChemConstraint){"equilibrium",CHEM_EQUILIBRIUM_CONSTRAINT,1,&a->reactant,&a->product,1,1,0,5,0,0,0};
    }
    a->state=chem_state_create_extended(&a->mechanism,builtin ? &a->builtin : NULL,(size_t)builtin,&a->model,1,cells);
    return a->state!=NULL;
}
static void aerosol_free(Aerosol *a) { chem_state_free(a->state); }
static int aerosol_metadata(int external)
{
    Aerosol a; if (!aerosol_init(&a,1,external,0,0,.01)) return 0;
    int ok=chem_cell_count(a.state)==1 && a.transfer.finalized==1 && !chem_is_algebraic(a.state,0) &&
           chem_is_algebraic(a.state,1)==external && chem_find_constraint(a.state,"mass")== (external ? 0 : -1);
    aerosol_free(&a); return ok;
}
static int aerosol_conservation(size_t cells,int solvent,int steps,double dt)
{
    Aerosol a; if (!aerosol_init(&a,cells,1,0,solvent,.1)) return 0;
    for (size_t c=0;c<cells;++c) {
        chem_species_values(a.state,0)[c]=.9-.2*c;
        chem_species_values(a.state,1)[c]=.1+.2*c;
        if (solvent) chem_species_values(a.state,2)[c]=1;
    }
    int ok=chem_update_rates(a.state)==CHEM_OK;
    ChemRosenbrock p=chem_rosenbrock_four_stage_dae();
    for (int i=0;i<steps;++i) {
        ChemStats stats; ok &= chem_solve(a.state,&p,dt,&stats)==CHEM_OK;
        for (size_t c=0;c<cells;++c)
            ok &= fabs(chem_species_values(a.state,0)[c]+chem_species_values(a.state,1)[c]-1)<1e-4;
    }
    aerosol_free(&a); return ok;
}
static int algebraic_row(void)
{
    Aerosol a; if (!aerosol_init(&a,1,1,0,0,10)) return 0;
    chem_species_values(a.state,0)[0]=.7; chem_species_values(a.state,1)[0]=.2;
    double f[2],j[4];
    int ok=chem_update_rates(a.state)==CHEM_OK && chem_evaluate(a.state,0,f,j)==CHEM_OK;
    ok &= fabs(f[0]+7)<1e-12 && fabs(f[1]+.1)<1e-12 && j[0]==10 && j[2]==-1 && j[3]==-1;
    aerosol_free(&a); return ok;
}
static int combined_constraints(void)
{
    const ChemSpecies sp[]={{"A_GAS",.8,1e-3,0,0},{"B",.1,1e-3,0,0},{"C",.5,1e-3,0,0},{"A_AQ",.2,1e-3,0,0}};
    const int reactants[]={0}; const ChemProduct products[]={{1,1}};
    const ChemReaction reaction={0,1,0,1,0,{.kind=CHEM_ARRHENIUS,.a=.05,.d=300}};
    const ChemMechanism m={4,1,0,sp,&reaction,reactants,products,NULL,NULL,298,101325};
    const ChemTerm b={1,1}, c={2,1};
    const ChemConstraint builtin={"eq",CHEM_EQUILIBRIUM_CONSTRAINT,2,&b,&c,1,1,0,5,0,0,0};
    const int deps[]={0,3}; Expression expr={.count=2,.species={0,3},.coefficients={1,1},.constant=-1,.parameter_term=-1};
    const ChemExternalConstraint ext={"mass",3,deps,2,evaluate_expression,&expr};
    const ChemJacobianEntry pat[]={{0,0},{3,0}}; Transfer transfer={.source=0,.product=3,.solvent=-1,.rate=.01};
    const ChemExternalModel model={.jacobian_pattern=pat,.jacobian_count=2,.constraints=&ext,.constraint_count=1,
        .finalize=transfer_finalize,.add_forcing=transfer_forcing,.subtract_jacobian=transfer_jacobian,.context=&transfer};
    ChemState *s=chem_state_create_extended(&m,&builtin,1,&model,1,1);
    if (!s) return 0;
    int ok=chem_is_algebraic(s,2)==1 && chem_is_algebraic(s,3)==1 && chem_find_constraint(s,"eq")==0 &&
           chem_find_constraint(s,"mass")==1 && chem_update_rates(s)==CHEM_OK;
    ChemStats stats; ChemRosenbrock p=chem_rosenbrock_four_stage_dae(); p.max_steps=10000;
    ok &= chem_solve(s,&p,10,&stats)==CHEM_OK;
    ok &= fabs(chem_species_values(s,0)[0]+chem_species_values(s,3)[0]-1)<1e-4 &&
          fabs(5*chem_species_values(s,1)[0]-chem_species_values(s,2)[0])<1e-4;
    chem_state_free(s); return ok;
}
static int builtin_pattern(void)
{
    Aerosol a; if (!aerosol_init(&a,1,0,1,0,.1)) return 0;
    int ok=chem_is_algebraic(a.state,1)==1 && a.transfer.finalized==1 && a.transfer.ps>=0;
    chem_species_values(a.state,0)[0]=.8; chem_species_values(a.state,1)[0]=.5; chem_species_values(a.state,2)[0]=.1;
    double f[3],j[9]; ok &= chem_update_rates(a.state)==CHEM_OK && chem_evaluate(a.state,0,f,j)==CHEM_OK;
    ok &= j[1*3+0]==0 && j[1*3+1]==1 && j[1*3+2]==-5;
    aerosol_free(&a); return ok;
}
static int finite_difference(int kind)
{
    Aerosol a; if (!aerosol_init(&a,2,kind!=0,0,0,kind==0 ? 3.5 : .1)) return 0;
    for (size_t c=0;c<2;++c) { chem_species_values(a.state,0)[c]=.8-.5*c; chem_species_values(a.state,1)[c]=.2+.5*c; }
    if (kind==2) { a.expression.coefficients[0]=2.5; a.expression.coefficients[1]=-1; a.expression.constant=0; }
    int ok=chem_update_rates(a.state)==CHEM_OK;
    for (size_t cell=0;cell<2;++cell) {
        double f[2],j[4]; ok &= chem_evaluate(a.state,cell,f,j)==CHEM_OK;
        for (size_t col=0;col<2;++col) {
            double *x=chem_species_values(a.state,col)+cell, old=*x,h=1e-6, plus[2],minus[2],unused[4];
            *x=old+h; ok &= chem_evaluate(a.state,cell,plus,unused)==CHEM_OK;
            *x=old-h; ok &= chem_evaluate(a.state,cell,minus,unused)==CHEM_OK; *x=old;
            for (size_t row=0;row<2;++row) ok &= fabs((plus[row]-minus[row])/(2*h)+j[row*2+col])<1e-8;
        }
    }
    aerosol_free(&a); return ok;
}

/* Kinetic and DAE networks use the same reactions as the reference tests. */
typedef struct {
    ChemSpecies species[4]; int reactants[5]; ChemProduct products[5]; ChemReaction reactions[5];
    ChemMechanism mechanism; int dependencies[3][4]; Expression expressions[3];
    ChemExternalConstraint constraints[3]; ChemExternalModel models[3];
    ChemTerm built_reactant,built_product; ChemConstraint built_constraint; ChemState *state;
} Network;
static int temperature_update(ChemState *s,void *context)
{
    (void)context; double *p=chem_parameter_values(s,0), *temperature=chem_temperatures(s);
    for (size_t c=0;c<chem_cell_count(s);++c) p[c]=2*exp(3000*(1/298.15-1/temperature[c]));
    return CHEM_OK;
}
static int network_init(Network *n,int mode,int noop_first,int temperature)
{
    /* mode: 0 kinetic, 1 simple ratio, 2 conservative, 3 two equilibria,
       4 built-in ratio, 5 kinetic two equilibria */
    static const char *names[]={"A","B","C","D"}, *param_names[]={"K_eq"};
    static const double defaults[]={0};
    memset(n,0,sizeof(*n));
    size_t count=(mode==3 || mode==5) ? 4 : 3;
    for (size_t i=0;i<count;++i) n->species[i]=(ChemSpecies){names[i],i==0 ? 1 : 0,1e-3,0,0};
    const int src[]={0,1,2,1,3}, dst[]={1,2,1,3,1}; const double rate[]={.1,10,2,8,8.0/3};
    size_t reactions=mode==0 ? 3 : mode==5 ? 5 : 1;
    for (size_t i=0;i<reactions;++i) {
        n->reactants[i]=src[i]; n->products[i]=(ChemProduct){dst[i],1};
        n->reactions[i]=(ChemReaction){i,1,i,1,0,{.kind=CHEM_ARRHENIUS,.a=rate[i],.d=300}};
    }
    n->mechanism=(ChemMechanism){count,reactions,temperature ? 1 : 0,n->species,n->reactions,n->reactants,n->products,
                                  temperature ? param_names : NULL,temperature ? defaults : NULL,298.15,101325};
    if (mode==4) {
        n->built_reactant=(ChemTerm){1,1}; n->built_product=(ChemTerm){2,1};
        n->built_constraint=(ChemConstraint){"BC",CHEM_EQUILIBRIUM_CONSTRAINT,2,&n->built_reactant,&n->built_product,1,1,0,5,0,0,0};
    }
    size_t num=0;
    if (mode==1 || mode==2 || mode==3) {
        n->dependencies[num][0]=1; n->dependencies[num][1]=2;
        n->expressions[num]=(Expression){.count=2,.species={1,2},.coefficients={5,-1},.parameter_term=temperature ? 0 : -1,.parameter=0};
        n->constraints[num]=(ChemExternalConstraint){"BC",2,n->dependencies[num],2,evaluate_expression,&n->expressions[num]};
        ++num;
    }
    if (mode==3) {
        n->dependencies[num][0]=1; n->dependencies[num][1]=3;
        n->expressions[num]=(Expression){.count=2,.species={1,3},.coefficients={3,-1},.parameter_term=-1};
        n->constraints[num]=(ChemExternalConstraint){"BD",3,n->dependencies[num],2,evaluate_expression,&n->expressions[num]};
        ++num;
    }
    if (mode==2 || mode==3) {
        for (size_t i=0;i<count;++i) n->dependencies[num][i]=(int)i;
        n->expressions[num]=(Expression){.count=count,.species={0,1,2,3},.coefficients={1,1,1,1},.constant=-1,.parameter_term=-1};
        n->constraints[num]=(ChemExternalConstraint){"mass",1,n->dependencies[num],count,evaluate_expression,&n->expressions[num]};
        ++num;
    }
    /* Separate model descriptors exercise composition and ordering. */
    if (noop_first) {
        if (num==3) return 0;
        for (size_t i=0;i<num;++i) n->models[i+1]=(ChemExternalModel){.constraints=&n->constraints[i],.constraint_count=1};
        n->models[0]=(ChemExternalModel){0};
        if (temperature) n->models[1].update=temperature_update;
        ++num;
    } else {
        for (size_t i=0;i<num;++i) n->models[i]=(ChemExternalModel){.constraints=&n->constraints[i],.constraint_count=1};
        if (temperature) n->models[0].update=temperature_update;
    }
    n->state=chem_state_create_extended(&n->mechanism,mode==4 ? &n->built_constraint : NULL,mode==4 ? 1 : 0,n->models,num,1);
    return n->state!=NULL;
}
static int advance(Network *n,int steps,double dt)
{
    ChemRosenbrock p=chem_rosenbrock_four_stage_dae();
    if (!chem_is_algebraic(n->state,1) && !chem_is_algebraic(n->state,2) &&
        (n->mechanism.species_count<4 || !chem_is_algebraic(n->state,3))) p=chem_rosenbrock_three_stage();
    for (int i=0;i<steps;++i) {
        ChemStats stats;
        if (chem_update_rates(n->state)!=CHEM_OK || chem_solve(n->state,&p,dt,&stats)!=CHEM_OK) return 0;
    }
    return 1;
}
static int convergence(int mode,double out[4])
{
    Network n; if (!network_init(&n,mode,0,0)) return 0;
    int ok=advance(&n,200,1);
    for (size_t i=0;i<n.mechanism.species_count;++i) out[i]=chem_species_values(n.state,i)[0];
    chem_state_free(n.state); return ok;
}
static int kinetic_conservative(int simple)
{
    double kinetic[4]={0},dae[4]={0}; int ok=convergence(0,kinetic) && convergence(simple ? 1 : 2,dae);
    if (simple) return ok && fabs(kinetic[2]/kinetic[1]-5)<1e-3 && fabs(dae[2]/dae[1]-5)<1e-3 &&
                        fabs(kinetic[0]+kinetic[1]+kinetic[2]-1)<1e-6 && fabs(dae[0]+dae[1]+dae[2]-6)<.1;
    for (size_t i=0;i<3;++i) ok &= fabs(kinetic[i]-dae[i])<1e-3;
    return ok && fabs(kinetic[0])<1e-6 && fabs(dae[1]-1.0/6)<1e-3 && fabs(dae[2]-5.0/6)<1e-3 &&
           fabs(kinetic[0]+kinetic[1]+kinetic[2]-1)<1e-6 && fabs(dae[0]+dae[1]+dae[2]-1)<1e-3;
}
static int multi_equilibrium(void)
{
    double kinetic[4]={0},dae[4]={0}; int ok=convergence(5,kinetic) && convergence(3,dae);
    const double expected[]={0,1.0/9,5.0/9,3.0/9}; double ks=0,ds=0;
    for (size_t i=0;i<4;++i) { ok &= fabs(kinetic[i]-expected[i])<1e-3 && fabs(dae[i]-expected[i])<1e-3; ks+=kinetic[i]; ds+=dae[i]; }
    return ok && fabs(ks-1)<1e-6 && fabs(ds-1)<1e-3 && fabs(5*dae[1]-dae[2])<1e-6 && fabs(3*dae[1]-dae[3])<1e-6;
}
static int compare_stepwise(void)
{
    Network built,external;
    if (!network_init(&built,4,0,0)) return 0;
    if (!network_init(&external,1,0,0)) { chem_state_free(built.state); return 0; }
    int ok=1;
    for (int step=0;step<100;++step) {
        ok &= advance(&built,1,.5) && advance(&external,1,.5);
        for (size_t i=0;i<3;++i) ok &= fabs(chem_species_values(built.state,i)[0]-chem_species_values(external.state,i)[0])<1e-8;
        ok &= fabs(5*chem_species_values(external.state,1)[0]-chem_species_values(external.state,2)[0])<1e-6;
    }
    chem_state_free(built.state); chem_state_free(external.state); return ok;
}
static int ordering(void)
{
    Network first,second;
    if (!network_init(&first,1,0,0)) return 0;
    if (!network_init(&second,1,1,0)) { chem_state_free(first.state); return 0; }
    int ok=advance(&first,200,1) && advance(&second,200,1);
    for (size_t i=0;i<3;++i) ok &= chem_species_values(first.state,i)[0]==chem_species_values(second.state,i)[0];
    ok &= fabs(5*chem_species_values(second.state,1)[0]-chem_species_values(second.state,2)[0])<1e-6;
    chem_state_free(first.state); chem_state_free(second.state); return ok;
}
static int build_multiple(void)
{
    Network n; if (!network_init(&n,1,0,0)) return 0;
    chem_state_free(n.state);
    n.species[3]=(ChemSpecies){"D",2,1e-3,0,0}; n.mechanism.species_count=4;
    ChemTerm a={0,1},d={3,1};
    ChemConstraint builtin={"AD",CHEM_EQUILIBRIUM_CONSTRAINT,3,&a,&d,1,1,0,2,0,0,0};
    n.state=chem_state_create_extended(&n.mechanism,&builtin,1,n.models,1,1);
    if (!n.state) return 0;
    ChemState *other=chem_state_create_extended(&n.mechanism,&builtin,1,n.models,1,1);
    if (!other) { chem_state_free(n.state); return 0; }
    ChemState *original=n.state; int ok=advance(&n,20,1); n.state=other; ok &= advance(&n,20,1);
    for (size_t i=0;i<4;++i) ok &= chem_species_values(original,i)[0]==chem_species_values(other,i)[0];
    ok &= fabs(2*chem_species_values(original,0)[0]-chem_species_values(original,3)[0])<1e-6 &&
          fabs(5*chem_species_values(original,1)[0]-chem_species_values(original,2)[0])<1e-6;
    chem_state_free(original); chem_state_free(other); return ok;
}
static int temperature_constraint(void)
{
    int ok=1;
    for (int i=0;i<2;++i) {
        Network n; if (!network_init(&n,1,0,1)) return 0;
        double temperature=i ? 350 : 298.15, expected=2*exp(3000*(1/298.15-1/temperature));
        chem_temperatures(n.state)[0]=temperature;
        ok &= advance(&n,200,1);
        double b=chem_species_values(n.state,1)[0], c=chem_species_values(n.state,2)[0];
        ok &= b>0 && fabs(c/b-expected)<1e-4;
        chem_state_free(n.state);
    }
    return ok;
}
static int validate_and_clone(void)
{
    Aerosol a; if (!aerosol_init(&a,2,1,0,0,.1)) return 0;
    chem_species_values(a.state,0)[0]=.8; chem_species_values(a.state,0)[1]=.6;
    chem_species_values(a.state,1)[0]=.2; chem_species_values(a.state,1)[1]=.4;
    ChemState *clone=chem_state_clone(a.state), *other=chem_state_create_like(a.state,3);
    int ok=clone && other && a.transfer.finalized==1;
    if (clone && other) {
        ChemExternalModel duplicate[]={a.model,a.model};
        ChemState *invalid=chem_state_create_extended(&a.mechanism,NULL,0,duplicate,2,1);
        ok &= invalid==NULL; chem_state_free(invalid);
        ChemJacobianEntry invalid_entry={99,0}; ChemExternalModel bad=a.model;
        bad.jacobian_pattern=&invalid_entry; bad.jacobian_count=1;
        invalid=chem_state_create_extended(&a.mechanism,NULL,0,&bad,1,1);
        ok &= invalid==NULL; chem_state_free(invalid);
        ChemExternalModel duplicate_species={.species=a.species,.species_count=1};
        ChemMechanismComposition composition={0};
        ok &= chem_mechanism_compose(&composition,&a.mechanism,&duplicate_species,1)==CHEM_INVALID;
        chem_mechanism_composition_destroy(&composition);
        aerosol_free(&a); a.state=NULL;
        ok &= chem_species_values(clone,0)[0]==.8 && chem_species_values(clone,0)[1]==.6 &&
              chem_cell_count(other)==3 && chem_update_rates(clone)==CHEM_OK;
        ChemStats stats; ChemRosenbrock p=chem_rosenbrock_four_stage_dae();
        ok &= chem_solve(clone,&p,2,&stats)==CHEM_OK &&
              fabs(chem_species_values(clone,0)[0]+chem_species_values(clone,1)[0]-1)<1e-4;
    }
    chem_state_free(clone); chem_state_free(other); aerosol_free(&a);
    return ok;
}
static int failures;
static void report(const char *group,const char *name,int ok)
{
    printf("external_model_constraints:%s.%s %s %s\n",group,name,ok ? "PASS" : "FAIL",
           ok ? "reference extension behavior matched" : "extension behavior differed"); failures+=!ok;
}
int main(void)
{
    const char *g="ExternalModelConstraints", *fd="ExternalModelFiniteDifferenceJacobian";
    if (!validate_and_clone()) { fprintf(stderr,"extension validation or state cloning failed\n"); ++failures; }
    report(g,"AddExternalModelWithConstraints",aerosol_metadata(1));
    report(g,"AddExternalModelProcessOnly",aerosol_metadata(0));
    report(g,"DAESolveEnforcesConservation",aerosol_conservation(1,0,20,10));
    report(g,"AlgebraicRowHoldsOnlyConstraintTerms",algebraic_row());
    report(g,"CombinedBuiltInAndExternalConstraints",combined_constraints());
    report(g,"AddExternalModelOnlyStandardRosenbrock",aerosol_metadata(0));
    report(g,"AddExternalModelConstraintsOnly",aerosol_conservation(1,0,1,50));
    report(g,"MultiGridCell",aerosol_conservation(3,0,1,50));
    report(g,"KineticVsConservativeConstraintConvergence",kinetic_conservative(0));
    report(g,"SimpleConstraintPreservesRatioNotMass",kinetic_conservative(1));
    report(g,"BuiltInVsExternalModelConstraintStepByStep",compare_stepwise());
    report(g,"MultiEquilibriumKineticVsComposedConstraints",multi_equilibrium());
    report(g,"ProcessJacobianElementInAlgebraicRowSurvivesFiltering",aerosol_conservation(1,1,10,10));
    report(g,"ProcessJacobianElementInBuiltInAlgebraicRowSurvivesFiltering",builtin_pattern());
    report(g,"TemperatureDependentConstraintParameter",temperature_constraint());
    report(g,"ConstraintActiveWhenProcessOnlyModelAddedFirst",ordering());
    report(g,"BuilderCanBuildMultipleSolvers",build_multiple());
    report(fd,"ProcessForcingJacobian",finite_difference(0));
    report(fd,"ConstraintResidualJacobian",finite_difference(1));
    report(fd,"EquilibriumConstraintModelJacobian",finite_difference(2));
    return failures ? 1 : 0;
}
