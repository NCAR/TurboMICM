/* SPDX-License-Identifier: Apache-2.0
 * Hand-built tutorial mechanisms, independent of the configuration generator. */
static ChemStatus hand_rates(int user,ChemMechanism **out)
{
    ChemSpecies sp[7]={{.name="A"},{.name="B"},{.name="C",.molecular_weight=.025,.diffusion=2.3e2},
                       {.name="D"},{.name="E"},{.name="F"},{.name="G"}};
    size_t lhs[10][2]={{0},{1},{1},{2},{3},{4,4},{5},{2},{0},{1}};
    ChemProduct rhs[10][1]={{{1,1}},{{2,1}},{{3,1}},{{4,1}},{{5,2}},{{6,1}},{{6,1}},{{6,1}},{{0,1}},{{0,0}}};
    ChemRateType types[]={CHEM_ARRHENIUS,CHEM_BRANCHED,CHEM_BRANCHED,CHEM_SURFACE,CHEM_TERNARY,CHEM_TROE,CHEM_TUNNELING,CHEM_USER,CHEM_USER,CHEM_USER};
    ChemReaction r[10]={0};
    for (size_t i=0;i<10;++i) {
        r[i].reactant_count=i==8?0:i==5?2:1;r[i].reactants=lhs[i];
        r[i].product_count=i==9?0:1;r[i].products=rhs[i];r[i].rate=chem_rate_default(types[i]);
    }
    r[0].rate.data.arrhenius=(ChemArrhenius){2.15e-4,0,110,300,0};
    r[1].rate.data.branched=(ChemBranched){.X=1.2,.Y=204.3,.a0=1e-3,.n=2,.branch=0};
    r[2].rate.data.branched=r[1].rate.data.branched;r[2].rate.data.branched.branch=1;
    r[3].rate.label="C surface";r[3].rate.data.surface=(ChemSurface){.diffusion=2.3e2,.molecular_weight=.025,.probability=.9};
    double conversion=1e-6*CHEM_AVOGADRO;
    r[4].rate.data.falloff=(ChemFalloff){1.2,2.3,302.3,2.6/conversion,-3.1,402.1,.9,1.2};
    r[5].rate.data.falloff=(ChemFalloff){1.8e-30*conversion*conversion,-3,0,2.8e-11*conversion,0,0,.6,1};
    r[6].rate.data.tunneling=(ChemTunneling){1.2,2.3,302.3};
    r[7].rate.label="my photolysis rate";r[8].rate.label="my emission rate";r[9].rate.label="my loss rate";
    return chem_mechanism_create(7,sp,user?10:7,r,out);
}
static ChemStatus hand_robertson(ChemMechanism **out)
{
    ChemSpecies sp[]={{.name="A"},{.name="B"},{.name="C"}};
    const size_t lhs[][2]={{0},{1,1},{1,2}};
    const ChemProduct rhs[][2]={{{1,1}},{{1,1},{2,1}},{{0,1},{2,1}}};
    const char *labels[]={"r1","r2","r3"};
    ChemReaction r[3]={0};
    for (size_t i=0;i<3;++i) {
        r[i].reactants=lhs[i];r[i].reactant_count=i?2:1;
        r[i].products=rhs[i];r[i].product_count=i?2:1;
        r[i].rate=chem_rate_default(CHEM_USER);r[i].rate.label=labels[i];
    }
    return chem_mechanism_create(3,sp,3,r,out);
}
