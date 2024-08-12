#include "../../../include/crystalPlasticity.h"
//constructor
template <int dim>
crystalPlasticity<dim>::crystalPlasticity(userInputParameters & _userInputs_cp):
multiPhysicsBVP<dim>(_userInputs_cp),
F(dim,dim),
F_tau(dim,dim),
FP_tau(dim,dim),
FE_tau(dim,dim),
T(dim,dim),
P(dim,dim)
{
    initCalled = false;

    //post processing
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("Eqv_stress");
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("Eqv_strain");
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("Twin");
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("output_Var1");
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("output_Var2");
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("output_Var3");
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("output_Var4");
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("output_Var5");
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("output_Var6");
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("output_Var7");
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("output_Var8");
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("output_Var9");
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("output_Var10");
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("output_Var11");
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("output_Var12");
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("output_Var13");
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("output_Var14");
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("output_Var15");
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("output_Var16");
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("output_Var17");
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("output_Var18");
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("output_Var19");
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("output_Var20");
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("output_Var21");
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("output_Var22");
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("output_Var23");
    multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("output_Var24");
    multiPhysicsBVP<dim>::numPostProcessedFields=27;

    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("meshGrain_ID");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var1");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var2");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var3");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var4");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var5");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var6");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var7");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var8");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var9");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var10");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var11");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var12");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var13");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var14");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var15");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var16");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var17");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var18");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var19");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var20");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var21");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var22");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var23");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var24");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var25");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var26");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var27");
    multiPhysicsBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var28");
    multiPhysicsBVP<dim>::numPostProcessedFieldsAtCellCenters=29; //grainID

}

#include "../../../include/crystalPlasticity_template_instantiations.h"
