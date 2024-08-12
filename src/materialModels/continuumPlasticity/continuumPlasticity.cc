#include "../../../include/continuumPlasticity.h"
//constructor
template <int dim>
continuumPlasticity<dim>::continuumPlasticity(userInputParameters_cp & _userInputs_cp)
:
multiPhysicsBVP<dim>(_userInputs_cp),
  enhStrain(this->FE,this->pcout, _userInputs_cp)
{
  //initialize "initCalled"
  initCalled = false;
  //initCalled = false;

  //post processing
  multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("Eqv_stress");
  multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("Eqv_strain");
  multiPhysicsBVP<dim>::postprocessed_solution_names.push_back("alpha");
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
  multiPhysicsBVP<dim>::numPostProcessedFieldsAtCellCenters=1; //grainID
}

#include "../../../include/continuumPlasticity_template_instantiations.h"
