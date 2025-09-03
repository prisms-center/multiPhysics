#include "../../../include/crystalPlasticity.h"

#include "../../../include/matrixFreePDE.h"

// constructor
template <int dim>
crystalPlasticity<dim>::crystalPlasticity(userInputParameters_pf<dim> _userInputs_pf,
                                          userInputParameters_cp     &_userInputs_cp,
                                          MatrixFreePDE<dim, 1>      &_pf_object)
  : MultiPhysicsBVP<dim, 1>(_userInputs_pf, _userInputs_cp)
  , F(dim, dim)
  , F_tau(dim, dim)
  , FP_tau(dim, dim)
  , FE_tau(dim, dim)
  , T(dim, dim)
  , P(dim, dim)
  , pf_object(_pf_object)
{
  initCalled = false;

  // post processing
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("Eqv_stress");
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("Eqv_strain");
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("Twin");
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("output_Var1");
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("output_Var2");
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("output_Var3");
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("output_Var4");
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("output_Var5");
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("output_Var6");
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("output_Var7");
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("output_Var8");
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("output_Var9");
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("output_Var10");
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("output_Var11");
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("output_Var12");
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("output_Var13");
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("output_Var14");
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("output_Var15");
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("output_Var16");
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("output_Var17");
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("output_Var18");
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("output_Var19");
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("output_Var20");
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("output_Var21");
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("output_Var22");
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("output_Var23");
  MultiPhysicsBVP<dim, 1>::postprocessed_solution_names.push_back("output_Var24");
  MultiPhysicsBVP<dim, 1>::numPostProcessedFields = 27;

  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "meshGrain_ID");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var1");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var2");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var3");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var4");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var5");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var6");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var7");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var8");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var9");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var10");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var11");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var12");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var13");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var14");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var15");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var16");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var17");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var18");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var19");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var20");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var21");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var22");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var23");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var24");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var25");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var26");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var27");
  MultiPhysicsBVP<dim, 1>::postprocessedFieldsAtCellCenters_solution_names.push_back(
    "outputCellCenters_Var28");
  MultiPhysicsBVP<dim, 1>::numPostProcessedFieldsAtCellCenters = 29; // grainID
}

// PF interfacing methods
template <int dim>
std::vector<const DoFHandler<dim> *> &
crystalPlasticity<dim>::getDofHandlersSet()
{
  return pf_object.getDofHandlersSet();
}

template <int dim>
std::vector<vectorType_pf *> &
crystalPlasticity<dim>::getSolutionSet()
{
  return pf_object.getSolutionSet();
}

template <int dim>
std::vector<const AffineConstraints<double> *> &
crystalPlasticity<dim>::getConstraintsDirichletSet()
{
  return pf_object.getConstraintsDirichletSet();
}

template <int dim>
std::vector<const AffineConstraints<double> *> &
crystalPlasticity<dim>::getConstraintsOtherSet()
{
  return pf_object.getConstraintsOtherSet();
}

template <int dim>
void
crystalPlasticity<dim>::getOutputResults()
{
  return pf_object.getOutputResults();
}

template <int dim>
double &
crystalPlasticity<dim>::getCurrentTime_pf()
{
  return pf_object.getCurrentTime();
}

template <int dim>
unsigned int &
crystalPlasticity<dim>::getCurrentIncrement_pf()
{
  return pf_object.getCurrentIncrement();
}

template <int dim>
unsigned int &
crystalPlasticity<dim>::getCurrentOutput()
{
  return pf_object.getCurrentOutput();
}

#include "../../../include/crystalPlasticity_template_instantiations.h"
