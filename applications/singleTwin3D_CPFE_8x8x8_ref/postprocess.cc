#include "customPDE.h"
// =================================================================================
// Set the attributes of the postprocessing variables
// =================================================================================
// This function is analogous to 'loadVariableAttributes' in 'equations.cc', but for
// the postprocessing expressions. It sets the attributes for each postprocessing
// expression, including its name, whether it is a vector or scalar (only scalars are
// supported at present), its dependencies on other variables and their derivatives,
// and whether to calculate an integral of the postprocessed quantity over the entire
// domain.

void variableAttributeLoader::loadPostProcessorVariableAttributes(){
	// Variable 0
	set_variable_name				(0,"mu_twV");
	set_variable_type				(0,SCALAR);

	set_dependencies_value_term_RHS(0, "n, grad(n)");
  set_dependencies_gradient_term_RHS(0, "");

  set_output_integral         	(0,true);

}

// =============================================================================================
// postProcessedFields: Set the postprocessing expressions
// =============================================================================================
// This function is analogous to 'explicitEquationRHS' and 'nonExplicitEquationRHS' in
// equations.cc. It takes in "variable_list" and "q_point_loc" as inputs and outputs two terms in
// the expression for the postprocessing variable -- one proportional to the test
// function and one proportional to the gradient of the test function. The index for
// each variable in this list corresponds to the index given at the top of this file (for
// submitting the terms) and the index in 'equations.cc' for assigning the values/derivatives of
// the primary variables.

template <int dim,int degree>
void customPDE<dim,degree>::postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
												const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// --- Getting the values and derivatives of the model variables ---

// The order parameter and its derivatives
scalarvalueType_pf n = variable_list.get_scalar_value(0);
scalargradType_pf nx = variable_list.get_scalar_gradient(0);

scalarvalueType_pf mu_twV = constV(delf_tw)*(4.0*n*(n-1.0)*(n-0.5));

// --- Setting the expressions for the terms in the postprocessing expressions ---
pp_variable_list.set_scalar_value_term_RHS(0, mu_twV);

}