#include "customPDE.h"

// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
// This function sets attributes for each variable/equation in the app. The
// attributes are set via standardized function calls. The first parameter for each
// function call is the variable index (starting at zero). The first set of
// variable/equation attributes are the variable name (any string), the variable
// type (SCALAR/VECTOR), and the equation type (EXPLICIT_TIME_DEPENDENT/
// TIME_INDEPENDENT/AUXILIARY). The next set of attributes describe the
// dependencies for the governing equation on the values and derivatives of the
// other variables for the value term and gradient term of the RHS and the LHS.
// The final pair of attributes determine whether a variable represents a field
// that can nucleate and whether the value of the field is needed for nucleation
// rate calculations.

void variableAttributeLoader::loadVariableAttributes(){
	// Variable 0 - Order Parameter
	set_variable_name				(0,"n");
	set_variable_type				(0,SCALAR);
	set_variable_equation_type		(0,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(0, "n, dndt, strain_df");
    set_dependencies_gradient_term_RHS(0, "grad(n)");
	
	// Variable 1 - Time Derivative of Order Parameter
	set_variable_name				(1,"dndt");
	set_variable_type				(1,SCALAR);
	set_variable_equation_type		(1,AUXILIARY);

    set_dependencies_value_term_RHS(1, "n, dndt, strain_df");
    set_dependencies_gradient_term_RHS(1, "grad(n)");

	// Variable 2 - Strain driving force
	set_variable_name				(2,"strain_df");
	set_variable_type				(2,SCALAR);
	set_variable_equation_type		(2,AUXILIARY);

    set_dependencies_value_term_RHS(2, "strain_df");
    set_dependencies_gradient_term_RHS(2, "");

}

// =============================================================================================
// explicitEquationRHS (needed only if one or more equation is explict time dependent)
// =============================================================================================
// This function calculates the right-hand-side of the explicit time-dependent
// equations for each variable. It takes "variable_list" as an input, which is a list
// of the value and derivatives of each of the variables at a specific quadrature
// point. The (x,y,z) location of that quadrature point is given by "q_point_loc".
// The function outputs two terms to variable_list -- one proportional to the test
// function and one proportional to the gradient of the test function. The index for
// each variable in this list corresponds to the index given at the top of this file.

template <int dim, int degree>
void customPDE<dim,degree>::explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double>> & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc) const {

	// --- Getting the values and derivatives of the model variables ---

	// The order parameter and its derivatives
	scalarvalueType_pf n = variable_list.get_scalar_value(0);
	scalargradType_pf nx = variable_list.get_scalar_gradient(0);

	// The time derivative of the order parameter
	scalarvalueType_pf dndt = variable_list.get_scalar_value(1);

	// The strain contribiution to the driving force
	scalarvalueType_pf strain_df= variable_list.get_scalar_value(2);

	// --- Setting the expressions for the terms in the governing equations ---

	scalarvalueType_pf mu_twV = constV(delf_tw)*(4.0*n*(n-1.0)*(n-0.5));
	scalargradType_pf kappagradn;
	kappagradn[0] = constV(K[0][0])*nx[0]+constV(K[0][1])*nx[1]+constV(K[0][2])*nx[2];
	kappagradn[1] = constV(K[1][0])*nx[0]+constV(K[1][1])*nx[1]+constV(K[1][2])*nx[2];
	kappagradn[2] = constV(K[2][0])*nx[0]+constV(K[2][1])*nx[1]+constV(K[2][2])*nx[2];

	//Outward Normal vector
	scalargradType_pf nvec = -nx/(std::sqrt(nx[0]*nx[0] + nx[1]*nx[1] + nx[2]*nx[2])+constV(regval));

	//Computing the outward mobility (L = grad(nvec) dot Ltens dot grad(nvec))
	scalarvalueType_pf L = constV(0.0);
	for(unsigned int i=0;i<dim;i++){
		for(unsigned int j=0;j<dim;j++){
			//Mobility tensor (rotated)
			L = L + nvec[i]*nvec[j]*Ltens[i][j];
		}
}

//Applying a filter to localize driving force to the twin boundary 
scalarvalueType_pf strain_df_filter = (1.0 - (2.0*n-1.0)*(2.0*n-1.0))*strain_df;

//Defining the value and gradient terms
scalarvalueType_pf eq_n = (n-constV(userInputs_pf.dtValue)*L*(mu_twV-strain_df_filter));
scalargradType_pf eqx_n = -(constV(userInputs_pf.dtValue)*L*kappagradn);

// --- Submitting the terms for the governing equations ---

variable_list.set_scalar_value_term_RHS(0,eq_n);
variable_list.set_scalar_gradient_term_RHS(0,eqx_n);

}

// =============================================================================================
// nonExplicitEquationRHS (needed only if one or more equation is time independent or auxiliary)
// =============================================================================================
// This function calculates the right-hand-side of all of the equations that are not
// explicit time-dependent equations. It takes "variable_list" as an input, which is
// a list of the value and derivatives of each of the variables at a specific
// quadrature point. The (x,y,z) location of that quadrature point is given by
// "q_point_loc". The function outputs two terms to variable_list -- one proportional
// to the test function and one proportional to the gradient of the test function. The
// index for each variable in this list corresponds to the index given at the top of
// this file.

template <int dim, int degree>
void customPDE<dim,degree>::nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
// --- Getting the values and derivatives of the model variables ---

// The order parameter and its derivatives
scalarvalueType_pf n = variable_list.get_scalar_value(0);
scalargradType_pf nx = variable_list.get_scalar_gradient(0);

// The time derivative of the order parameter
scalarvalueType_pf dndt = variable_list.get_scalar_value(1);

// The strain contribiution to the driving force
scalarvalueType_pf strain_df= variable_list.get_scalar_value(2);

// --- Setting the expressions for the terms in the governing equations ---

scalarvalueType_pf mu_twV = constV(delf_tw)*(4.0*n*(n-1.0)*(n-0.5));
scalargradType_pf kappagradn;
kappagradn[0] = constV(K[0][0])*nx[0]+constV(K[0][1])*nx[1]+constV(K[0][2])*nx[2];
kappagradn[1] = constV(K[1][0])*nx[0]+constV(K[1][1])*nx[1]+constV(K[1][2])*nx[2];
kappagradn[2] = constV(K[2][0])*nx[0]+constV(K[2][1])*nx[1]+constV(K[2][2])*nx[2];

//Outward Normal vector
scalargradType_pf nvec = -nx/(std::sqrt(nx[0]*nx[0] + nx[1]*nx[1] + nx[2]*nx[2])+constV(regval));

//Computing the outward mobility (L = grad(nvec) dot Ltens dot grad(nvec))
scalarvalueType_pf L = constV(0.0);
for(unsigned int i=0;i<dim;i++){
	for(unsigned int j=0;j<dim;j++){
		//Mobility tensor (rotated)
		L = L + nvec[i]*nvec[j]*Ltens[i][j];;
	}
}

//Applying a filter to localize driving force to the twin boundary 
scalarvalueType_pf strain_df_filter = (1.0 - (2.0*n-1.0)*(2.0*n-1.0))*strain_df;

scalarvalueType_pf eq_dndt = -L*(mu_twV-strain_df_filter);
scalargradType_pf eqx_dndt = -L*kappagradn;

// --- Submitting the terms for the governing equations ---

variable_list.set_scalar_value_term_RHS(1,eq_dndt);
variable_list.set_scalar_gradient_term_RHS(1,eqx_dndt);

variable_list.set_scalar_value_term_RHS(2,strain_df);

}

// =============================================================================================
// equationLHS (needed only if at least one equation is time independent)
// =============================================================================================
// This function calculates the left-hand-side of time-independent equations. It
// takes "variable_list" as an input, which is a list of the value and derivatives of
// each of the variables at a specific quadrature point. The (x,y,z) location of that
// quadrature point is given by "q_point_loc". The function outputs two terms to
// variable_list -- one proportional to the test function and one proportional to the
// gradient of the test function -- for the left-hand-side of the equation. The index
// for each variable in this list corresponds to the index given at the top of this
// file. If there are multiple elliptic equations, conditional statements should be
// sed to ensure that the correct residual is being submitted. The index of the field
// being solved can be accessed by "this->currentFieldIndex".

template <int dim, int degree>
void customPDE<dim,degree>::equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
}