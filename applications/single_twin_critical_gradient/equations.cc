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
    set_variable_name                (0,"n");
    set_variable_type                (0,SCALAR);
    set_variable_equation_type        (0,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(0, "n, grad(n), dndt");
    set_dependencies_gradient_term_RHS(0, "");
    
    // Variable 1 - Time Derivative of Order Parameter
    set_variable_name                (1,"dndt");
    set_variable_type                (1,SCALAR);
    set_variable_equation_type        (1,AUXILIARY);

    set_dependencies_value_term_RHS(1, "n, grad(n), strain_df");
    set_dependencies_gradient_term_RHS(1, "grad(n)");

    // Variable 2 - Strain driving force
    set_variable_name                (2,"strain_df");
    set_variable_type                (2,SCALAR);
    set_variable_equation_type        (2,AUXILIARY);

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
void customPDE<dim,degree>::explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double>> &variable_list,
        dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc) const {

    // The time derivative of the order parameter is calculated as an AUXILIARY field,
    // in nonExplicitEquationRHS(...) below. Then, here we use that value to set the
    // change in the order parameter, including capping the order parameter between
    // zero and one.

    // Get the value of the order parameter
    scalarvalueType_pf n = variable_list.get_scalar_value(0);
    scalargradType_pf nx = variable_list.get_scalar_gradient(0);

    // Get the time derivative of the order parameter (calculated as AUXILIARY field)
    scalarvalueType_pf dndt = variable_list.get_scalar_value(1);

    // Prevent the order parameter from decreasing (no detwinning)
    dndt = std::max(dndt, constV(0.0));

    // Only allow twinning to initiate if the critical gradient is exceeded
    for (unsigned int v = 0; v < nx[0].size(); v++)
      {
        if (n[v] <= 0.0 && nx.norm()[v] < critical_grad)
          {
            dndt[v] = 0.0;
          }
      }

    // Calculate the new order parameter value
    scalarvalueType_pf new_n = n + constV(userInputs_pf.dtValue) * dndt;

    // Restrict the new order parameter to be between zero and one
    new_n = std::min(new_n, constV(1.0));
    new_n = std::max(new_n, constV(0.0));

    // --- Submitting the terms for the governing equations ---
    variable_list.set_scalar_value_term_RHS(0,new_n);
}

// Rotation conversion function, copied from crystalPlasticity<dim>::odfpoint()
template <int dim, int degree>
void customPDE<dim, degree>::rodrigues_to_rotmat(dealii::Tensor<2,dim> &orientationMatrix,
        dealii::Tensor<1,dim> r) const {    
    double rdotr = 0.0;

    for (unsigned int i = 0; i < dim; i++)
      {
        rdotr = rdotr + r[i]*r[i];
      }
    
    double term1 = 1.0 - rdotr;
    double term2 = 1.0 + rdotr;

    orientationMatrix.clear();
    orientationMatrix[0][0] = 1.0;
    orientationMatrix[1][1] = 1.0;
    orientationMatrix[2][2] = 1.0;

    for (unsigned int i = 0; i < dim; i++)
      {
        orientationMatrix[i][i] = orientationMatrix[i][i]*term1;
      }
    
    for (unsigned int i = 0; i < dim; i++)
      {
        for (unsigned int j = 0; j < dim; j++)
          {
            orientationMatrix[i][j] = orientationMatrix[i][j] + 2.0*r[i]*r[j];
          }
      }
    
    orientationMatrix[0][1] = orientationMatrix[0][1] - 2.0*r[2];
    orientationMatrix[0][2] = orientationMatrix[0][2] + 2.0*r[1];
    orientationMatrix[1][2] = orientationMatrix[1][2] - 2.0*r[0];
    orientationMatrix[1][0] = orientationMatrix[1][0] + 2.0*r[2];
    orientationMatrix[2][0] = orientationMatrix[2][0] - 2.0*r[1];
    orientationMatrix[2][1] = orientationMatrix[2][1] + 2.0*r[0];

    for(unsigned int i=0;i<dim;i++)
      {
        for(unsigned int j=0;j<dim;j++)
          {
            orientationMatrix[i][j] = orientationMatrix[i][j]/term2;
          }
      }
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
void customPDE<dim,degree>::nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double>> &variable_list,
        dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc) const {
    // --- Getting the values and derivatives of the model variables ---

    // The order parameter and its derivatives
    scalarvalueType_pf n = variable_list.get_scalar_value(0);
    scalargradType_pf nx = variable_list.get_scalar_gradient(0);

    // The strain contribiution to the driving force
    scalarvalueType_pf strain_df= variable_list.get_scalar_value(2);

    // --- Setting the expressions for the terms in the governing equations ---

    scalarvalueType_pf mu_twV = constV(delf_tw)*(4.0*n*(n-1.0)*(n-0.5));
    scalargradType_pf kappagradn;
    scalarvalueType_pf L = constV(0.0);

    //Outward Normal vector
    scalargradType_pf nvec = -nx/(std::sqrt(nx[0]*nx[0] + nx[1]*nx[1] + nx[2]*nx[2])+constV(regval));

    // q_point_loc is vectorized, but CPFE functions are not.
    // TODO: created vectorized versions of the relevant CPFE functions to fix this

    // For now, unroll the vectorization
    for (unsigned int v = 0; v < q_point_loc[0].size(); v++)
      {
        // Get the materialID (a.k.a. grain ID) from CPFE for this point
        double coords[3] = { q_point_loc[0][v], q_point_loc[1][v], q_point_loc[2][v] };
        unsigned int materialID = this->cpfe_orientations->getMaterialID(coords);

        dealii::Tensor<2,dim> K = Kij_map[materialID];
        dealii::Tensor<2,dim> Ltens = Lij_map[materialID];

        kappagradn[0][v] = K[0][0]*nx[0][v] + K[0][1]*nx[1][v] + K[0][2]*nx[2][v];
        kappagradn[1][v] = K[1][0]*nx[0][v] + K[1][1]*nx[1][v] + K[1][2]*nx[2][v];
        kappagradn[2][v] = K[2][0]*nx[0][v] + K[2][1]*nx[1][v] + K[2][2]*nx[2][v];

        for (unsigned int i = 0; i < dim; i++)
          {
            for (unsigned int j = 0; j < dim; j++)
              {
                L[v] += nvec[i][v] * nvec[j][v] * Ltens[i][j];
              }
          }
      }
    
    L = std::max(L, constV(minL));

    //Applying a filter to localize driving force to the twin boundary 
    scalarvalueType_pf strain_df_filter = 1.5*(1.0 - (2.0*n-1.0)*(2.0*n-1.0))*strain_df;

    scalarvalueType_pf eq_dndt = -L*(mu_twV-strain_df_filter);
    scalargradType_pf eqx_dndt = -L*(kappagradn);

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
