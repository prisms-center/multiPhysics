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

  // During seeding - recreate the initial condition, multiplied by a scalar factor
	if (this->seedingIncrement < this->userInputs_cp.stepsForSeeding)
	  {
		  double seeding_scale_factor = (this->seedingIncrement / ((double)this->userInputs_cp.stepsForSeeding));
			scalarType_pf seed_n = constV(0.0);
			scalargradType_pf seed_nx = variable_list.get_scalar_gradient(0);
			// unroll vectorization
			for (unsigned int v = 0; v < q_point_loc[0].size(); v++)
			  {
					// TODO: fix this.
					// Duplicate the entire ICs_and_BCs function because calling it in this context is disallowed by const-correctness.
					double scalar_IC;
					dealii::Point<dim, double> p;
					p[0] = q_point_loc[0][v];
					p[1] = q_point_loc[1][v];
					p[2] = q_point_loc[2][v];
					double center[3] = {0.5,0.5,0.0}; // TODO: change the center to reflect the actual domain size
					dealii::Tensor<1, dim> nc;
					dealii::Tensor<1, dim> n;
					dealii::Tensor<1, dim> n_int;
					dealii::Tensor<1, dim> td_int;
					dealii::Tensor<1, dim> tn_int;
					nc.clear(); n.clear(); n_int.clear(); td_int.clear(); tn_int.clear();
					double dist, edist;
					double b0=a0*std::sqrt((1.0-ecc*ecc));
					double nXc, nYc, nZc, nX, nY, nZ, nZ_reg;
					double pi = 3.14159265358979323846;
					scalar_IC = 0;

					//Calculating distance from center of the system
					dist = 0.0;
					for (unsigned int dir = 0; dir < dim; dir++){
						dist += (p[dir]-center[dir]*userInputs_pf.domain_size[dir])*(p[dir]-center[dir]*userInputs_pf.domain_size[dir]);
					}
					dist = std::sqrt(dist);
						
					//Elliptical seed
					//Calculating distance from center to perimeter of the ellipse
					//Components of the normal vector from the center to point p
					nc[0] = (p[0]-center[0]*userInputs_pf.domain_size[0])/(dist + 1.0e-7);
					nc[1] = (p[1]-center[1]*userInputs_pf.domain_size[1])/(dist + 1.0e-7);
					nc[2] = (p[2]-center[2]*userInputs_pf.domain_size[2])/(dist + 1.0e-7);

					//Rotating td and tn according to Euler angles
					// Euler angles in radians (ZXZ convention)
					double phi1 = pi*euler_angs[0]/180.0; // e.g., 0.785398 for 45 degrees
					double Phi  = pi*euler_angs[1]/180.0; // e.g., 1.0472   for 60 degrees
					double phi2 = pi*euler_angs[2]/180.0; // e.g., 0.523599 for 30 degrees

					// Step 1: Build the rotation matrices Rz(phi1), Rx(Phi), Rz(phi2)
					dealii::Tensor<2, dim> Rz1, Rx, Rz2;
					Rz1.clear(); Rx.clear(); Rz2.clear();

					// Initialize identity diagonals first
					for (unsigned int i = 0; i < dim; ++i)
					  {
							Rz1[i][i] = 1.0;
							Rx[i][i] = 1.0;
							Rz2[i][i] = 1.0;
					  }

					// Rz(phi1)
					Rz1[0][0] = std::cos(phi1); Rz1[0][1] = -std::sin(phi1);
					Rz1[1][0] = std::sin(phi1); Rz1[1][1] =  std::cos(phi1);

					// Rx(Phi)
					Rx[1][1] = std::cos(Phi); Rx[1][2] = -std::sin(Phi);
					Rx[2][1] = std::sin(Phi); Rx[2][2] =  std::cos(Phi);

					// Rz(phi2)
					Rz2[0][0] = std::cos(phi2); Rz2[0][1] = -std::sin(phi2);
					Rz2[1][0] = std::sin(phi2); Rz2[1][1] =  std::cos(phi2);

					// Step 2: Compute Q = Rz1 * Rx * Rz2
					dealii::Tensor<2, dim> Q_temp, Q_1;
					Q_temp.clear(); Q_1.clear();

					// Q_temp = Rz1 * Rx
					for (unsigned int i = 0; i < dim; ++i)
							for (unsigned int j = 0; j < dim; ++j)
									for (unsigned int k = 0; k < dim; ++k)
											Q_temp[i][j] += Rz1[i][k] * Rx[k][j];

					// Q1 = Q_temp * Rz2
					for (unsigned int i = 0; i < dim; ++i)
							for (unsigned int j = 0; j < dim; ++j)
									for (unsigned int k = 0; k < dim; ++k)
											Q_1[i][j] += Q_temp[i][k] * Rz2[k][j];
					

					//Rotated normal vector, tn_int
					for (unsigned int i = 0; i < dim; ++i)
							for (unsigned int j = 0; j < dim; ++j)
									tn_int[i] += Q_1[i][j] * tn[j];

					//Rotated direction vector, td_int
					for (unsigned int i = 0; i < dim; ++i)
							for (unsigned int j = 0; j < dim; ++j)
									td_int[i] += Q_1[i][j] * td[j];


					//Rotated twin direction and twin normal vectors
					dealii::Tensor<1, dim> e_X = td_int;
					dealii::Tensor<1, dim> e_Y = tn_int;
					dealii::Tensor<1, dim> e_Z;

					// Compute e_Z = e_X cross e_Y
					e_Z[0] = e_X[1]*e_Y[2] - e_X[2]*e_Y[1];
					e_Z[1] = e_X[2]*e_Y[0] - e_X[0]*e_Y[2];
					e_Z[2] = e_X[0]*e_Y[1] - e_X[1]*e_Y[0];

					// Normalize e_Z
					double norm_eZ = std::sqrt(e_Z*e_Z);
					for (unsigned int i = 0; i < dim; ++i)
							e_Z[i] /= norm_eZ;

					// Construct rotation matrix Q as 3 column vectors
					dealii::Tensor<2, dim> Q;
					for (unsigned int i = 0; i < dim; ++i) {
							Q[i][0] = e_X[i];
							Q[i][1] = e_Y[i];
							Q[i][2] = e_Z[i];
					}
					//Intermediate vector, n_int
					for (unsigned int i = 0; i < 3; ++i)
							for (unsigned int j = 0; j < 3; ++j)
									n[i] += Q[j][i] * nc[j];

					//Rotated unit vector with respect to the twin plane    
					nX = n[0];
					nY = n[1];
					nZ = n[2];

					//Distance from center of the ellipse to a point in the ellipse that intersects the line defined by (nX,nY,nZ)
					//Regularized nZ such that nX^2 + nY^2 + nZ^2 = 1
					nZ_reg = std::sqrt(1.0-nX*nX-nY*nY);
					edist = 1.0/std::sqrt((nX/a0)*(nX/a0) + (nY/b0)*(nY/b0) + (nZ_reg/a0)*(nZ_reg/a0));
													
					scalar_IC =  0.5*(1.0-std::tanh((dist-edist)/(1.0*del0*std::sqrt(edist/a0))));

					// Set the value
					seed_n[v] = seeding_scale_factor * scalar_IC;
				}
			
			variable_list.set_scalar_value_term_RHS(0,seed_n);
			variable_list.set_scalar_gradient_term_RHS(0,constV(0.0)*seed_nx);
		}
	else
	  {
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

      // Computing the outward mobility (L = grad(nvec) dot Ltens dot grad(nvec))
      scalarvalueType_pf L = constV(0.0);
      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              // Mobility tensor (rotated)
              L = L + nvec[i] * nvec[j] * Ltens[i][j];
            }
        }

			//Applying a filter to localize driving force to the twin boundary 
			scalarvalueType_pf strain_df_filter = 1.5*(1.0 - (2.0*n-1.0)*(2.0*n-1.0))*strain_df;

			//Defining the value and gradient terms
			scalarvalueType_pf eq_n = (n-constV(userInputs_pf.dtValue)*L*(mu_twV-strain_df_filter));
			scalargradType_pf eqx_n = -(constV(userInputs_pf.dtValue)*L*kappagradn);

			// --- Submitting the terms for the governing equations ---

			variable_list.set_scalar_value_term_RHS(0,eq_n);
			variable_list.set_scalar_gradient_term_RHS(0,eqx_n);
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
scalarvalueType_pf strain_df_filter = 1.5*(1.0 - (2.0*n-1.0)*(2.0*n-1.0))*strain_df;

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
