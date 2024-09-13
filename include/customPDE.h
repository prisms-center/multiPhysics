#ifndef CUSTOMPDE_H
#define CUSTOMPDE_H

#include "multiPhysicsBVP.h"

template <int dim, int degree>
class customPDE: public MultiPhysicsBVP<dim,degree>
{
public:
    // Original Constructor
		//customPDE(userInputParameters<dim> _userInputs): MatrixFreePDE<dim,degree>(_userInputs) , userInputs(_userInputs) {};
		//New Constructor
		//customPDE(userInputParameters_pf<dim> _userInputs_pf, userInputParameters_cp & _userInputs_cp);
    // Constructor definition inside the class
    customPDE(userInputParameters_pf<dim> _userInputs_pf, userInputParameters_cp & _userInputs_cp)
        : MultiPhysicsBVP<dim, degree>(_userInputs_pf, _userInputs_cp), userInputs_pf(_userInputs_pf) 
    {
      // Optional: Additional constructor logic can go here
      double cth = std::cos(th);
      double sth = std::sin(th);
    
      //Rotation Matrix
      double R[2][2] = {{cth,-sth},{sth,cth}};
      
      //Gradient energy coefficient in the reference frame of the parent phase
      //dealii::Tensor<2,dim> R = ((cth,-sth),(sth,sth));
      for (unsigned int m=0;m<2;m++){
        for (unsigned int n=0;n<2;n++){
          K[m][n]=0.0;
          for(unsigned int i=0;i<2;i++){
            for(unsigned int j=0;j<2;j++){
              K[m][n] = K[m][n] + R[m][i]*R[n][j]*Kij_tp[i][j];
            }
          }
        }
      }
      //Average equilibrium interface width
      del0 = std::sqrt(2.0*(K[0][0]+K[1][1])/delf_tw);
    }
    // Function to set the initial conditions (in ICs_and_BCs.h)
    void setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC);

    // Function to set the non-uniform Dirichlet boundary conditions (in ICs_and_BCs.h)
    void setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC);

		// Override run function
		void run() override {}

		void getElementalValues(FEValues<dim>& fe_values, unsigned int dofs_per_cell, unsigned int num_quad_points, FullMatrix<double>& elementalJacobian, Vector<double>&elementalResidual) override {}

	// ================================================================
	// Methods specific to this subclass
	// ================================================================

  // Function to set postprocessing expressions (in postprocess.h)
	//#ifdef POSTPROCESS_FILE_EXISTS
  void
  postProcessedFields(
    const variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
    variableContainer<dim, degree, dealii::VectorizedArray<double>> &pp_variable_list,
    const dealii::Point<dim, dealii::VectorizedArray<double>>        q_point_loc) const;
	//#endif

private:
	#include "typeDefs.h"

	const userInputParameters_pf<dim> userInputs_pf;

	// Function to set the RHS of the governing equations for explicit time dependent equations (in equations.cc)
    void explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

    // Function to set the RHS of the governing equations for all other equations (in equations.cc)
    void nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

	// Function to set the LHS of the governing equations (in equations.cc)
		void equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;


	// ================================================================
	// Model constants specific to this subclass
	// ================================================================

    double L = userInputs_pf.get_model_constant_double("L");
    dealii::Tensor<2,dim> Kij_tp = userInputs_pf.get_model_constant_rank_2_tensor("Kij_tp");
    double delf_tw = userInputs_pf.get_model_constant_double("delf_tw");
    double th = userInputs_pf.get_model_constant_double("th");
    double l0 = userInputs_pf.get_model_constant_double("l0");
    double a0 = userInputs_pf.get_model_constant_double("a0");
    double ecc = userInputs_pf.get_model_constant_double("ecc");

    dealii::Tensor<2,2> K;
    double cth;
    double sth;
		
    //Average equilibrium interface width
    double del0;
	// ================================================================
};

#endif