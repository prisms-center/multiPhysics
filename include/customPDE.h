#include "multiPhysicsBVP.h"

template <int dim, int degree>
class customPDE: public MultiPhysicsBVP<dim,degree>
{
public:
    // Original Constructor
		//customPDE(userInputParameters<dim> _userInputs): MatrixFreePDE<dim,degree>(_userInputs) , userInputs(_userInputs) {};
		//New Constructor
		//customPDE(userInputParameters_pf<dim> _userInputs_pf): MultiPhysicsBVP(_userInputs_pf, _userInputs_cp) , userInputs_pf(_userInputs_pf) {};
    customPDE(userInputParameters_pf<dim> _userInputs_pf): userInputs_pf(_userInputs_pf) {};
    // Function to set the initial conditions (in ICs_and_BCs.h)
    void setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC);

    // Function to set the non-uniform Dirichlet boundary conditions (in ICs_and_BCs.h)
    void setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC);

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
	// Methods specific to this subclass
	// ================================================================


	// ================================================================
	// Model constants specific to this subclass
	// ================================================================

	double MnV = userInputs_pf.get_model_constant_double("MnV");
	double KnV = userInputs_pf.get_model_constant_double("KnV");

	// ================================================================

};
