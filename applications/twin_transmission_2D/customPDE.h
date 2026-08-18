#ifndef CUSTOMPDE_H
#define CUSTOMPDE_H

#include "matrixFreePDE.h"

// PRISMS-Plasticity headers
#include "userInputParameters_cp.h"

template <int dim, int degree>
class customPDE : public MatrixFreePDE<dim,degree>
{
public:
    // Constructor
    customPDE(userInputParameters_pf<dim> _userInputs_pf,
              userInputParameters_cp & _userInputs_cp)
        : MatrixFreePDE<dim, degree>(_userInputs_pf, _userInputs_cp)
        , userInputs_pf(_userInputs_pf)
        , userInputs_cp(_userInputs_cp)
    {
    //Equilibrium half-interface width in the propagation direction (i.e., the twin plane direction)
    del0 = l/2.0;
    }

    // Function to set the initial conditions (in ICs_and_BCs.h)
    void setInitialCondition(const dealii::Point<dim> &p,
                             const unsigned int index,
                             double & scalar_IC,
                             dealii::Vector<double> & vector_IC) override;

    // Function to set the non-uniform Dirichlet boundary conditions (in ICs_and_BCs.h)
    void setNonUniformDirichletBCs(const dealii::Point<dim> &p,
                                   const unsigned int index,
                                   const unsigned int direction,
                                   const double time,
                                   double & scalar_BC,
                                   dealii::Vector<double> & vector_BC) override;

private:
    #include "typeDefs.h"

    const userInputParameters_pf<dim> userInputs_pf;
    const userInputParameters_cp& userInputs_cp;

    // Function to set the RHS of the governing equations for explicit time dependent equations (in equations.cc)
    void explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double>> & variable_list,
        dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc) const override;

    // Function to set the RHS of the governing equations for all other equations (in equations.cc)
    void nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
        dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const override;

    // Function to set the LHS of the governing equations (in equations.cc)
    void equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
        dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const override;

    // Function to set postprocessing expressions (in postprocess.h)
    void postProcessedFields(const variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
        variableContainer<dim, degree, dealii::VectorizedArray<double>> &pp_variable_list,
        const dealii::Point<dim, dealii::VectorizedArray<double>>        q_point_loc) const override;

    // ================================================================
    // Methods specific to this subclass
    // ================================================================

    // ================================================================
    // Model constants specific to this subclass
    // ================================================================

    dealii::Tensor<2,dim> Lij_tp = userInputs_pf.get_model_constant_rank_2_tensor("Lij_tp");
    dealii::Tensor<2,dim> Kij_tp = userInputs_pf.get_model_constant_rank_2_tensor("Kij_tp");
    dealii::Tensor<1,dim> td = userInputs_pf.get_model_constant_rank_1_tensor("td");
    dealii::Tensor<1,dim> tn = userInputs_pf.get_model_constant_rank_1_tensor("tn");
    dealii::Tensor<1,dim> euler_angs = userInputs_pf.get_model_constant_rank_1_tensor("euler_angs");
    double U = userInputs_pf.get_model_constant_double("U");
    double l = userInputs_pf.get_model_constant_double("l");
    double interface_width = userInputs_pf.get_model_constant_double("interface_width");
    double critical_grad = userInputs_pf.get_model_constant_double("critical_grad");
    double a0 = userInputs_pf.get_model_constant_double("a0");
    double ecc = userInputs_pf.get_model_constant_double("ecc");
    double regval = userInputs_pf.get_model_constant_double("regval");
    double minL = userInputs_pf.get_model_constant_double("minL");
    double fbr_mu = userInputs_pf.get_model_constant_double("fbr_mu");
    double alpha = userInputs_pf.get_model_constant_double("alpha");

    //Declaring constants
    //Grad. energy coefficient and mobility tensors
    dealii::Tensor<2,dim> K;
    dealii::Tensor<2,dim> Ltens;

    double del0; //Equilibrium half-interface width in the propagation direction (i.e., the twin plane direction)   

    // A function for this subclass to convert Rodrigues vectors to rotation matrices
    void rodrigues_to_rotmat(dealii::Tensor<2,dim> &orientationMatrix, dealii::Tensor<1,dim> r) const;

    // Containers to store the K and L tensors for each grain ID
    // These must be mutable so as to be assigned during initialization
    mutable std::map<unsigned int, dealii::Tensor<2,dim>> Kij_map = {};
    mutable std::map<unsigned int, dealii::Tensor<2,dim>> Lij_map = {};

};

// Custom hyperbolic tangent function for dealii vectorized arrays
template <typename Number, std::size_t width>
inline ::dealii::VectorizedArray<Number, width>
tanh(const ::dealii::VectorizedArray<Number, width> &x)
{
  ::dealii::VectorizedArray<Number, width> out;
  for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size(); ++i)
    out[i] = std::tanh(x[i]);
  return out;
}

#endif