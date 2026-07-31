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
        //Average equilibrium interface width
        del0 = std::sqrt(2.0*Kij_tp[1][1]/delf_tw);
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
    double delf_tw = userInputs_pf.get_model_constant_double("delf_tw");
    double l0 = userInputs_pf.get_model_constant_double("l0");
    double a0 = userInputs_pf.get_model_constant_double("a0");
    double ecc = userInputs_pf.get_model_constant_double("ecc");
    double regval = userInputs_pf.get_model_constant_double("regval");
    double minL = userInputs_pf.get_model_constant_double("minL");

    //Declaring constants
    //Grad. energy coefficient and mobility tensors
    dealii::Tensor<2,dim> K;
    dealii::Tensor<2,dim> Ltens;

    //Average equilibrium interface width
    double del0;

    // A function for this subclass to convert Rodrigues vectors to rotation matrices
    void rodrigues_to_rotmat(dealii::Tensor<2,dim> &orientationMatrix, dealii::Tensor<1,dim> r) const;

    // Containers to store the K and L tensors for each grain ID
    // These must be mutable so as to be assigned during initialization
    mutable std::map<unsigned int, dealii::Tensor<2,dim>> Kij_map = {};
    mutable std::map<unsigned int, dealii::Tensor<2,dim>> Lij_map = {};

};

#endif