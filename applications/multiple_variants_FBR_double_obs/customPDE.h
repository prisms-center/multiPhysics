#ifndef CUSTOMPDE_H
#define CUSTOMPDE_H

#include <iostream>
#include <fstream>

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
        del0 = l / 2.0;

        // Read the CPFE twin normal and direction files.
        // TODO: access these values from CPFE instead of reading the files in two places
        // For now, it's easier to just read the files again here, copying code from crystalPlasticity::init

        if (userInputs_cp.numTwinSystems1 != num_twin_variants)
          {
            std::cout << "ERROR: number of twin systems in the CPFE parameters file is not equal to the declared number of twin variants in customPFE.h" << std::endl;
            std::cout << "customPFE::num_twin_variants=" << num_twin_variants << std::endl;
            std::cout << "userInputs_cp.numTwicSystems1=" << userInputs_cp.numTwinSystems1 << std::endl;
            exit(1);
          }
        
        dealii::Vector<double> t_init(dim);
        tn_list.resize(num_twin_variants, t_init);
        td_list.resize(num_twin_variants, t_init);

        std::string line;
        double n_norm, m_norm;

        //open data file to read twin normals
        std::ifstream twinNormalsDataFile(userInputs_cp.twinNormalsFile1);
        //read data
        unsigned int id = 0;
        if (twinNormalsDataFile.is_open())
          {
            while (getline(twinNormalsDataFile, line) && id < num_twin_variants)
              {
                std::stringstream ss(line);
                ss >> tn_list[id][0];
                ss >> tn_list[id][1];
                ss >> tn_list[id][2];
                n_norm = 0;
                n_norm = n_norm + tn_list[id][0] * tn_list[id][0];
                n_norm = n_norm + tn_list[id][1] * tn_list[id][1];
                n_norm = n_norm + tn_list[id][2] * tn_list[id][2];
                n_norm = sqrt(n_norm);
                tn_list[id][0] = tn_list[id][0] / n_norm;
                tn_list[id][1] = tn_list[id][1] / n_norm;
                tn_list[id][2] = tn_list[id][2] / n_norm;
                id++;
              }
          }
        else
          {
            std::cout << "Unable to open twin normals file\n";
            exit(1);
          }

        //open data file to read twin directions
        std::ifstream twinDirectionsDataFile(userInputs_cp.twinDirectionsFile1);
        //read data
        id = 0;
        if (twinDirectionsDataFile.is_open())
          {
            //read data
            while (getline(twinDirectionsDataFile, line) && id < num_twin_variants)
              {
                std::stringstream ss(line);
                ss >> td_list[id][0];
                ss >> td_list[id][1];
                ss >> td_list[id][2];
                m_norm = 0 ;
                m_norm = m_norm + td_list[id][0] * td_list[id][0];
                m_norm = m_norm + td_list[id][1] * td_list[id][1];
                m_norm = m_norm + td_list[id][2] * td_list[id][2];
                m_norm = sqrt(m_norm) ;
                td_list[id][0] = td_list[id][0] / m_norm;
                td_list[id][1] = td_list[id][1] / m_norm;
                td_list[id][2] = td_list[id][2] / m_norm;
                id++;
              }
          }
        else
          {
            std::cout << "Unable to open twin directions file\n";
            exit(1);
          }
        
        for (unsigned int i = 0; i < num_twin_variants; i++)
          {
            this->pcout << "tn_list[" << i << "][0] = " << tn_list[i][0] << std::endl;
            this->pcout << "tn_list[" << i << "][1] = " << tn_list[i][1] << std::endl;
            this->pcout << "tn_list[" << i << "][2] = " << tn_list[i][2] << std::endl;
          }

        for (unsigned int i = 0; i < num_twin_variants; i++)
          {
            this->pcout << "td_list[" << i << "][0] = " << td_list[i][0] << std::endl;
            this->pcout << "td_list[" << i << "][1] = " << td_list[i][1] << std::endl;
            this->pcout << "td_list[" << i << "][2] = " << td_list[i][2] << std::endl;
          }

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

    //Declaring constants
    //Grad. energy coefficient and mobility tensors
    dealii::Tensor<2,dim> Lij_tp = userInputs_pf.get_model_constant_rank_2_tensor("Lij_tp");
    dealii::Tensor<2,dim> Kij_tp = userInputs_pf.get_model_constant_rank_2_tensor("Kij_tp");
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
    double twin_interaction_coef = userInputs_pf.get_model_constant_double("twin_interaction_coefficient");

    // Twin normals and directions
    std::vector<dealii::Vector<double>> tn_list;
    std::vector<dealii::Vector<double>> td_list;

    //Average equilibrium interface width
    double del0;

    // A function for this subclass to convert Rodrigues vectors to rotation matrices
    void rodrigues_to_rotmat(dealii::Tensor<2,dim> &orientationMatrix, dealii::Tensor<1,dim> r) const;

    // Containers to store the K and L tensors for each grain ID
    // These must be mutable so as to be assigned during initialization
    mutable std::map<unsigned int, std::vector<dealii::Tensor<2,dim>>> Kij_map = {};
    mutable std::map<unsigned int, std::vector<dealii::Tensor<2,dim>>> Lij_map = {};

    // TODO: Change this to not be hard-coded.
    // Requires PRISMS-PF 4.0, so that the variable attributes need not be known at compile-time
    unsigned int num_twin_variants = 2;
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