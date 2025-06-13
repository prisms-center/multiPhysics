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
        dealii::Tensor<1, dim> e_X = td;
        dealii::Tensor<1, dim> e_Y = tn;
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

        // Compute Q^T
        dealii::Tensor<2, dim> Q_T;
        for (unsigned int i = 0; i < dim; ++i)
            for (unsigned int j = 0; j < dim; ++j)
                Q_T[i][j] = Q[j][i];

        //Compute Ltens_ccref = Q * Lij_tp * Q^T
        dealii::Tensor<2, dim> temp1;
        dealii::Tensor<2, dim> Ltens_ccref;
        temp1.clear();
        for (unsigned int i = 0; i < dim; ++i)
            for (unsigned int j = 0; j < dim; ++j)
                for (unsigned int k = 0; k < dim; ++k)
                    temp1[i][j] += Q[i][k] * Lij_tp[k][j];

        Ltens_ccref.clear();
        for (unsigned int i = 0; i < dim; ++i)
            for (unsigned int j = 0; j < dim; ++j)
                for (unsigned int k = 0; k < dim; ++k)
                    Ltens_ccref[i][j] += temp1[i][k] * Q_T[k][j];
        
        //Gradient energy coefficient tensor
        // Compute K_ccref = Q * Kij_tp * Q^T
        dealii::Tensor<2, dim> temp2;
        dealii::Tensor<2, dim> K_ccref;
        temp2.clear();
        for (unsigned int i = 0; i < dim; ++i)
            for (unsigned int j = 0; j < dim; ++j)
                for (unsigned int k = 0; k < dim; ++k)
                    temp2[i][j] += Q[i][k] * Kij_tp[k][j];
                    
        K_ccref.clear();
        for (unsigned int i = 0; i < dim; ++i)
            for (unsigned int j = 0; j < dim; ++j)
                for (unsigned int k = 0; k < dim; ++k)
                    K_ccref[i][j] += temp2[i][k] * Q_T[k][j];  // Q^T[j][k] = Q[k][j]

        //Tensors in the simulation coordinate system
        // Euler angles in radians (ZXZ convention)
        double pi = 3.14159265358979323846;
        double phi1 = pi*euler_angs[0]/180.0; // e.g., 0.785398 for 45 degrees
        double Phi  = pi*euler_angs[1]/180.0; // e.g., 1.0472   for 60 degrees
        double phi2 = pi*euler_angs[2]/180.0; // e.g., 0.523599 for 30 degrees

        // Step 1: Build the rotation matrices Rz(phi1), Rx(Phi), Rz(phi2)
        dealii::Tensor<2, dim> Rz1, Rx, Rz2;
        Rz1.clear(); Rx.clear(); Rz2.clear();

        // Initialize identity diagonals first
        for (unsigned int i = 0; i < dim; ++i) {
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

        // Q = Q_temp * Rz2
        for (unsigned int i = 0; i < dim; ++i)
            for (unsigned int j = 0; j < dim; ++j)
                for (unsigned int k = 0; k < dim; ++k)
                    Q_1[i][j] += Q_temp[i][k] * Rz2[k][j];

        //Compute Q_1^T
        dealii::Tensor<2, dim> Q_1_T;
        for (unsigned int i = 0; i < dim; ++i)
            for (unsigned int j = 0; j < dim; ++j)
                Q_1_T[i][j] = Q_1[j][i];

        std::cout << "Ltens_ccref" << std::endl;
        std::cout << "(";
        for (unsigned int i = 0; i < dim; ++i)
        {
            std::cout << "(";
            for (unsigned int j = 0; j < dim; ++j)
            {
                std::cout << Ltens_ccref[i][j];
                if (j < dim - 1)
                    std::cout << ",";
            }
            std::cout << ")";
            if (i < dim - 1)
                std::cout << ",";
        }

        std::cout << "Rotation Matrix Q_1" << std::endl;
        std::cout << "(";
        for (unsigned int i = 0; i < dim; ++i)
        {
            std::cout << "(";
            for (unsigned int j = 0; j < dim; ++j)
            {
                std::cout << Q_1[i][j];
                if (j < dim - 1)
                    std::cout << ",";
            }
            std::cout << ")";
            if (i < dim - 1)
                std::cout << ",";
        }
        std::cout << ")" << std::endl;        

        // Transform Ltens_ccref to Ltens and K_ccref to K using

        //Mobility, LTens
        dealii::Tensor<2, dim> temp3;
        temp3.clear(); Ltens.clear();

        for (unsigned int i = 0; i < dim; ++i)
            for (unsigned int j = 0; j < dim; ++j)
                for (unsigned int k = 0; k < dim; ++k)
                    temp3[i][j] += Q_1[i][k] * Ltens_ccref[k][j];  

        // Ltens= temp3 * Q_1
        for (unsigned int i = 0; i < dim; ++i)
            for (unsigned int j = 0; j < dim; ++j)
                for (unsigned int k = 0; k < dim; ++k)
                    Ltens[i][j] += temp3[i][k] * Q_1_T[k][j];

        //Grad Energy coefficient tensor
        dealii::Tensor<2, dim> temp4;
        temp4.clear(); K.clear();
        for (unsigned int i = 0; i < dim; ++i)
            for (unsigned int j = 0; j < dim; ++j)
                for (unsigned int k = 0; k < dim; ++k)
                    temp4[i][j] += Q_1[i][k] * K_ccref[k][j];  

        // Ltens= temp3 * Q_1
        for (unsigned int i = 0; i < dim; ++i)
            for (unsigned int j = 0; j < dim; ++j)
                for (unsigned int k = 0; k < dim; ++k)
                    K[i][j] += temp4[i][k] * Q_1_T[k][j];

        std::cout << "Tensor K (in simulation coordinates): ";
        std::cout << "(";
        for (unsigned int i = 0; i < dim; ++i)
        {
            std::cout << "(";
            for (unsigned int j = 0; j < dim; ++j)
            {
                std::cout << K[i][j];
                if (j < dim - 1)
                    std::cout << ",";
            }
            std::cout << ")";
            if (i < dim - 1)
                std::cout << ",";
        }
        std::cout << ")" << std::endl;

        std::cout << "Tensor Ltens (in simulation coordinates)" << std::endl;
        std::cout << "(";
        for (unsigned int i = 0; i < dim; ++i)
        {
            std::cout << "(";
            for (unsigned int j = 0; j < dim; ++j)
            {
                std::cout << Ltens[i][j];
                if (j < dim - 1)
                    std::cout << ",";
            }
            std::cout << ")";
            if (i < dim - 1)
                std::cout << ",";
        }
        std::cout << ")" << std::endl;

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

    //double L = userInputs_pf.get_model_constant_double("L");
    dealii::Tensor<2,dim> Lij_tp = userInputs_pf.get_model_constant_rank_2_tensor("Lij_tp");
    dealii::Tensor<2,dim> Kij_tp = userInputs_pf.get_model_constant_rank_2_tensor("Kij_tp");
    dealii::Tensor<1,dim> td = userInputs_pf.get_model_constant_rank_1_tensor("td");
    dealii::Tensor<1,dim> tn = userInputs_pf.get_model_constant_rank_1_tensor("tn");
    dealii::Tensor<1,dim> euler_angs = userInputs_pf.get_model_constant_rank_1_tensor("euler_angs");
    double delf_tw = userInputs_pf.get_model_constant_double("delf_tw");
    //double th = userInputs_pf.get_model_constant_double("th");
    double l0 = userInputs_pf.get_model_constant_double("l0");
    double a0 = userInputs_pf.get_model_constant_double("a0");
    double regval = userInputs_pf.get_model_constant_double("regval");

    dealii::Tensor<2,dim> K;
    dealii::Tensor<2,dim> Ltens;
		
    //Average equilibrium interface width
    double del0;
	// ================================================================
};

#endif
