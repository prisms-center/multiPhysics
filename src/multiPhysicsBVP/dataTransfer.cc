// Data transfer and interpolation methods for multiPhysicsBVP class

#include "../../include/multiPhysicsBVP.h"
#include "customPDE.h"

#include <deal.II/fe/fe_tools.h>
#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

template <int dim, int degree>
void MultiPhysicsBVP<dim, degree>::interpolate_order_parameter(
    customPDE<dim, 1> &pf_obj,
    const DoFHandler<dim> &dofHandler_Scalar,
    const Quadrature<dim> &quadrature,
    std::vector<std::vector<std::vector<double>>> &twinfraction_iter1,
    FEValues<dim> &fe_values) 
{
    // Obtaining owned and relevant dofs from PF solution
    const auto &owned_dofs = pf_obj.getDofHandlersSet()[0]->locally_owned_dofs();
    IndexSet relevant_dofs;
    DoFTools::extract_locally_relevant_dofs(*pf_obj.getDofHandlersSet()[0], relevant_dofs);

    // Creating a temporary vector to store the current solutionSet data
    dealii::LinearAlgebra::distributed::Vector<double> temp_vector;
    temp_vector.reinit(owned_dofs, relevant_dofs, MPI_COMM_WORLD);

    // Copy data from pf_obj.getSolutionSet()[0] into the temporary vector
    temp_vector = *pf_obj.getSolutionSet()[0];
    temp_vector.update_ghost_values(); // Ensure ghost values are updated

    // Create interpolation function from PF solutionSet
    Functions::FEFieldFunction<dim, vectorType_pf> fe_function_1(
        *pf_obj.getDofHandlersSet()[0], temp_vector);

    // Declare CPFE solution vectors (ghosted and non-ghosted)
    IndexSet own_dofs = dofHandler_Scalar.locally_owned_dofs();
    IndexSet locally_relevant_dofs;
    DoFTools::extract_locally_relevant_dofs(dofHandler_Scalar, locally_relevant_dofs);
    vectorType_cp non_ghosted_solution_cp;
    vectorType_cp solution_cp;

    // Reinitialize vectors for interpolation
    solution_cp.reinit(own_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    non_ghosted_solution_cp.reinit(own_dofs, MPI_COMM_WORLD);

    // Perform the interpolation into non_ghosted_solution_cp
    VectorTools::interpolate(dofHandler_Scalar, fe_function_1, non_ghosted_solution_cp);

    // Copy data to the ghosted vector and compress
    solution_cp = non_ghosted_solution_cp;
    solution_cp.compress(VectorOperation::insert);

    // Interpolation into quadrature points of the CPFE mesh
    std::vector<double> solution_values_cp(quadrature.size());

    // Create a mapping from active_cell_index() to a contiguous local index
    std::map<unsigned int, unsigned int> active_to_local_index;
    unsigned int local_index = 0;
    for (const auto &cell : dofHandler_Scalar.active_cell_iterators())
    {
        if (cell->is_locally_owned())
        {
            active_to_local_index[cell->active_cell_index()] = local_index++;
        }
    }

    // Resize twinfraction_iter1
    twinfraction_iter1.resize(local_index,
                              std::vector<std::vector<double>>(quadrature.size(),
                                                               std::vector<double>(userInputs_cp.numTwinSystems1, 0.0)));

    // Perform interpolation into quadrature values
    for (const auto &cell : dofHandler_Scalar.active_cell_iterators())
    {
        if (cell->is_locally_owned())
        {
            fe_values.reinit(cell);
            fe_values.get_function_values(solution_cp, solution_values_cp);

            // Use the local index from the mapping
            unsigned int local_index = active_to_local_index[cell->active_cell_index()];

            for (unsigned int q = 0; q < quadrature.size(); ++q)
            {
                twinfraction_iter1[local_index][q][0] = solution_values_cp[q];
            }
        }
    }
}


template <int dim, int degree>
void MultiPhysicsBVP<dim, degree>::interpolate_twin_energy(
    customPDE<dim, 1> &pf_obj,
    const DoFHandler<dim> &dofHandler_Scalar) 
{
  // Obtaining own and locally_relevant dofs from CPFE solution
  IndexSet own_dofs = dofHandler_Scalar.locally_owned_dofs();
  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dofHandler_Scalar, locally_relevant_dofs);

  //Creating temporal vector to store the current solutionSet data
  vectorType_cp temp_vector_cp;
  temp_vector_cp.reinit(own_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
  temp_vector_cp = *postFieldsWithGhosts[3];
  temp_vector_cp.update_ghost_values(); // Ensure ghost values are updated

  // Define object fe_function_3 with CPFE vector
  Functions::FEFieldFunction<dim,vectorType_cp> fe_function_3(
      dofHandler_Scalar, temp_vector_cp);
  pcout << "\nCreated fe_function_3 object " << std::endl;

  //Obtaining owned dofs of pf_obj.getDofHandlersSet()[2]
  const auto &owned_dofs = pf_obj.getDofHandlersSet()[2]->locally_owned_dofs();

  //Obtaining relevant dofs of pf_obj.getDofHandlersSet()[2]
  IndexSet relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(*pf_obj.getDofHandlersSet()[2], relevant_dofs);
  //Declare temporary vector solution_pf
  vectorType_pf solution_pf;

  //Reinit solution_pf
  solution_pf.reinit(owned_dofs, relevant_dofs, MPI_COMM_WORLD);
  pcout << "\nReinit of solution_pf successful " << std::endl;

  // Interpolate values from CPFE into solution_pf using fe_function_3
  VectorTools::interpolate(*pf_obj.getDofHandlersSet()[2], fe_function_3, 
                            solution_pf); 

  //Copy values solution_pf into *pf_obj.getSolutionSet()[2]
  *pf_obj.getSolutionSet()[2] = solution_pf;
  pcout << "\nInterpolation of twin energy complete" << std::endl;
  // ***** End of of twin energy complete! ******
}



#include "../../include/multiPhysicsBVP_template_instantiations.h"