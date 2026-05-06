#include "../../../include/crystalPlasticity.h"

template <int dim>
void
crystalPlasticity<dim>::reorient_active_zone()
{
  unsigned int CheckBufferRegion, dimBuffer;
  double lowerBuffer, upperBuffer;
  Point<dim> pnt2;
  QGauss<dim>        quadrature(this->userInputs_cp.quadOrder);
  FEValues<dim>      fe_values(this->FE,
                               quadrature,
                               update_quadrature_points | update_gradients |
                                 update_JxW_values);
  const unsigned int num_quad_points = quadrature.size();
  const unsigned int dofs_per_cell   = this->FE.dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  // loop over elements
  unsigned int                                   cellID = 0;
  typename DoFHandler<dim>::active_cell_iterator cell   = this->dofHandler.begin_active(),
                                                 endc   = this->dofHandler.end();
  for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          // loop over quadrature points
          cell->set_user_index(fe_values.get_cell()->user_index());
          cell->get_dof_indices(local_dof_indices);

          //////////Buffer layer feature/////////////
          if (this->userInputs_cp.flagBufferLayer)
            {
              pnt2        = cell->center();
              dimBuffer   = this->userInputs_cp.dimBufferLayer;
              lowerBuffer = this->userInputs_cp.lowerBufferLayer;
              upperBuffer = this->userInputs_cp.upperBufferLayer;
              if ((pnt2[dimBuffer] >= lowerBuffer) && (pnt2[dimBuffer] <= upperBuffer))
                {
                  CheckBufferRegion = 1;
                }
              else
                {
                  CheckBufferRegion = 0;
                }
            }
          else
            {
              CheckBufferRegion = 1;
            }
          /////////////////////////////////////////////////////

          for (unsigned int q = 0; q < num_quad_points; ++q)
            {
              // Reorient the active zone
              this->reoriented_zone[cellID][q] = this->reoriented_zone[cellID][q] || this->active_zone[cellID][q];
            }
          cellID++;
        }
    }
}

template <int dim>
int
crystalPlasticity<dim>::atr_calc(double active_zone_threshold)
{
  unsigned int CheckBufferRegion, dimBuffer;
  double lowerBuffer, upperBuffer;
  unsigned int local_quad_points_in_new_atr = 0;
  Point<dim> pnt2;
  QGauss<dim>        quadrature(this->userInputs_cp.quadOrder);
  FEValues<dim>      fe_values(this->FE,
                               quadrature,
                               update_quadrature_points | update_gradients |
                                 update_JxW_values);
  const unsigned int num_quad_points = quadrature.size();
  const unsigned int dofs_per_cell   = this->FE.dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  // loop over elements
  unsigned int                                   cellID = 0;
  typename DoFHandler<dim>::active_cell_iterator cell   = this->dofHandler.begin_active(),
                                                 endc   = this->dofHandler.end();
  for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          // loop over quadrature points
          cell->set_user_index(fe_values.get_cell()->user_index());
          cell->get_dof_indices(local_dof_indices);

          //////////Buffer layer feature/////////////
          if (this->userInputs_cp.flagBufferLayer)
            {
              pnt2        = cell->center();
              dimBuffer   = this->userInputs_cp.dimBufferLayer;
              lowerBuffer = this->userInputs_cp.lowerBufferLayer;
              upperBuffer = this->userInputs_cp.upperBufferLayer;
              if ((pnt2[dimBuffer] >= lowerBuffer) && (pnt2[dimBuffer] <= upperBuffer))
                {
                  CheckBufferRegion = 1;
                }
              else
                {
                  CheckBufferRegion = 0;
                }
            }
          else
            {
              CheckBufferRegion = 1;
            }
          /////////////////////////////////////////////////////

          for (unsigned int q = 0; q < num_quad_points; ++q)
            {
              std::vector<double> ttwinvf(n_twin_systems);
              std::vector<double> ttwinvf1(n_twin_systems);
              ttwinvf  = twinfraction_conv[cellID][q];
              ttwinvf1 = this->twinfraction_iter1[cellID][q];

              // TODO: this assumes only 1 twin system
              bool atr_old = ttwinvf[0] > active_zone_threshold;
              bool atr_new = ttwinvf1[0] > active_zone_threshold;

              // The new active zone should include points above the threshold
              // but exclude points above the threshold in the old active zone.
              this->active_zone[cellID][q] = atr_new && !atr_old;
              if (this->active_zone[cellID][q])
                {
                  local_quad_points_in_new_atr++;
                }
            }

          cellID++;
        }
    }

    // TODO: need to return the global number of points in the active zone
    // requires MPI communication
    return Utilities::MPI::sum(local_quad_points_in_new_atr, this->mpi_communicator);
}
