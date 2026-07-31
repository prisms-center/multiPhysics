#include "../../../include/crystalPlasticity.h"

template <int dim> void crystalPlasticity<dim>::copy_current_op_to_old() {
    twinfraction_conv = this->twinfraction_iter1;
}

template <int dim> void crystalPlasticity<dim>::update_twin_mask() {
    double local_op_change = 0.0;
    double delta_orderparam, tr;

    Vector<double> rot1(dim), rot_new(dim);
    FullMatrix<double> rotmat(dim, dim), rotmat_twin(dim, dim), temp1(dim, dim), FP_t(dim, dim);
    FullMatrix<double> Fstar_crystal(dim, dim);
    FullMatrix<double> Fstar_sample(dim, dim), Fstar(dim, dim);
    FullMatrix<double> Fe_new(dim, dim), Fp_new(dim, dim), Fp_new_inv(dim, dim);
    rotmat = 0.0;

    QGauss<dim> quadrature(this->userInputs_cp.quadOrder);
    FEValues<dim> fe_values(this->FE, quadrature,
                          update_quadrature_points | update_gradients |
                              update_JxW_values);
    const unsigned int num_quad_points = quadrature.size();
    const unsigned int dofs_per_cell = this->FE.dofs_per_cell;
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    if (this->userInputs_cp.flagTaylorModel)
      {
        if (initCalled == false)
          {
            if (this->userInputs_cp.enableAdvancedTwinModel)
              {
                init2(num_quad_points);
              } 
            else 
              {
                init(num_quad_points);
              }
          }
      }
    
    // Determine the number of slip systems
    unsigned int n_slip_systemsWOtwin = this->userInputs_cp.numSlipSystems1;

    // Copy "single phase" slip/twin directions/normals
    m_alpha = m_alpha_SinglePhase;
    n_alpha = n_alpha_SinglePhase;

    // Loop over elements
    unsigned int cellID = 0;
    typename DoFHandler<dim>::active_cell_iterator cell = this->dofHandler
                                                            .begin_active(),
                                                 endc = this->dofHandler.end();

    for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);
            // loop over quadrature points
            cell->set_user_index(fe_values.get_cell()->user_index());
            cell->get_dof_indices(local_dof_indices);

            Vector<double> Ulocal(dofs_per_cell);

            if (!this->userInputs_cp.flagTaylorModel)
              {
                for (unsigned int i = 0; i < dofs_per_cell; i++)
                  {
                    Ulocal[i] = this->solutionWithGhosts[local_dof_indices[i]];
                  }
              }

            for (unsigned int q = 0; q < num_quad_points; ++q)
              {
                delta_orderparam = std::max(this->twinfraction_iter1[cellID][q][0] - twinfraction_conv[cellID][q][0], 0.0);

                local_op_change += delta_orderparam;

                if (delta_orderparam > 0.0)
                  {
                    // Get deformation gradient
                    F = 0.0;
                    if (this->userInputs_cp.flagTaylorModel)
                      {
                        F = this->Fprev;
                      }
                    else 
                      {
                        for (unsigned int d = 0; d < dofs_per_cell; ++d)
                          {
                            unsigned int i =
                                fe_values.get_fe().system_to_component_index(d).first;
                            for (unsigned int j = 0; j < dim; ++j) {
                            F[i][j] +=
                                Ulocal(d) *
                                fe_values.shape_grad(
                                    d, q)[j]; // u_{i,j}= U(d)*N(d)_{,j}, where d is the DOF
                                                // correonding to the i'th dimension
                            }
                          }
                        for (unsigned int i = 0; i < dim; ++i)
                          {
                            F[i][i] += 1;
                          }
                      }

                    // Update Fp and Fe
                    FP_t = Fp_conv[cellID][q];

                    // Get the current orientation for this point
                    for (unsigned int i = 0; i < dim; i++)
                      {
                        rot1 = rot_conv[cellID][q];
                      }
                    odfpoint(rotmat, rot1);
                    
                    // Calculate Fstar in the crystal:
                    Fstar_crystal.reinit(dim, dim);
                    Fstar_crystal(0,0) = 1.0;
                    Fstar_crystal(1,1) = 1.0;
                    Fstar_crystal(2,2) = 1.0;

                    for (unsigned int i = 0; i < 3; i++)
                      {
                        for (unsigned int j = 0; j < 3; j++)
                          {
                            Fstar_crystal[i][j] += this->userInputs_cp.twinShear1 * m_alpha[n_slip_systemsWOtwin][i] * n_alpha[n_slip_systemsWOtwin][j];
                          }
                      }

                    // Convert Fstar to sample coords
                    rotmat.mmult(temp1, Fstar_crystal);
                    temp1.mTmult(Fstar_sample, rotmat);

                    // Set Fstar based on the order parameter
                    Fstar = Fstar_sample;
                    Fstar *= delta_orderparam;
                    Fstar(0,0) += (1 - delta_orderparam);
                    Fstar(1,1) += (1 - delta_orderparam);
                    Fstar(2,2) += (1 - delta_orderparam);

                    // Adjust Fe and Fp
                    FP_t.mmult(Fp_new, Fstar);
                    Fp_new_inv.invert(Fp_new);
                    F.mmult(Fe_new, Fp_new_inv);

                    Fp_conv[cellID][q] = Fp_new;
                    Fe_conv[cellID][q] = Fe_new;
                  }
                
                // Reorient if the twin threshold has been reached
                if (!this->twin_mask[cellID][q][0]
                    && this->twinfraction_iter1[cellID][q][0] > this->userInputs_cp.MPtwinLowerThresholdFraction1)
                  {
                    // Get the current orientation for this point
                    for (unsigned int i = 0; i < dim; i++)
                      {
                        rot1 = rot_conv[cellID][q];
                      }
                    odfpoint(rotmat, rot1);

                    temp1.reinit(dim, dim);
                    temp1[0][0] = 1.0;
                    temp1[1][1] = 1.0;
                    temp1[2][2] = 1.0;
                    for (unsigned int i = 0; i < 3; i++)
                      {
                        for (unsigned int j = 0; j < 3; j++)
                          {
                            temp1[i][j] -= 2*n_alpha[n_slip_systemsWOtwin][i]*n_alpha[n_slip_systemsWOtwin][j];
                          }
                      }
                    rotmat.mmult(rotmat_twin, temp1);
                    tr = rotmat_twin.trace();

                    rot_new = 0.0;
                    rot_new[0] = (-1 / (1 + tr)) * (rotmat_twin(1,2) - rotmat_twin(2,1));
                    rot_new[1] = (-1 / (1 + tr)) * (rotmat_twin(2,0) - rotmat_twin(0,2));
                    rot_new[2] = (-1 / (1 + tr)) * (rotmat_twin(0,1) - rotmat_twin(1,0));

                    // Very large Rodrigues vector norm leads to NaN or Inf, so cap the norm at 10000
                    double rnew_Norm, max_rnew_Norm;
                    max_rnew_Norm = 10000;
                    rnew_Norm = sqrt(rot_new(0)*rot_new(0) + rot_new(1)*rot_new(1) + rot_new(2)*rot_new(2));

                    if (rnew_Norm > max_rnew_Norm)
                      {
                        rot_new(0) = rot_new(0) * max_rnew_Norm / rnew_Norm;
                        rot_new(1) = rot_new(1) * max_rnew_Norm / rnew_Norm;
                        rot_new(2) = rot_new(2) * max_rnew_Norm / rnew_Norm;
                      }

                    for (unsigned int i = 0; i < dim; i++)
                      {
                        rot_conv[cellID][q][i] = rot_new[i];
                        rot_iter[cellID][q][i] = rot_new[i];
                      }
                    
                    this->twin_mask[cellID][q][0] = true;
                  }
              }

            cellID++;

          }
      }
    
    double total_op_change = Utilities::MPI::sum(local_op_change, this->mpi_communicator);
    this->pcout << "Integrated change in the order parameter: " << total_op_change << std::endl;

}

template <int dim> void crystalPlasticity<dim>::update_twin_df() {
    QGauss<dim> quadrature(this->userInputs_cp.quadOrder);
    const unsigned int num_quad_points = quadrature.size();

    if (this->userInputs_cp.flagTaylorModel)
      {
        if (initCalled == false)
          {
            if (this->userInputs_cp.enableAdvancedTwinModel)
              {
                init2(num_quad_points);
              } 
            else 
              {
                init(num_quad_points);
              }
          }
      }
  
    // Loop over elements
    unsigned int cellID = 0;
    typename DoFHandler<dim>::active_cell_iterator cell = this->dofHandler
                                                            .begin_active(),
                                                 endc = this->dofHandler.end();

    for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
          {
            for (unsigned int q = 0; q < num_quad_points; ++q)
              {
                this->postprocessValues(cellID, q, 3, 0) = energy[cellID][q][0]; 
              }
            cellID++;
          }
      }
    
    MultiPhysicsBVP<dim, 1>::projection();
}

#include "../../../include/crystalPlasticity_template_instantiations.h"