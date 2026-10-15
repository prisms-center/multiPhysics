#include "../../../include/crystalPlasticity.h"

template <int dim> void crystalPlasticity<dim>::copy_current_op_to_old() {
    twinfraction_conv = this->twinfraction_iter1;
}

template <int dim> void crystalPlasticity<dim>::update_twin_mask() {
    double local_op_change = 0.0;
    double delta_orderparam, tr;

    Vector<double> rot1(dim), rot_new(dim);
    FullMatrix<double> rotmat(dim, dim), rotmat_twin(dim, dim), temp1(dim, dim), FP_t(dim, dim);
    FullMatrix<double> Lstar_crystal(dim, dim);
    FullMatrix<double> Lstar_sample(dim, dim), Fstar(dim, dim);
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

              for (unsigned int idx = 0; idx < this->userInputs_cp.numTwinSystems1; idx++)
                {
                  for (unsigned int q = 0; q < num_quad_points; ++q)
                    {
                      delta_orderparam = std::max(this->twinfraction_iter1[cellID][q][idx] - twinfraction_conv[cellID][q][idx], 0.0);

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

                          // Get the parent orientation for this point
                          for (unsigned int i = 0; i < dim; i++)
                            {
                              rot1 = rot_original[cellID][q];
                            }
                          odfpoint(rotmat, rot1);
                          
                          // Calculate Lstar in the crystal:
                          Lstar_crystal.reinit(dim, dim);
                          
                          for (unsigned int i = 0; i < 3; i++)
                            {
                              for (unsigned int j = 0; j < 3; j++)
                                {
                                  Lstar_crystal[i][j] += this->userInputs_cp.twinShear1 * m_alpha[n_slip_systemsWOtwin + idx][i] * n_alpha[n_slip_systemsWOtwin + idx][j];
                                }
                            }

                          // Convert Lstar to sample coords
                          rotmat.mmult(temp1, Lstar_crystal);
                          temp1.mTmult(Lstar_sample, rotmat);

                          // Set Fstar based on the order parameter
                          Fstar = Lstar_sample;
                          Fstar *= delta_orderparam;
                          Fstar(0,0) += 1.0;
                          Fstar(1,1) += 1.0;
                          Fstar(2,2) += 1.0;

                          // Adjust Fe and Fp
                          Fstar.mmult(Fp_new, FP_t);
                          Fp_new_inv.invert(Fp_new);
                          F.mmult(Fe_new, Fp_new_inv);

                          Fp_conv[cellID][q] = Fp_new;
                          Fe_conv[cellID][q] = Fe_new;
                        }
                      
                      // Reorient if the twin threshold has been reached
                      if (!this->twin_mask[cellID][q][0]
                          && this->twinfraction_iter1[cellID][q][idx] > this->userInputs_cp.MPtwinLowerThresholdFraction1)
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
                                  temp1[i][j] -= 2*n_alpha[n_slip_systemsWOtwin + idx][i]*n_alpha[n_slip_systemsWOtwin + idx][j];
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
                for (unsigned int i = 0; i < this->userInputs_cp.numTwinSystems1; i++)
                  {
                    this->postprocessValues(cellID, q, 26 - i, 0) = energy[cellID][q][i]; 
                  }
              }
            cellID++;
          }
      }
    
    MultiPhysicsBVP<dim, 1>::projection();
}

template <int dim> void crystalPlasticity<dim>::calculate_twin_rss() {
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
                // calculate resolved shear stress
                Vector<double> rss(n_Tslip_systems);
                Vector<double> m1(dim), n1(dim);
                FullMatrix<double> schmidtensor(dim, dim);
                FullMatrix<double> temprot(dim, dim);
                FullMatrix<double> temp(dim, dim);
                rss = 0.0;
                unsigned int n_slip_systemsWOtwin = this->userInputs_cp.numSlipSystems1;
                // TODO: fix the following line to consider multiple twin systems
                bool isTwinned = twinfraction_iter[cellID][q][0] >= this->userInputs_cp.MPtwinLowerThresholdFraction1;

                FullMatrix<double> rotmat(dim, dim), rotmat_twin(dim, dim);
                Vector<double> rot1(dim), rot_twin(dim);
                for (unsigned int i = 0; i < dim; i++)
                  {
                    rot1[i] = rot_iter[cellID][q][i];
                  }
                odfpoint(rotmat, rot1);
                Vector<double> quatprod(4), quat1(4), quat2(4);
                rot_twin = 0.0;
                rod2quat(quat2, rot1);
                quat1[0] = 0.0;
                quat1[1] = n_alpha[n_slip_systemsWOtwin][0];
                quat1[2] = n_alpha[n_slip_systemsWOtwin][1];
                quat1[3] = n_alpha[n_slip_systemsWOtwin][2];
                quatproduct(quatprod, quat2, quat1);
                quat2rod(quatprod, rot_twin);
                odfpoint(rotmat_twin, rot_twin);

                for (unsigned int i = 0; i < n_Tslip_systems; i++)
                  {
                    schmidtensor = 0.0;
                    for (unsigned int j = 0; j < dim; j++)
                      {
                        m1(j) = m_alpha[i][j];
                        n1(j) = n_alpha[i][j];
                      }
                    temp = 0.0;
                    temprot = 0.0;
                    for (unsigned int j = 0; j < dim; j++)
                      {
                        for (unsigned int k = 0; k < dim; k++)
                          {
                            temp[j][k] = m1(j) * n1(k);
                          }
                      }
                    // Convert schmid tensor to sample coordinates
                    if (isTwinned && i < n_slip_systemsWOtwin)
                      {
                        // If we are inside the twin, the orientation matrix should
                        // change, but only for the slip systems.
                        rotmat_twin.mmult(temprot, temp);
                        temprot.mTmult(temp, rotmat_twin);
                      }
                    else
                      {
                        // If we are outside the twin, OR we are considering the shear on
                        // the twin system, we should use the original grain orientation.
                        rotmat.mmult(temprot, temp);
                        temprot.mTmult(temp, rotmat);
                      }

                    for (unsigned int j = 0; j < dim; j++)
                      {
                        for (unsigned int k = 0; k < dim; k++)
                          {
                            schmidtensor[j][k] = temp[j][k];
                          }
                      }

                    for (unsigned int j = 0; j < dim; j++)
                      {
                        for (unsigned int k = 0; k < dim; k++)
                          {
                            rss(i) += T[j][k] * schmidtensor[j][k];
                          }
                      }

                    // TODO: change this function to not calculate all the slip system RSS
                    // also, modify this to work with multiple twin variants
                    this->postprocessValues(cellID, q, 6, 0) = rss[12];
                  }
              }
            cellID++;
          }
      }
    
    MultiPhysicsBVP<dim, 1>::projection();
}

#include "../../../include/crystalPlasticity_template_instantiations.h"