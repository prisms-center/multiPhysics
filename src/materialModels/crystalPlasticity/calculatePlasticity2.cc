#include "../../../include/crystalPlasticity.h"

//////////////////////////////////////////////////////////////////////////
// calculatePlasticity.cc numerically integrates the constitive model.
// This calculatePlasticity.cc is based on the following crystal plasticity model:
// L. Anand, M. Kothari, A computational procedure for rate independent crystal
// plasticity,
//  J. Mech. Phys. Solids, 44 (1996), pp. 525-558.
//
// The following modifications are made to incorporate twinning, where the
// twinned region and twin growth rate is determined by a coupled
// phase-field model:
//   0. Multi-phase capability was removed. This can be re-added in
//      a future version.
//   1. A twin orientation matrix is introduced, which is substituted
//      when rotating the stiffness tensor and slip system Schmid tensor
//      in twinned regions.
//   2. The twin characteristic shear is enfored inside twinned regions,
//      and zero outside twinned regions, based on the PF order parameter.
//
// This copy, calculatePlasticity2, is used to compute the numeric
// tangent modulus and should only be called from calculatePlasticity
//////////////////////////////////////////////////////////////////////////

template <int dim>
void
crystalPlasticity<dim>::calculatePlasticity2(unsigned int cellID,
                                             unsigned int quadPtID,
                                             unsigned int StiffnessCalFlag)
{
  // Determine the number of slip and twin systems
  unsigned int n_slip_systemsWOtwin = this->userInputs_cp.numSlipSystems1;
  n_Tslip_systems = n_slip_systemsWOtwin;
  n_twin_systems = 0;
  if (this->userInputs_cp.enableTwinning1) {
    n_Tslip_systems += this->userInputs_cp.numTwinSystems1;
    n_twin_systems = this->userInputs_cp.numTwinSystems1;
  }
  
  // Initialize temporary vectors and tensors for intermediate steps
  std::vector<double> ttwinvf(n_twin_systems);
  std::vector<double> ttwinvf1(n_twin_systems);
  double delta_orderparam; // TODO: change to std::vector<double>

  // twinfraction_conv contains the phase field order parameter from
  // the previous time increment.
  // twinfraction_iter1 contains the phase field order parameter,
  // interpolated from the PF mesh during interpolate_order_parameter(),
  // see dataTransfer.cc
  ttwinvf  = twinfraction_conv[cellID][quadPtID];
  ttwinvf1 = this->twinfraction_iter1[cellID][quadPtID];
  for (unsigned int i = 0; i < n_twin_systems; i++)
    {
      if (ttwinvf1[i] < 0.0) // order parameter shouldn't be negative
        ttwinvf1[i] = 0.0;
    }
  
  // Copy "single phase" slip/twin directions/normals
  m_alpha = m_alpha_SinglePhase;
  n_alpha = n_alpha_SinglePhase;

  // *************** Declare local variables *************** //
  // Determinants of F and FE at time tau
  double det_F_tau, det_FE_tau;
  
  FullMatrix<double> temp(dim, dim), temp1(dim, dim), temp2(dim, dim), temp3(dim, dim),
    temp4(dim, dim); // Temporary matrices for intermediate steps
  FullMatrix<double> T_tau(dim, dim);
  FullMatrix<double> Fpn_inv(dim, dim), FE_tau_trial(dim, dim), F_trial(dim, dim),
    CE_tau_trial(dim, dim), FP_t2(dim, dim), Ee_tau_trial(dim, dim),
    CE_tau(dim, dim);
  temp  = 0;

  FullMatrix<double> PK1_Stiff(dim * dim, dim * dim), P_tau(dim, dim);
  FullMatrix<double> T_star_tau(dim, dim);
  FullMatrix<double> T_star_tau_trial(dim, dim);

  // Slip resistance and hardening moduli
  initialHardeningModulus.reinit(n_slip_systemsWOtwin);
  saturationStress.reinit(n_slip_systemsWOtwin);
  powerLawExponent.reinit(n_slip_systemsWOtwin);
  initialHardeningModulusTwin.reinit(n_twin_systems);
  saturationStressTwin.reinit(n_twin_systems);
  powerLawExponentTwin.reinit(n_twin_systems);
  for(unsigned int i=0;i<n_slip_systemsWOtwin;i++)
    {
      initialHardeningModulus[i] = this->userInputs_cp.initialHardeningModulus1[i];
      saturationStress[i] = this->userInputs_cp.saturationStress1[i];
      powerLawExponent[i] = this->userInputs_cp.powerLawExponent1[i];
    }
  for(unsigned int i=0;i<n_twin_systems;i++)
    {
      initialHardeningModulusTwin[i] = this->userInputs_cp.initialHardeningModulusTwin1[i];
      saturationStressTwin[i] = this->userInputs_cp.saturationStressTwin1[i];
      powerLawExponentTwin[i] = this->userInputs_cp.powerLawExponentTwin1[i];
    }

  Vector<double> s_beta(n_Tslip_systems), h_beta(n_Tslip_systems);

  // Change in plastic deformation gradient
  FullMatrix<double> del_FP(dim, dim);
  double det_FP_tau;

  // Latent hardening
  FullMatrix<double> h_alpha_beta_t(n_Tslip_systems, n_Tslip_systems);
  q.reinit(n_slip_systemsWOtwin, n_slip_systemsWOtwin);
  q = q_phase1;

  // Tolerance
  double tol1 = this->userInputs_cp.modelStressTolerance;

  // Resolved shear stress
  Vector<double> resolved_shear_tau_trial(n_Tslip_systems);
  Vector<double> resolved_shear_tau(n_Tslip_systems);

  // Slip increments
  Vector<double> x_beta(n_Tslip_systems);
  Vector<double> x_beta_old(n_Tslip_systems);

  // LHS matrix of Ax=b linear system for shear increments
  // (consistency condition)
  FullMatrix<double> A(n_Tslip_systems, n_Tslip_systems);
  FullMatrix<double> A_PA;

  // RHS of Ax=b linear system for shear increments
  Vector<double> b(n_Tslip_systems);

  // Potentially active slip systems
  unsigned int n_PA = 0;
  Vector<double> PA, PA_temp(1);
  Vector<double> active;

  // Elastic and Plastic deformation gradients from time t
  // Time t is current step, time tau is next step which we are solving for
  FullMatrix<double> FE_t(dim, dim), FP_t(dim, dim);

  // Slip resistance
  Vector<double> s_alpha_t(n_Tslip_systems);
  Vector<double> s_alpha_tau(n_Tslip_systems);

  // Rotation of the parent crystal at the current point (Rodrigues vector)
  Vector<double> rot1(dim);
  FullMatrix<double> rotmat(dim, dim);
  rotmat = 0.0;

  // Rotation of the twin at the current point (Rodrigues vector)
  Vector<double> rot_twin(dim);
  FullMatrix<double> rotmat_twin(dim, dim);
  rotmat_twin = 0.0;

  // Schmid Tensors
  FullMatrix<double> SCHMID_TENSOR(n_Tslip_systems * dim, dim);
  
  // B = symm(FE_tau_trial' * FE_tau_trial * S_alpha)
  FullMatrix<double> B(n_Tslip_systems * dim, dim);

  // Vectors for slip/twin direction and normal (respectively) used in intermediate calculations
  Vector<double> m1(dim), n1(dim);

  // Elastic Modulus
  FullMatrix<double> Dmat2(2 * dim, 2 * dim), TM(dim * dim, dim * dim);
  FullMatrix<double> ElasticityTensor(2 * dim, 2 * dim);

  Vector<double> vec1(6), vec2(9);
  vec1(0) = 0;
  vec1(1) = 5;
  vec1(2) = 4;
  vec1(3) = 1;
  vec1(4) = 3;
  vec1(5) = 2;
  vec2(0) = 0;
  vec2(1) = 5;
  vec2(2) = 4;
  vec2(3) = 5;
  vec2(4) = 1;
  vec2(5) = 3;
  vec2(6) = 4;
  vec2(7) = 3;
  vec2(8) = 2;



  // ******** Actual constitutive calculations start here ******** //
  // Deformation gradients
  // Use the perturbed deformation gradient
  F_tau = F_per;

  FE_t = Fe_conv[cellID][quadPtID];
  FP_t = Fp_conv[cellID][quadPtID];

  // Slip resistance and backstress from last step (current state)
  for (unsigned int i = 0; i < n_Tslip_systems; i++)
    {
      s_alpha_t[i] = s_alpha_conv[cellID][quadPtID][i];

      // Veera's changes Feb 18
      /*
      if (i >= n_slip_systemsWOtwin
        && ttwinvf1[i - n_slip_systemsWOtwin] < this->userInputs_cp.MPtwinLowerThresholdFraction1)
        {
          // TODO: determine if this is working correctly, or whether it needs to
          // be done differently.
          s_alpha_t[i] = 1e9;
        }*/
    }
  
  // Parent and twin orientations
  // TODO: switch to using rotnew_conv and rotnew_iter and reorient2()?
  //       Check whether reorient is being correctly performed, either
  //       in this function or in updateAfterIncrement()
  for (unsigned int i = 0; i < dim; i++)
    {
      rot1[i] = rot_conv[cellID][quadPtID][i];
    }
  
  odfpoint(rotmat, rot1);

  // Calculate rot_twin (for the first twin system only, right now).
  // Assuming the crystal lattice is centrosymmetric, the twin misorientation
  // relationship (reflection across twin plane) is equivalent to 180 degrees
  // rotation about the twin plane normal.
  // The quaternion for a 180deg rotation about an axis n_i is simply:
  //    q = [0, n_1, n_2, n_3]
  // Calculate this quaternion, multiply it to the parent orientation quaternion
  // using the Hamilton product of quaternions, then convert back.
  // TODO: modify to allow multiple twin systems
  // TODO: calculate twin misorientation during initialization so it's not
  //       recomputed at every point at every time step
  
  /*  NOTE: This commented code corresponds to the above description,
   *        where the twin misorientation is considered to be a 180 degree
   *        rotation about the twin plane normal.
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
  */

  /*   Alternative twin misorientation: consider the twin misorientation to
   *        be a reflection across the twin plane (improper rotation).
   */
  FullMatrix<double> temprot(3,3);
  temprot[0][0] = 1.0;
  temprot[1][1] = 1.0;
  temprot[2][2] = 1.0;
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      temprot[i][j] -= 2*n_alpha[n_slip_systemsWOtwin][i]*n_alpha[n_slip_systemsWOtwin][j];
    }
  }
  rotmat.mmult(rotmat_twin, temprot);

  // Copy the elastic stiffness tensor from user inputs
  // TODO: can this be moved to initialization, since it does not change w/ time?
  elasticStiffnessMatrix.reinit(2*dim, 2*dim);
  for (unsigned int i = 0; i < 6; i++)
    {
      for (unsigned int j = 0; j < 6; j++)
        {
          elasticStiffnessMatrix[i][j] = this->userInputs_cp.elasticStiffness1[i][j];
        }
    }
  
  // Determine if the current point is inside a twin
  // TODO: if there are multiple twin systems, this needs to be re-worked
  // Right now, if any twin systems have ttwinvf>threshold, the voxel is
  // considered twinned and the orientation for twin system 1 is used
  bool isTwinned = false;
  for (unsigned int i = 0; i < n_twin_systems; i++)
    {
      if (ttwinvf1[i] > this->userInputs_cp.MPtwinLowerThresholdFraction1)
        {
          isTwinned = true;
        }
    }

  // Rotate the elastic moduli based on the twin or parent grain orientation
  if (isTwinned)
    {
      elasticmoduli(Dmat2, rotmat_twin, elasticStiffnessMatrix);
    }
  else
    {
      elasticmoduli(Dmat2, rotmat, elasticStiffnessMatrix);
    }
  
  Dmat.reinit(6, 6);
  Dmat = 0.0;

  for (unsigned int i = 0; i < 6; i++)
    {
      for (unsigned int j = 0; j < 6; j++)
        {
          Dmat[i][j] = Dmat2[i][j];
        }
    }
  
  for (unsigned int i = 0; i < 6; i++)
    {
      for (unsigned int j = 3; j < 6; j++)
        {
          Dmat[i][j] = 2 * Dmat[i][j];
        }
    }

  FP_tau = FP_t;
  Fpn_inv = 0.0;
  Fpn_inv.invert(FP_t);
  s_alpha_tau = s_alpha_t;
  FE_tau_trial = 0.0;
  F_trial = 0.0;
  F_tau.mmult(FE_tau_trial, Fpn_inv);
  F_trial = FE_tau_trial;

  temp.reinit(dim, dim);
  temp = 0.0;
  temp = FE_tau_trial;
  FE_tau_trial.Tmmult(CE_tau_trial, temp);
  Ee_tau_trial = CE_tau_trial;
  temp = IdentityMatrix(dim);

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          Ee_tau_trial[i][j] = 0.5 * (Ee_tau_trial[i][j] - temp[i][j]);
        }
    }
  
  //% % % % % STEP 1 % % % % %
  // Calculate the Schmid tensors
  for (unsigned int i = 0; i < n_Tslip_systems; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          m1(j) = m_alpha[i][j];
          n1(j) = n_alpha[i][j];
        }

      temp  = 0.0;
      temp2 = 0.0;
      for (unsigned int j = 0; j < dim; j++)
        {
          for (unsigned int k = 0; k < dim; k++)
            {
              temp[j][k] = m1(j) * n1(k);
            }
        }
      // convert the Schmid tensor to Sample coordinates
      if (isTwinned && i < n_slip_systemsWOtwin)
        {
          // If we are inside the twin, the orientation matrix should
          // change, but only for the slip systems.
          rotmat_twin.mmult(temp2, temp);
          temp2.mTmult(temp, rotmat_twin);
        }
      else
        {
          // If we are outside the twin, OR we are considering the shear on
          // the twin system, we should use the original grain orientation.
          rotmat.mmult(temp2, temp);
          temp2.mTmult(temp, rotmat);
        }
      
      for (unsigned int j = 0; j < dim; j++)
        {
          for (unsigned int k = 0; k < dim; k++)
            {
              SCHMID_TENSOR[dim * i + j][k] = temp[j][k];
            }
        }
      CE_tau_trial.mmult(temp2, temp);
      temp2.symmetrize();
      for (unsigned int j = 0; j < dim; j++)
        {
          for (unsigned int k = 0; k < dim; k++)
            {
              B[dim * i + j][k] = 2 * temp2[j][k];
            }
        }
    }
  
  //% % % % % STEP 2 % % % % %
  // Calculate the trial stress T_star_tau_trial
  Vector<double> tempv1(6), tempv2(6);
  tempv1 = 0.0;
  Dmat.vmult(tempv1, vecform(Ee_tau_trial));
  matform(T_star_tau_trial, tempv1);
  T_star_tau.equ(1.0, T_star_tau_trial);

  // If no plastic deformation occured from t to tau=t+dt, the trial values
  // are already correct. To determine whether plastic deformation occurred,
  // calculate the resolved shear stress on each slip system and compare to
  // the current slip resistance (CRSS + hardening).
  
  //% % % % % STEP 3 % % % % %
  // Calculate the trial resolved shear stress resolved_shear_tau_trial for each slip
  // system

  resolved_shear_tau_trial = 0.0;

  for (unsigned int i = 0; i < n_Tslip_systems; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          for (unsigned int k = 0; k < dim; k++)
            {
              resolved_shear_tau_trial(i) +=
                T_star_tau_trial[j][k] * SCHMID_TENSOR[dim * i + j][k];
            }
        }
    }

  // In the advanced twinning version, here would be an adjustment to
  // resolved_shear_tau_trial for the twin systems, where
  // resolved_shear_tau_trial[j] is set to zero if the order parameter
  // is below the MP twin threshold. This is not needed in this version,
  // since we are instead adjusting the shear increments x_beta, below.
  
  x_beta_old         = 0.0;
  unsigned int iter1 = 1;
  unsigned int flag2 = 0;
  resolved_shear_tau = resolved_shear_tau_trial;

  // Veera's changes Feb 18
  x_beta = 0.0;
  delta_orderparam = ttwinvf1[0] - ttwinvf[0];
  bool twinflag = delta_orderparam > 0.0; // TODO: make this work for detwinning (negative delta_orderparam)
  // </end>

  // ********* Start Nonlinear iteration for Slip increments ********* //

  while (iter1)
    {
      x_beta = 0.0;

      if (iter1 > this->userInputs_cp.modelMaxSlipSearchIterations)
        {
          flag2 = 1;
          break;
        }

      // Skip slip calculation if twinning is happening
      if (!twinflag)
        {
          //% % % % % STEP 4 % % % % %
          n_PA = 0; // Number of active slip systems
          b.reinit(n_Tslip_systems);
          // Determine the set set of the n potentially active slip systems
          // Phil change 3/10: only consider slip. Twinning is determined by the PF model
          for (unsigned int i = 0; i < n_slip_systemsWOtwin; i++)
            {
              b(i) = fabs(resolved_shear_tau(i)) - s_alpha_tau(i);
              if (b(i) >= tol1)
                {
                  // TODO: decide if PA should be changed to std::vector<double>,
                  // which would eliminate the need to manually update the size
                  // every time a new value is appended
                  if (n_PA == 0)
                    {
                      n_PA = n_PA + 1;
                      PA.reinit(n_PA);
                      PA(0) = i;
                    }
                  else
                    {
                      PA_temp = PA;
                      n_PA    = n_PA + 1;
                      PA.reinit(n_PA);
                      for (unsigned int j = 0; j < (n_PA - 1); j++)
                        {
                          PA(j) = PA_temp(j);
                        }
                      PA(n_PA - 1) = i;
                      PA_temp.reinit(n_PA); //%%%%% Potentially active slip systems
                    }
                }
            }

          // If there are no potentially active slip systems, there is no plastic
          // deformation, and the trial values are correct already.
          if (n_PA == 0)
            break;

          Vector<double> b_PA(n_PA);
          b_PA.reinit(n_PA);

          for (unsigned int i = 0; i < n_PA; i++)
            {
              b_PA(i) = b(PA(i));
            }

          if ((b_PA.linfty_norm()) < tol1)
            break;

          if (iter1 > 1)
            {
              temp.reinit(dim, dim);
            }

          //	resolved_shear_tau_trial = resolved_shear_tau;
          FP_t2   = FP_tau;
          Fpn_inv = 0.0;
          Fpn_inv.invert(FP_t2);
          FE_tau = 0.0;
          F_tau.mmult(FE_tau, Fpn_inv);
          temp.reinit(dim, dim);
          temp = 0.0;
          temp = FE_tau;
          FE_tau.Tmmult(CE_tau, temp);

          //% % % % % STEP 5 % % % % %
          // Calculate the shear increments from the consistency condition
          s_beta = s_alpha_tau;

          // Single slip hardening rate
          for (unsigned int i = 0; i < n_slip_systemsWOtwin; i++)
            {
              h_beta(i) = initialHardeningModulus[i] *
                          pow((1 - s_beta(i) / saturationStress[i]), powerLawExponent[i]);
            }

          // Twin hardening rate
          for (unsigned int i = 0; i < n_twin_systems; i++)
            {
              h_beta(n_slip_systemsWOtwin + i) =
                initialHardeningModulusTwin[i] *
                pow((1 - s_beta(n_slip_systemsWOtwin + i) / saturationStressTwin[i]),
                    powerLawExponentTwin[i]);
            }

          for (unsigned int i = 0; i < n_Tslip_systems; i++)
            {
              for (unsigned int j = 0; j < n_Tslip_systems; j++)
                {
                  h_alpha_beta_t[i][j] = q[i][j] * h_beta(j);
                  A[i][j]              = h_alpha_beta_t[i][j];
                }
            }

          for (unsigned int i = 0; i < n_Tslip_systems; i++)
            {
              temp1.reinit(dim, dim);
              temp1 = 0.0;
              for (unsigned int k = 0; k < dim; k++)
                {
                  for (unsigned int l = 0; l < dim; l++)
                    {
                      temp1[k][l] = SCHMID_TENSOR(dim * i + k, l);
                    }
                }
              for (unsigned int j = 0; j < n_Tslip_systems; j++)
                {
                  temp.reinit(dim, dim);
                  temp = 0.0;
                  for (unsigned int k = 0; k < dim; k++)
                    {
                      for (unsigned int l = 0; l < dim; l++)
                        {
                          temp[k][l] = SCHMID_TENSOR(dim * j + k, l);
                        }
                    }
                  temp2.reinit(dim, dim);
                  CE_tau.mmult(temp2, temp);
                  temp2.symmetrize();
                  tempv1 = 0.0;
                  Dmat.vmult(tempv1, vecform(temp2));
                  temp3 = 0.0;
                  matform(temp3, tempv1);
                  // temp3 is symm in Matlab line 94 of constitutive.m

                  for (unsigned int k = 0; k < dim; k++)
                    {
                      for (unsigned int l = 0; l < dim; l++)
                        {
                          if ((resolved_shear_tau(i) * resolved_shear_tau(j)) < 0.0)
                            A[i][j] -= temp1[k][l] * temp3[k][l];
                          else
                            A[i][j] += temp1[k][l] * temp3[k][l];
                        }
                    }
                }
            }
          
          // change A and b to account for already enforced twin shear, allow slip
          // in the twin
          // TODO: make this work when there are multiple twin systems
          

          // Modified slip system search for adding corrective term
          inactive_slip_removal(active,
                                x_beta_old,
                                x_beta,
                                n_PA,
                                n_Tslip_systems,
                                PA,
                                b,
                                A,
                                A_PA,
                                n_slip_systemsWOtwin);
        } // End slip calculation, which only occurs when not twinning
      else
        {
          // Twinning is occuring. Set x_beta
          x_beta(n_slip_systemsWOtwin) = delta_orderparam * this->userInputs_cp.twinShear1;
          
          //	resolved_shear_tau_trial = resolved_shear_tau;
          FP_t2   = FP_tau;
          Fpn_inv = 0.0;
          Fpn_inv.invert(FP_t2);
          FE_tau = 0.0;
          F_tau.mmult(FE_tau, Fpn_inv);
          temp.reinit(dim, dim);
          temp = 0.0;
          temp = FE_tau;
          FE_tau.Tmmult(CE_tau, temp);
        }

      // Step 6: Update the plastic deformation gradient at time tau=t+dt,
      // using the linearized incremental flow rule,
      // Equation (21) of the PRISMS-Plasticity paper.
      temp.reinit(dim, dim);
      del_FP.reinit(dim, dim);
      del_FP = 0.0;
      for (unsigned int i = 0; i < n_Tslip_systems; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              for (unsigned int k = 0; k < dim; k++)
                {
                  temp[j][k] = SCHMID_TENSOR[dim * i + j][k];
                }
            }

          temp2.reinit(dim, dim);
          temp.mmult(temp2, FP_t2);
          for (unsigned int j = 0; j < dim; j++)
            {
              for (unsigned int k = 0; k < dim; k++)
                {
                  if (resolved_shear_tau(i) > 0)
                    FP_tau[j][k] = FP_tau[j][k] + x_beta(i) * temp2[j][k];
                  else
                    FP_tau[j][k] = FP_tau[j][k] - x_beta(i) * temp2[j][k];
                }
            }
        }

      // Step 7: Normalize the plastic deformation gradient such that
      // det(FP_tau) = 1. NOTE: This was commented out in the MATLAB version
      // of the code, so I have commented it out here too.
      /*
      det_FP_tau = FP_tau.determinant();
      for (unsigned int j = 0; j < dim; j++)
        {
          for (unsigned int k = 0; k < dim; k++)
            {
              FP_tau[j][k] = FP_tau[j][k] * pow(det_FP_tau, -1 / 3);
            }
        }
      */

      // % % % % % STEP 8 % % % % %

      for (unsigned int i = 0; i < n_Tslip_systems; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              for (unsigned int k = 0; k < dim; k++)
                {
                  temp[j][k] = SCHMID_TENSOR[dim * i + j][k];
                }
            }

          temp2.reinit(dim, dim);
          CE_tau.mmult(temp2, temp);
          temp3.reinit(dim, dim);
          for (unsigned int j = 0; j < dim; j++)
            {
              for (unsigned int k = 0; k < dim; k++)
                {
                  temp3[j][k] = temp2[j][k] + temp2[k][j];
                }
            }
          temp4.reinit(dim, dim);
          tempv1 = 0.0;
          Dmat.vmult(tempv1, vecform(temp3));
          matform(temp4, tempv1);
          for (unsigned int j = 0; j < dim; j++)
            {
              for (unsigned int k = 0; k < dim; k++)
                {
                  if (resolved_shear_tau_trial(i) > 0)
                    T_star_tau[j][k] = T_star_tau[j][k] - 0.5 * x_beta(i) * temp4[j][k];
                  else
                    T_star_tau[j][k] = T_star_tau[j][k] + 0.5 * x_beta(i) * temp4[j][k];
                }
            }
        }

      resolved_shear_tau = 0.0;
      for (unsigned int i = 0; i < n_Tslip_systems; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              for (unsigned int k = 0; k < dim; k++)
                {
                  resolved_shear_tau(i) +=
                    T_star_tau[j][k] * SCHMID_TENSOR[dim * i + j][k];
                }
            }
        }


      // % % % % % STEP 9 % % % % %

      double h1 = 0;
      for (unsigned int i = 0; i < n_Tslip_systems; i++)
        {
          h1 = 0;
          for (unsigned int j = 0; j < n_Tslip_systems; j++)
            {
              h1 = h1 + h_alpha_beta_t(i, j) * x_beta(j);
            }
          s_alpha_tau(i) = s_alpha_tau(i) + h1;
        }

      for (unsigned int i = 0; i < n_slip_systemsWOtwin; i++)
        {
          if (s_alpha_tau(i) > saturationStress[i])
            {
              s_alpha_tau(i) = saturationStress[i];
            }
        }

      for (unsigned int i = 0; i < n_twin_systems; i++)
        {
          if ((s_alpha_tau(n_slip_systemsWOtwin + i) > (saturationStressTwin[i])) &&
              (twin_conv[cellID][quadPtID] != 1.0))
            {
              s_alpha_tau(n_slip_systemsWOtwin + i) = saturationStressTwin[i];
            }
        }

      iter1 = iter1 + 1;
    }
  
  // ********** End Nonlinear iteration for Slip increments ********** //

  Fpn_inv = 0.0;
  Fpn_inv.invert(FP_tau);
  FE_tau = 0.0;
  F_tau.mmult(FE_tau, Fpn_inv);
  temp.reinit(dim, dim);
  det_FE_tau = FE_tau.determinant();
  FE_tau.mmult(temp, T_star_tau);
  temp.equ(1.0 / det_FE_tau, temp);
  temp.mTmult(T_tau, FE_tau);

  det_F_tau = F_tau.determinant();
  temp.invert(F_tau);
  T_tau.mTmult(P_tau, temp);
  P_tau.equ(det_F_tau, P_tau);
  

  // **************** FINAL UPDATES **************** //

  P_per.reinit(dim,dim);
  P_per=P_tau;

  sres_tau.reinit(n_Tslip_systems);
  sres_tau = s_alpha_tau;

}

#include "../../../include/crystalPlasticity_template_instantiations.h"
