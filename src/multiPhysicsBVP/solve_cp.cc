// solve method for multiPhysicsBVP class
#include <deal.II/fe/fe_tools.h>
#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include "../../include/matrixFreePDE.h"
#include "../../include/multiPhysicsBVP.h"

// Loop over increments and solve each increment in PF and CPFE

/*
 PTR Model - Twin Active Zone Scheme
 Loop procedure:
  0. Place the twin seed and set the initial active zone.
  1. Solve a CPFE step. This performs the FEM calculation to reach
     mechanical equilibrium based on the applied velocity gradient.
     Repeatedly calls constitutive.m in a loop to obtain the Cauchy stress
     and tangent modulus (dT/dF), and only updates the displacement.
  2. Update the CPFE data. This calculates and saves various quantities by
     calling constitutive.m once at each integration point. The saved
     quantities include slip resistances, von Mises stress, etc. This also
     updates the twin volume fraction.
    2a. If the reorientation criteria is met, and reorientation occurs
         during updatedata(), then go to Step 3.
    2b. Otherwise, go back to step 1 and solve another CPFE step.
  3. In the CPFE mesh, copy the current order parameter to the old order
     parameter field.
  4. Transfer the driving force for twinning to phase field and solve the
     phase field equations to evolve the order parameter.
  5. Solve phase field steps (how many? when to stop?)
  6. Transfer the updated order parameter from the phase field mesh to the
    current order parameter field in the CPFE mesh.
  7. Update the active twinning region based on the new order parameter.
  8. Go back to Step 1 and continue solving CPFE
*/

template <int dim, int degree>
void
MultiPhysicsBVP<dim, degree>::solve_cp()
{
  /***************************************************************************
   * STEP 0: Twin Seeding                                                     *
   ***************************************************************************/

  // Section for first phase field step solution BEGINS
  //   Accessing pf_object through the virtual function
  auto &pf_obj = this->get_pf_object();
  pcout << "\npf_obj access successful " << std::endl;

  // log time
  computing_timer_cp.enter_subsection("multiPhysicsBVP: solve");
  pcout << "\nsolving PF (first step)";

  pcout << "\ncurrentIncrement_pf = " << pf_obj.getCurrentIncrement() << "\n\n";

  // Do an initial solve to set the elliptic fields
  pf_obj.solveIncrement(true);

  // Apply constraints and update ghost values
  for (unsigned int fieldIndex = 0; fieldIndex < pf_obj.fields.size(); fieldIndex++)
    {
      pf_obj.getConstraintsDirichletSet()[fieldIndex]->distribute(
        *pf_obj.getSolutionSet()[fieldIndex]);
      pf_obj.getConstraintsOtherSet()[fieldIndex]->distribute(
        *pf_obj.getSolutionSet()[fieldIndex]);
      pf_obj.getSolutionSet()[fieldIndex]->update_ghost_values();
    }

  // Output Result for initial conditions
  pf_obj.getOutputResults();
  // currentOutput++;
  pf_obj.getCurrentOutput() += 1;

  // Increase the current increment from 0 to 1 now that the initial conditions have been
  // output
  // currentIncrement_pf++;
  pf_obj.getCurrentIncrement() += 1;

  // Cycle up to the proper output counter
  while (userInputs_pf.outputTimeStepList.size() > 0 &&
         userInputs_pf.outputTimeStepList[pf_obj.getCurrentOutput()] <
           pf_obj.getCurrentIncrement())
    {
      // currentOutput++;
      pf_obj.getCurrentOutput() += 1;
    }

  // time stepping
  pcout << "\nTime stepping parameters: timeStep: " << userInputs_pf.dtValue
        << "  timeFinal: " << userInputs_pf.finalTime
        << "  timeIncrements: " << userInputs_pf.totalIncrements_pf << "\n";
  // Section for first phase field step solution ENDS

  pcout << "begin solve... CPFE\n\n";

  bool success;
  // load increments
  unsigned int successiveIncs = 0;

  // Setting up interpolation from PF to CPFE mesh
  QGauss<dim>         quadrature(userInputs_cp.quadOrder);
  FEValues<dim>       fe_values(FE_Scalar,
                                quadrature,
                                update_values | update_gradients | update_JxW_values);
  const unsigned int  dofs_per_cell   = FE_Scalar.dofs_per_cell;
  const unsigned int  num_quad_points = quadrature.size();
  unsigned int        num_local_cells = triangulation_cp.n_locally_owned_active_cells();
  std::vector<double> twin_init(userInputs_cp.numTwinSystems1);
  for (unsigned int i = 0; i < userInputs_cp.numTwinSystems1; i++)
    {
      twin_init[i] = 0.0;
    }
  twinfraction_iter1.resize(num_local_cells,
                            std::vector<std::vector<double>>(num_quad_points, twin_init));
  dtwinfraction_iter1.resize(num_local_cells,
                             std::vector<std::vector<double>>(num_quad_points,
                                                              twin_init));

  active_zone.resize(num_local_cells, std::vector<bool>(num_quad_points));
  reoriented_zone.resize(num_local_cells, std::vector<bool>(num_quad_points));
  twin_vf_conv.resize(num_local_cells,
                      std::vector<std::vector<double>>(num_quad_points, twin_init));
  twin_vf_iter.resize(num_local_cells,
                      std::vector<std::vector<double>>(num_quad_points, twin_init));

  // TODO: transfer PF mesh to CPFE here
  // ***** Interpolation of order parameter "n" from PF mesh into
  // twin volume fraction CPFE mesh ******
  interpolate_order_parameter(pf_obj,
                              dofHandler_Scalar,
                              quadrature,
                              twinfraction_iter1,
                              fe_values);
  pcout << "\nInterpolation of n complete" << std::endl;
  pcout << "\nInterpolation of dndt disabled" << std::endl;

  // TODO: set the active zone here
  atr_calc(userInputs_cp.active_zone_threshold);

  // CPFE time-stepping loop STARTS
  currentIncrement_cp = 0;
  while (currentIncrement_cp < totalIncrements_cp)
    {
      /***********************************************************************
       * STEP 1: Solve a CPFE step (solveNonLinearSystem)                     *
       ***********************************************************************/

      pcout << "\nCPFE: Current increment number = " << currentIncrement_cp
            << ", Current time = " << currentIncrement_cp * delT
            << ", Time increment = " << delT << std::endl;

      // call updateBeforeIncrement
      updateBeforeIncrement();

      if (!userInputs_cp.flagTaylorModel)
        {
          // solve time increment
          success = solveNonLinearSystem();
        }

      /***********************************************************************
       * STEP 2: Update CPFE data (updateAfterIncrement)                      *
       ***********************************************************************/

      // call updateAfterIncrement if solve was successful
      if ((success) || (userInputs_cp.flagTaylorModel))
        {
          updateAfterIncrement();

          if (currentIncrement_cp * delT >= timeBeforeC)
            {
              // ***** Interpolation of twin energy from CPFE mesh to PF mesh ******
              interpolate_twin_energy(pf_obj, dofHandler_Scalar);
              pcout << "\nInterpolation of twin energy complete" << std::endl;
            }
          // update totalLoadFactor
          totalLoadFactor += loadFactorSetByModel;

          // increase loadFactorSetByModel, if succesiveIncForIncreasingTimeStep
          // satisfied.
          successiveIncs++;
          // output results to file
          computing_timer_cp.enter_subsection("postprocess");

          //////////////////////TabularOutput Start///////////////
          std::vector<unsigned int> tabularTimeInputIncInt;
          std::vector<double>       tabularTimeInputInc;
          if (userInputs_cp.tabularOutput)
            {
              tabularTimeInputInc = userInputs_cp.tabularTimeOutput;
              for (unsigned int i = 0; i < userInputs_cp.tabularTimeOutput.size(); i++)
                {
                  tabularTimeInputInc[i] = tabularTimeInputInc[i] / delT;
                }
              tabularTimeInputIncInt.resize(userInputs_cp.tabularTimeOutput.size(), 0);
              /// Converting to an integer always rounds down, even if the fraction
              /// part is 0.99999999.
              // Hence, I add 0.1 to make sure we always get the correct integer.
              for (unsigned int i = 0; i < userInputs_cp.tabularTimeOutput.size(); i++)
                {
                  tabularTimeInputIncInt[i] = int(tabularTimeInputInc[i] + 0.1);
                }
            }
          //////////////////////TabularOutput Finish///////////////
          if (((!userInputs_cp.tabularOutput) &&
               ((currentIncrement_cp + 1) % userInputs_cp.skipOutputSteps == 0)) ||
              ((userInputs_cp.tabularOutput) &&
               (std::count(tabularTimeInputIncInt.begin(),
                           tabularTimeInputIncInt.end(),
                           (currentIncrement_cp + 1)) == 1)))
            {
              if (userInputs_cp.writeOutput)
                output();
            }
          computing_timer_cp.leave_subsection("postprocess");
        }
      else
        {
          successiveIncs = 0;
        }

      /***********************************************************************
       * STEP 2a/b: Check the reorientation criteria                          *
       ***********************************************************************/

      cout << "Active zone average twin vf = " << atr_avg_twin_vf << std::endl;

      if (atr_avg_twin_vf > userInputs_cp.reorient_threshold)
        {
          /*******************************************************************
           * STEP 2b: Reorient the active zone                                *
           *******************************************************************/
          cout << "Reorienting the active zone, then solving PF" << std::endl;
          reorient_active_zone();

          /*******************************************************************
           * STEP 3: Copy current to old order parameter in CPFE mesh         *
           *******************************************************************/

          // No need to do anything here. This step is handled during CPFE
          // updateAfterIncrement(), where twinfraction_iter1 is copied into
          // twinfraction_conv

          /*******************************************************************
           * STEP 4: Transfer driving force to PF mesh                        *
           *******************************************************************/

          // ***** Interpolation of twin energy from CPFE mesh to PF mesh ******
          interpolate_twin_energy(pf_obj, dofHandler_Scalar);
          pcout << "\nInterpolation of twin energy complete" << std::endl;

          /*******************************************************************
           * STEP 5: Solve PF equations                                       *
           *******************************************************************/
          bool         new_atr_empty = true;
          unsigned int pf_loop_iter  = 1;
          while (new_atr_empty && pf_loop_iter < userInputs_cp.max_pf_loop_iters)
            {
              pcout << "Evolving PF equations - performing "
                   << userInputs_pf.increments_pftocpfe << " steps" << std::endl;
              for (unsigned int n = 0; n < userInputs_pf.increments_pftocpfe; n++)
                {
                  // Phase-Field regular step STARTS
                  // increment current time
                  pf_obj.getCurrentTime() += userInputs_pf.dtValue;
                  if (pf_obj.getCurrentIncrement() % userInputs_pf.skip_print_steps == 0)
                    {
                      pcout << "\ntime increment PF:" << pf_obj.getCurrentIncrement()
                            << "  time: " << pf_obj.getCurrentTime() << "\n";
                      // pcout << "\ncurrent output PF:" << pf_obj.getCurrentOutput() <<
                      // "\n";
                    }

                  // solve time increment
                  pf_obj.solveIncrement(false);

                  // if (userInputs_pf.outputTimeStepList[pf_obj.getCurrentOutput()] ==
                  // pf_obj.getCurrentIncrement_pf()) { Apply constraints and update ghost
                  // values
                  for (unsigned int fieldIndex = 0; fieldIndex < pf_obj.fields.size();
                       fieldIndex++)
                    {
                      pf_obj.getConstraintsDirichletSet()[fieldIndex]->distribute(
                        *pf_obj.getSolutionSet()[fieldIndex]);
                      pf_obj.getConstraintsOtherSet()[fieldIndex]->distribute(
                        *pf_obj.getSolutionSet()[fieldIndex]);
                      pf_obj.getSolutionSet()[fieldIndex]->update_ghost_values();
                    }

                  if (userInputs_pf.outputTimeStepList[pf_obj.getCurrentOutput()] ==
                      pf_obj.getCurrentIncrement())
                    {
                      // Output Results
                      pf_obj.getOutputResults();
                      pf_obj.getCurrentOutput() += 1;
                    }

                  pf_obj.getCurrentIncrement() += 1;
                  // Phase-Field regular step ENDS
                }

              /***************************************************************
               * STEP 6: Transfer new order parameter from PF to CPFE mesh    *
               ***************************************************************/

              // ***** Interpolation of order parameter "n" from PF mesh into
              // twin volume fraction CPFE mesh ******
              interpolate_order_parameter(pf_obj,
                                          dofHandler_Scalar,
                                          quadrature,
                                          twinfraction_iter1,
                                          fe_values);
              pcout << "\nInterpolation of n complete" << std::endl;
              pcout << "\nInterpolation of dndt disabled" << std::endl;

              /***************************************************************
               * STEP 7: Update the active zone                               *
               ***************************************************************/

              int num_atr_points = atr_calc(userInputs_cp.active_zone_threshold);

              // Check the number of points in the new active zone

              pcout << "New active zone contains " << num_atr_points << " elements."
                   << std::endl;

              if (num_atr_points > 0)
                {
                  // If there are points in the new ATR, exit the while-loop
                  new_atr_empty = false;
                }
              else
                {
                  pcout << "No points in the new active zone: phase field has not evolved "
                          "enough.\n"
                       << "Need to continue running phase field." << std::endl;
                }

              pf_loop_iter++;
            }

          if (pf_loop_iter >= userInputs_cp.max_pf_loop_iters)
            {
              pcout << "Reached the maximum number of phase-field loop iterations without "
                      "the active zone changing."
                   << std::endl;
              break;
            }
        }

      /***********************************************************************
       * STEP 2a: No reorientation. Continue the loop                         *
       ***********************************************************************/

      currentIncrement_cp += 1;
    }
}

#include "../../include/multiPhysicsBVP_template_instantiations.h"
