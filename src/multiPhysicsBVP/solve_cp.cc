// solve method for multiPhysicsBVP class
#include <deal.II/fe/fe_tools.h>
#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include "../../include/matrixFreePDE.h"
#include "../../include/multiPhysicsBVP.h"

// Loop over increments and solve each increment in PF and CPFE

template <int dim, int degree>
void
MultiPhysicsBVP<dim, degree>::solve_cp()
{
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

  // CPFE time-stepping loop STARTS
  currentIncrement_cp          = 0;
  int cp_increment_switch_flag = 1;
  // Flag to track whether seeding is complete
  bool seeding_complete = false;
  while (currentIncrement_cp < totalIncrements_cp)
    {
      if (!seeding_complete && (currentIncrement_cp * delT >= seedingT))
        {
          // Replace CPFE timestep with a scaled-down one
          delT = delT_pf_adjust / ((double) userInputs_cp.stepsForSeeding);

          pcout << "Starting twin seeding loop. Using reduced timestep delT = " << delT
                << std::endl;

          // Perform the seeding loop.
          for (unsigned int seedingstep = 0; seedingstep < userInputs_cp.stepsForSeeding;
               seedingstep++)
            {
              // 1. increment PF
              // Phase-Field seeding step STARTS
              pcout << "seeding increment PF:" << pf_obj.getSeedingIncrement() << "\n";

              // solve time increment (equations.cc must handle the seeding steps)
              pf_obj.solveIncrement(false);

              // Apply constraints and update ghost values
              for (unsigned int fieldIndex = 0; fieldIndex < pf_obj.fields.size();
                   fieldIndex++)
                {
                  pf_obj.getConstraintsDirichletSet()[fieldIndex]->distribute(
                    *pf_obj.getSolutionSet()[fieldIndex]);
                  pf_obj.getConstraintsOtherSet()[fieldIndex]->distribute(
                    *pf_obj.getSolutionSet()[fieldIndex]);
                  pf_obj.getSolutionSet()[fieldIndex]->update_ghost_values();
                }

              pf_obj.getSeedingIncrement() += 1;
              // Phase-Field seeding step ENDS

              // 2. interpolate order parameter
              interpolate_order_parameter(pf_obj,
                                          dofHandler_Scalar,
                                          quadrature,
                                          twinfraction_iter1,
                                          fe_values);
              pcout << "\nInterpolation of n complete" << std::endl;

              // 3. solve CPFE w/ smaller delT
              pcout << "Seeding increment CPFE with delT = " << delT << std::endl;
              // call updateBeforeIncrement
              updateBeforeIncrement();

              if (!userInputs_cp.flagTaylorModel)
                {
                  // solve time increment
                  success = solveNonLinearSystem();
                }

              // call updateAfterIncrement if solve was successful
              if ((success) || (userInputs_cp.flagTaylorModel))
                {
                  updateAfterIncrement();

                  // update totalLoadFactor
                  totalLoadFactor += loadFactorSetByModel;

                  // increase loadFactorSetByModel, if succesiveIncForIncreasingTimeStep
                  // satisfied.
                  successiveIncs++;
                }
              else
                {
                  successiveIncs = 0;
                }
            }
          // 4. Loop complete, change the relevant delT, continue from the top
          pcout << "PF Twin Seeding complete. Resuming solve loop." << std::endl;
          seeding_complete = true;

          if (cp_increment_switch_flag == 1)
            {
              cp_increment_switch_flag = 2;
              delT                     = delT_pf_adjust;
              currentIncrement_cp      = std::round(seedingT / delT);
              totalIncrements_cp       = std::round(totalT / delT);

              // NOTE: the entire seeding loop corresponds to a single CPFE increment,
              //       using the post-seeding delT. Increment here, after the timestep
              //       is adjusted.
              currentIncrement_cp += 1;

              pcout << "\nCPFE time increment switched to " << delT_pf_adjust
                    << std::endl;
              pcout << "\nFrom this point on, the current increment number is calculated "
                       "using the new time increment (Delta t)"
                    << std::endl;
              pcout << "\nCPFE: Current increment number = " << currentIncrement_cp
                    << ", Current time = " << currentIncrement_cp * delT << std::endl;
            }
        }
      else
        {
          pcout << "\nCPFE: Current increment number = " << currentIncrement_cp
                << ", Current time = " << currentIncrement_cp * delT
                << ", Time increment = " << delT << std::endl;
          if ((currentIncrement_cp * delT >= seedingT))
            {
              // ***** Interpolation of order parameter "n" from PF mesh into
              // twin volume fraction CPFE mesh ******
              interpolate_order_parameter(pf_obj,
                                          dofHandler_Scalar,
                                          quadrature,
                                          twinfraction_iter1,
                                          fe_values);
              pcout << "\nInterpolation of n complete" << std::endl;
              pcout << "\nInterpolation of dndt disabled" << std::endl;
            }
          // call updateBeforeIncrement
          updateBeforeIncrement();

          if (!userInputs_cp.flagTaylorModel)
            {
              // solve time increment
              success = solveNonLinearSystem();
            }

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
                  for (unsigned int i = 0; i < userInputs_cp.tabularTimeOutput.size();
                       i++)
                    {
                      tabularTimeInputInc[i] = tabularTimeInputInc[i] / delT;
                    }
                  tabularTimeInputIncInt.resize(userInputs_cp.tabularTimeOutput.size(),
                                                0);
                  /// Converting to an integer always rounds down, even if the fraction
                  /// part is 0.99999999.
                  // Hence, I add 0.1 to make sure we always get the correct integer.
                  for (unsigned int i = 0; i < userInputs_cp.tabularTimeOutput.size();
                       i++)
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
          if (currentIncrement_cp * delT >= timeBeforeC)
            {
              for (unsigned int pf_step = 0; pf_step < userInputs_pf.increments_pftocpfe;
                   pf_step++)
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
            }
          currentIncrement_cp += 1;
        }
    }
}

#include "../../include/multiPhysicsBVP_template_instantiations.h"
