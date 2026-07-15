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
  std::vector<bool> twin_init_bool(userInputs_cp.numTwinSystems1);
  for (unsigned int i = 0; i < userInputs_cp.numTwinSystems1; i++)
    {
      twin_init[i] = 0.0;
      twin_init_bool[i] = false;
    }
  twinfraction_iter1.resize(num_local_cells,
                            std::vector<std::vector<double>>(num_quad_points, twin_init));
  dtwinfraction_iter1.resize(num_local_cells,
                             std::vector<std::vector<double>>(num_quad_points,
                                                              twin_init));
  twin_mask.resize(num_local_cells,
                             std::vector<std::vector<bool>>(num_quad_points,
                                                              twin_init_bool));

  // CPFE time-stepping loop STARTS
  currentIncrement_cp          = 0;
  int cp_increment_switch_flag = 1;
  // Flag to track whether seeding is complete
  bool seeding_complete = false;
  double seeding_point_rss = 0.0;
  while (currentIncrement_cp < totalIncrements_cp)
    {
      pcout << "\nCPFE: Current increment number = " << currentIncrement_cp
                << ", Current time = " << currentIncrement_cp * delT
                << ", Time increment = " << delT << std::endl;
      
      // Solve CPFE nonlinear problem
      // call updateBeforeIncrement
      updateBeforeIncrement();

      if (!userInputs_cp.flagTaylorModel)
        {
          // solve time increment
          success = solveNonLinearSystem();
        }
      
      // Check CRSS at seeding point
      // TO-DO: this needs to be re-worked for multiple twin seeds and multiple variants (effort: hard)
      if (!seeding_complete && seeding_point_rss > userInputs_cp.initialSlipResistanceTwin1[0])
        {
          pcout << "\nTwin seeding time reached. Placing seed in CPFE mesh" << std::endl;
          // Place the seed if not already placed, by doing the following (seed already present in PF mesh)
          // Transfer PF data to CPFE
          interpolate_order_parameter(pf_obj,
                                      dofHandler_Scalar,
                                      quadrature,
                                      twinfraction_iter1,
                                      fe_values);
          pcout << "\nInterpolation of n complete" << std::endl;

          // Update twin mask and Fstar
          pcout << "\nUpdating twin mask" << std::endl;
          update_twin_mask();

          // Solve CPFE nonlinear problem
          pcout << "\nResolving mechanical equilibrium..." << std::endl;
          if (!userInputs_cp.flagTaylorModel)
            {
              // solve time increment
              success = solveNonLinearSystem();
            }

          seeding_complete = true;
        }
        
      // PF solve loop:
      if (seeding_complete)
        {
          // If the coupling time has been reached, solve N pf steps
          for (unsigned int pf_step = 0; pf_step < userInputs_pf.increments_pftocpfe; pf_step++)
            {
              // Transfer CPFE data to PF
              interpolate_twin_energy(pf_obj, dofHandler_Scalar);
              pcout << "\nInterpolation of twin energy complete" << std::endl;

              // Solve PF equations to evolve the twin
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

              // In CPFE mesh: copy current OP to old OP, then transfer PF data to CPFE
              copy_current_op_to_old();
              
              interpolate_order_parameter(pf_obj,
                                        dofHandler_Scalar,
                                        quadrature,
                                        twinfraction_iter1,
                                        fe_values);
              pcout << "\nInterpolation of n complete" << std::endl;
              pcout << "\nInterpolation of dndt disabled" << std::endl;
              
              // Solve CPFE nonlinear problem
              pcout << "\nResolving mechanical equilibrium..." << std::endl;
              if (!userInputs_cp.flagTaylorModel)
                {
                  // solve time increment
                  success = solveNonLinearSystem();
                }

              // TODO: Check whether the order parameter changed.
              // If no evolution occured, we can exit this loop and skip forward by
              // the appropriate number of timesteps/output steps
            }
        }
        
      // Commit plastic state: UpdateAfterIncrement
      if ((success) || (userInputs_cp.flagTaylorModel))
        {
          updateAfterIncrement();

          if (!seeding_complete)
            {
              // Evaluate the twin CRSS at the seeding point
              IndexSet own_dofs = dofHandler_Scalar.locally_owned_dofs();
              IndexSet locally_relevant_dofs;
              DoFTools::extract_locally_relevant_dofs(dofHandler_Scalar, locally_relevant_dofs);
              vectorType_cp twin_rss;
              twin_rss.reinit(own_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
              twin_rss = *postFieldsWithGhosts[6]; // Index 6 needs to be the twin RSS. TO-DO: change this so it's not reliant on user-defined output variables
              twin_rss.update_ghost_values();
              Utilities::MPI::RemotePointEvaluation<dim, dim> rpe;
              dealii::Point<dim> eval_point;
              eval_point[0] = userInputs_pf.twin_nucleation_point[0] * userInputs_cp.span[0];
              eval_point[1] = userInputs_pf.twin_nucleation_point[1] * userInputs_cp.span[1];
              eval_point[2] = userInputs_pf.twin_nucleation_point[2] * userInputs_cp.span[2];
              std::vector<Point<dim>> evaluation_points = { eval_point };
              MappingQ1<dim,dim> mapping;
              const std::vector<double> evaluated_values = VectorTools::point_values<1>(mapping, dofHandler_Scalar, twin_rss, evaluation_points, rpe);
              seeding_point_rss = evaluated_values[0];
              pcout << "Twin RSS at seeding point = " << seeding_point_rss << std::endl;
            }

          if (seeding_complete)
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

      // TO-DO: fix this to work with variable seed time
      //        (i.e., seeding based on CRSS, such that we don't know seedingT in advance)
      // Check if we just hit the seeding time. If so, adjust the CPFE timestep
      /*
      if (cp_increment_switch_flag == 1 && currentIncrement_cp * delT >= seedingT)
        {
          cp_increment_switch_flag = 2;
          delT                     = delT_pf_adjust;
          currentIncrement_cp      = std::round(seedingT / delT);
          totalIncrements_cp       = std::round(totalT / delT);

          pcout << "\nCPFE time increment switched to " << delT_pf_adjust
                << std::endl;
          pcout << "\nFrom this point on, the current increment number is calculated "
                    "using the new time increment (Delta t)"
                << std::endl;
          pcout << "\nCPFE: Current increment number = " << currentIncrement_cp
                << ", Current time = " << currentIncrement_cp * delT << std::endl;
        }
      */
      
      // CPFE increment complete
      currentIncrement_cp += 1;

    }

}

#include "../../include/multiPhysicsBVP_template_instantiations.h"
