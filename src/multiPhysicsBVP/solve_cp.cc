// solve method for multiPhysicsBVP class
#include "../../include/multiPhysicsBVP.h"
#include "customPDE.h"

#include <deal.II/fe/fe_tools.h>
#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

// Loop over increments and solve each increment in PF and CPFE

template <int dim, int degree> void MultiPhysicsBVP<dim, degree>::solve_cp() {

  //Section for first phase field step solution BEGINS
  //  Accessing pf_object through the virtual function
  auto &pf_obj = this->get_pf_object();
  pcout << "\npf_obj access successful " << std::endl;

  //log time
  computing_timer_pf.enter_subsection("multiPhysicsBVP: solve");
  pcout << "\nsolving PF (first step)";

  pcout << "\ncurrentIncrement_pf = " << pf_obj.getCurrentIncrement_pf() << "\n\n";

  // For any nonlinear equation, set the initial guess as the solution to Laplace's equations
  //generatingInitialGuess = true;
  //pf_obj.getSetNonlinearEqInitialGuess();
  //generatingInitialGuess = false;

  // Do an initial solve to set the elliptic fields
  pf_obj.getSolveIncrement(true);
  
  //Apply constraints and update ghost values
  for(unsigned int fieldIndex=0; fieldIndex<pf_obj.fields.size(); fieldIndex++){     
      pf_obj.getConstraintsDirichletSet()[fieldIndex]->distribute(*pf_obj.getSolutionSet()[fieldIndex]);
      pf_obj.getConstraintsOtherSet()[fieldIndex]->distribute(*pf_obj.getSolutionSet()[fieldIndex]);
      pf_obj.getSolutionSet()[fieldIndex]->update_ghost_values();
  }

  //Output Result for initial conditions
  pf_obj.getOutputResults();
  //currentOutput++;
  pf_obj.getCurrentOutput() += 1;

  //Output initial condition checkpoint (uncomment later)
  //if (userInputs_pf.checkpointTimeStepList[currentCheckpoint] == currentIncrement_pf) {
  //    save_checkpoint();
  //    currentCheckpoint++;
  //}

  // Increase the current increment from 0 to 1 now that the initial conditions have been output
  //currentIncrement_pf++;
  pf_obj.getCurrentIncrement_pf() += 1;

  // Cycle up to the proper output counter
  while (userInputs_pf.outputTimeStepList.size() > 0 && userInputs_pf.outputTimeStepList[pf_obj.getCurrentOutput()] < pf_obj.getCurrentIncrement_pf()){
      //currentOutput++;
      pf_obj.getCurrentOutput() += 1;
  }
// Cycle up to the proper checkpoint counter (uncomment later)
  //while (userInputs_pf.checkpointTimeStepList.size() > 0 && userInputs_pf.checkpointTimeStepList[currentCheckpoint] < currentIncrement_pf){
  //    currentCheckpoint++;
  //}
  
  //time stepping
  pcout << "\nTime stepping parameters: timeStep: " << userInputs_pf.dtValue << "  timeFinal: " << userInputs_pf.finalTime << "  timeIncrements: " << userInputs_pf.totalIncrements << "\n";
  //Section for first phase field step solution ENDS

  pcout << "begin solve... CPFE\n\n";

  bool success;
  // load increments
  unsigned int successiveIncs = 0;
  // cp_twinfraction_addition
  // local variables
  QGauss<dim> quadrature(userInputs_cp.quadOrder);
  FEValues<dim> fe_values(FE_Scalar, quadrature,
                          update_values | update_gradients | update_JxW_values);
  const unsigned int dofs_per_cell = FE_Scalar.dofs_per_cell;
  const unsigned int num_quad_points = quadrature.size();
  unsigned int num_local_cells =
      triangulation_cp.n_locally_owned_active_cells();
  std::vector<double> twin_init(userInputs_cp.numTwinSystems1);
  for (unsigned int i = 0; i < userInputs_cp.numTwinSystems1; i++) {
    twin_init[i] = 0.0;
  }
  twinfraction_iter1.resize(num_local_cells, std::vector<std::vector<double>>(
                                                 num_quad_points, twin_init));
  dtwinfraction_iter1.resize(num_local_cells, std::vector<std::vector<double>>(
                                                  num_quad_points, twin_init));

  // Needed for interpolation below
  IndexSet own_dofs = dofHandler_Scalar.locally_owned_dofs();
  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dofHandler_Scalar, locally_relevant_dofs);

  if (userInputs_cp.enableAdaptiveTimeStepping) {
    for (; totalLoadFactor < totalIncrements;) {
      ++currentIncrement_cp;
      loadFactorSetByModel =
          std::min(loadFactorSetByModel, totalIncrements - totalLoadFactor);
      pcout << "\nincrement: " << currentIncrement_cp << std::endl;
      char buffer[100];
      sprintf(buffer,
              "current load factor: %12.6e\ntotal load factor:   %12.6e\n",
              loadFactorSetByModel, totalLoadFactor);
      pcout << buffer;

      // call updateBeforeIncrement, if any
      updateBeforeIncrement();

      if (!userInputs_cp.flagTaylorModel) {
        // solve time increment
        success = solveNonLinearSystem();
      }

      // call updateAfterIncrement, if any
      if ((success) || (userInputs_cp.flagTaylorModel)) {
        updateAfterIncrement();

        // update totalLoadFactor
        totalLoadFactor += loadFactorSetByModel;

        // increase loadFactorSetByModel, if succesiveIncForIncreasingTimeStep
        // satisfied.
        successiveIncs++;

        if (successiveIncs >= userInputs_cp.succesiveIncForIncreasingTimeStep) {
          loadFactorSetByModel *= userInputs_cp.adaptiveLoadIncreaseFactor;
          char buffer1[100];
          sprintf(buffer1,
                  "current increment increased. Restarting increment with "
                  "loadFactorSetByModel: %12.6e\n",
                  loadFactorSetByModel);
          pcout << buffer1;
        }
        computing_timer_cp.enter_subsection("postprocess");

        if (currentIncrement_cp % userInputs_cp.skipOutputSteps == 0)
          if (userInputs_cp.writeOutput)
            output();
        computing_timer_cp.leave_subsection("postprocess");

      } else {
        successiveIncs = 0;
      }
    }
    char buffer[100];
    sprintf(buffer, "\nfinal load factor  : %12.6e\n", totalLoadFactor);
    pcout << buffer;
  } else
    for (; currentIncrement_cp < totalIncrements; ++currentIncrement_cp) {
      pcout << "\nincrement CPFE: " << currentIncrement_cp << ", Time CPFE: "<< currentIncrement_cp*delT << std::endl;

      // Interpolate twin fraction and twinfraction change from PF mesh into
      // CPFE mesh.

      // ***** Interpolation of twin seed (Order parameter "n" ******
      // solutionSet[0]) in phase-field mesh into variable dtwinfraction_iter in
      // CPFE mesh
      // Define object fe_function_1 for phase field data given dof handler and
      // the given solution vector. . This object can evaluate the finite
      // element solution at given points in the domain.

      if ((currentIncrement_cp*delT >=  timeBeforeS) || (currentIncrement_cp == 0)){
        pcout << "\n Passing order parameter from phase field as twinfraction_iter1 " << std::endl; 
        Functions::FEFieldFunction<dim, vectorType_pf> fe_function_1(
            *pf_obj.getDofHandlersSet()[0], *pf_obj.getSolutionSet()[0]);
        pcout << "\nCreated fe_function_1 object " << std::endl;
        // Interpolate into the CPFE domain
        vectorType_cp non_ghosted_solution_cp; // No pointer, just a regular object
        vectorType_cp solution_cp;
        //vectorType_cp non_ghosted_solution_cp;
        //non_ghosted_solution_cp.reinit(solution_cp.locally_owned_elements(), MPI_COMM_WORLD);
        //non_ghosted_solution_cp = solution_cp; // Copy the contents without ghost elements      
        solution_cp.reinit(own_dofs, locally_relevant_dofs, mpi_communicator);
        non_ghosted_solution_cp.reinit(own_dofs, mpi_communicator);
        pcout << "\nreinit of solution_cp successful " << std::endl;
        //VectorTools::interpolate(dofHandler, fe_function_1,
        //                        solution_cp); // Works but not in parallel
        VectorTools::interpolate(dofHandler_Scalar, fe_function_1,
                                non_ghosted_solution_cp); // Works but not in parallel
        solution_cp = non_ghosted_solution_cp;
        solution_cp.compress(VectorOperation::insert);                 
        pcout << "\nInterpolated into solution_cp " << std::endl;
        std::vector<double> solution_values_cp(quadrature.size());
        pcout << "\nDefined solution_values_cp for the cell" << std::endl;
        typename DoFHandler<dim>::active_cell_iterator
            cell = dofHandler_Scalar.begin_active(),
            endc = dofHandler_Scalar.end(), dg_cell = dofHandler_Scalar.begin_active();
        pcout << "\nDefined cell iterator" << std::endl;
        unsigned int cellID = 0;
        for (; cell != endc; ++cell) {
          //pcout << "\nInside cell loop" << std::endl;
          fe_values.reinit(cell);
          //pcout << "\nReinit cell successful" << std::endl;
          fe_values.get_function_values(solution_cp, solution_values_cp);
          for (unsigned int q = 0; q < quadrature.size(); ++q) {
            unsigned int i = Utilities::MPI::this_mpi_process(mpi_communicator);
            // std::cout << "Processor " << i
            //           << ", solution_values_cp=" << solution_values_cp[q]
            //           << std::endl;
            twinfraction_iter1[cellID][q][0] =
                solution_values_cp[q]; // Passing the solution for twin of index
                                      // =0
          }
          if (cell->is_locally_owned()) {
            cellID++;
          }
        }
        pcout << "\nInterpolation of n complete" << std::endl;
        // ****** End of Interpolation of "n" ******

        // ++++++ Interpolation of Order parameter time derivative  "dndt" ++++++
        // solutionSet[1]) in phase-field mesh into variable dtwinfraction_iter in
        // CPFE mesh
        // Define object fe_function_1 for phase field data given dof handler and
        // the given solution vector. . This object can evaluate the finite
        // element solution at given points in the domain.
        pcout << "\n Passing dndt from phase field as dtwinfraction_iter1 " << std::endl; 
        Functions::FEFieldFunction<dim, vectorType_pf> fe_function_2(
            *pf_obj.getDofHandlersSet()[1], *pf_obj.getSolutionSet()[1]);
        pcout << "\nCreated fe_function_2 object " << std::endl;
        // Interpolate into the CPFE domain
        vectorType_cp solution_cp_2; // No pointer, just a regular object
        vectorType_cp non_ghosted_solution_cp_2;
        solution_cp_2.reinit(own_dofs, locally_relevant_dofs, mpi_communicator);
        non_ghosted_solution_cp_2.reinit(own_dofs, mpi_communicator);
        pcout << "\nreinit of solution_cp_2 successful " << std::endl;
        //VectorTools::interpolate(dofHandler, fe_function_1,
        //                        solution_cp); // Works but not in parallel
        VectorTools::interpolate(dofHandler_Scalar, fe_function_2,
                                non_ghosted_solution_cp_2); // Works but not in parallel
        solution_cp_2 = non_ghosted_solution_cp_2;
        solution_cp_2.compress(VectorOperation::insert); 
        pcout << "\nInterpolated into solution_cp_2 " << std::endl;
        std::vector<double> solution_values_cp_2(quadrature.size());
        pcout << "\nDefined solution_values_cp_2 for the cell" << std::endl;
        typename DoFHandler<dim>::active_cell_iterator
            cell2 = dofHandler_Scalar.begin_active(),
            endc2 = dofHandler_Scalar.end(), dg_cell2 = dofHandler_Scalar.begin_active();
        pcout << "\nDefined cell iterator 2" << std::endl;
        unsigned int cellID2 = 0;
        for (; cell2 != endc2; ++cell2) {
          fe_values.reinit(cell2);
          fe_values.get_function_values(solution_cp_2, solution_values_cp_2);
          for (unsigned int q = 0; q < quadrature.size(); ++q) {
            unsigned int i = Utilities::MPI::this_mpi_process(mpi_communicator);
            // std::cout << "Processor " << i
            //           << ", solution_values_cp_2=" << solution_values_cp_2[q]
            //           << std::endl;
            dtwinfraction_iter1[cellID2][q][0] =
                solution_values_cp_2[q]; // Passing the solution for twin of index
                                        // = 0
          }
          if (cell2->is_locally_owned()) {
            cellID2++;
          }
        }
        pcout << "\nInterpolation of dndt complete" << std::endl;
        // ++++++ End of Interpolation of "dndt" ++++++
      }
      
      if (userInputs_cp.enableIndentationBCs) {
        MultiPhysicsBVP<dim, degree>::updateBeforeIncrement();
        if (!userInputs_cp.continuum_Isotropic)
          updateBeforeIncrement();
      } else {
        updateBeforeIncrement();
      }

      if (!userInputs_cp.flagTaylorModel) {
        // solve time increment
        success = solveNonLinearSystem();
      }

      // call updateAfterIncrement, if any
      if ((success) || (userInputs_cp.flagTaylorModel)) {
        updateAfterIncrement();

        if (currentIncrement_cp*delT >=  timeBeforeS){
          // Interpolate twin fraction and twinfraction change from CPFE mesh into
          // PF mesh.

          // ***** Interpolation of twin energy ******
          //
          // Define object fe_function_3 with CPFE vector
          pcout << "\n Passing twin energy from CPFE into Phase Field " << std::endl; 
          Functions::FEFieldFunction<dim,vectorType_cp> fe_function_3(
              dofHandler_Scalar, *postFieldsWithGhosts[3]);
          pcout << "\nCreated fe_function_3 object " << std::endl;
          // Interpolate into the PF domain
          VectorTools::interpolate(*pf_obj.getDofHandlersSet()[2], fe_function_3,
                                  *pf_obj.getSolutionSet()[2]); // Works but not in parallel
          pcout << "\n Interpolation of twin driving force complete " << std::endl;
          // ***** End of interpolation of twin energy ******
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
        std::vector<double> tabularTimeInputInc;
        if (userInputs_cp.tabularOutput) {
          tabularTimeInputInc = userInputs_cp.tabularTimeOutput;
          for (unsigned int i = 0; i < userInputs_cp.tabularTimeOutput.size();
               i++) {
            tabularTimeInputInc[i] = tabularTimeInputInc[i] / delT;
          }
          tabularTimeInputIncInt.resize(userInputs_cp.tabularTimeOutput.size(),
                                        0);
          /// Converting to an integer always rounds down, even if the fraction
          /// part is 0.99999999.
          // Hence, I add 0.1 to make sure we always get the correct integer.
          for (unsigned int i = 0; i < userInputs_cp.tabularTimeOutput.size();
               i++) {
            tabularTimeInputIncInt[i] = int(tabularTimeInputInc[i] + 0.1);
          }
        }
        //////////////////////TabularOutput Finish///////////////
        if (((!userInputs_cp.tabularOutput) &&
             ((currentIncrement_cp + 1) % userInputs_cp.skipOutputSteps ==
              0)) ||
            ((userInputs_cp.tabularOutput) &&
             (std::count(tabularTimeInputIncInt.begin(),
                         tabularTimeInputIncInt.end(),
                         (currentIncrement_cp + 1)) == 1))) {
          if (userInputs_cp.writeOutput)
            output();
        }
        computing_timer_cp.leave_subsection("postprocess");
      } else {
        successiveIncs = 0;
      }
      if (currentIncrement_cp*delT >=  timeBeforeS){
        for (unsigned int pf_step = 0; pf_step < userInputs_pf.increments_pftocpfe; pf_step++){
          //Phase-Field regular step STARTS
          //increment current time
          pf_obj.getCurrentTime_pf() += userInputs_pf.dtValue;
          if (pf_obj.getCurrentIncrement_pf()%userInputs_pf.skip_print_steps==0){
              pcout << "\ntime increment PF:" << pf_obj.getCurrentIncrement_pf() << "  time: " << pf_obj.getCurrentTime_pf() << "\n";
              pcout << "\ncurrent output PF:" << pf_obj.getCurrentOutput() << "\n";
          }

          //check and perform adaptive mesh refinement
          //adaptiveRefine(currentIncrement_pf);

          // Update the list of nuclei (if relevant)
          //updateNucleiList();

          // If grain reassignment is activated, reassign grains
          //if (userInputs_pf.grain_remapping_activated and (currentIncrement_pf%userInputs_pf.skip_grain_reassignment_steps == 0 or currentIncrement_pf == 0) ) {
          //    reassignGrains();
          //}

          //solve time increment
          //pcout << "\n PF time-stepping temporarily disabled"  << std::endl;
          pf_obj.getSolveIncrement(false);

          //if (userInputs_pf.outputTimeStepList[pf_obj.getCurrentOutput()] == pf_obj.getCurrentIncrement_pf()) {
          //Apply constraints and update ghost values
          for(unsigned int fieldIndex=0; fieldIndex<pf_obj.fields.size(); fieldIndex++){     
            pf_obj.getConstraintsDirichletSet()[fieldIndex]->distribute(*pf_obj.getSolutionSet()[fieldIndex]);
            pf_obj.getConstraintsOtherSet()[fieldIndex]->distribute(*pf_obj.getSolutionSet()[fieldIndex]);
            pf_obj.getSolutionSet()[fieldIndex]->update_ghost_values();
          }

          //Output Results
          pf_obj.getOutputResults();
          /*
          if (userInputs_pf.print_timing_with_output && currentIncrement_pf < userInputs_pf.totalIncrements){
              computing_timer_pf.print_summary();
          }
          */
          //currentOutput++;
          pf_obj.getCurrentOutput() += 1;
          //}
          /*
          // Create a checkpoint (on the proper increments)
          if (userInputs_pf.checkpointTimeStepList[currentCheckpoint] == currentIncrement_pf) {
              save_checkpoint();
              currentCheckpoint++;
          }
          */
          //currentIncrement_pf++;
          pf_obj.getCurrentIncrement_pf() += 1;
          //Phase-Field regular step ENDS
        }
      }
    }
}
#include "../../include/multiPhysicsBVP_template_instantiations.h"
