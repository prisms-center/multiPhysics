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

  // Increase the current increment from 0 to 1 now that the initial conditions have been output
  //currentIncrement_pf++;
  pf_obj.getCurrentIncrement_pf() += 1;

  // Cycle up to the proper output counter
  while (userInputs_pf.outputTimeStepList.size() > 0 && userInputs_pf.outputTimeStepList[pf_obj.getCurrentOutput()] < pf_obj.getCurrentIncrement_pf()){
      //currentOutput++;
      pf_obj.getCurrentOutput() += 1;
  }
  
  //time stepping
  pcout << "\nTime stepping parameters: timeStep: " << userInputs_pf.dtValue << "  timeFinal: " << userInputs_pf.finalTime << "  timeIncrements: " << userInputs_pf.totalIncrements_pf << "\n";
  //Section for first phase field step solution ENDS

  pcout << "begin solve... CPFE\n\n";

  bool success;
  // load increments
  unsigned int successiveIncs = 0;

  //Setting up interpolation from PF to CPFE mesh
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

  if (userInputs_cp.enableAdaptiveTimeStepping) {
    for (; totalLoadFactor < totalIncrements_cp;) {
      ++currentIncrement_cp;
      loadFactorSetByModel =
          std::min(loadFactorSetByModel, totalIncrements_cp - totalLoadFactor);
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
    for (; currentIncrement_cp < totalIncrements_cp; ++currentIncrement_cp) {
      pcout << "\nincrement CPFE: " << currentIncrement_cp << ", Time CPFE: "<< currentIncrement_cp*delT << std::endl;

      if ((currentIncrement_cp*delT >=  seedingT)){
        // ***** Interpolation of order parameter "n" from PF mesh into
        // twin volume fraction CPFE mesh ******
        interpolate_order_parameter(pf_obj, dofHandler_Scalar, quadrature, twinfraction_iter1, fe_values);
        pcout << "\nInterpolation of n complete" << std::endl;
        
        pcout << "\nInterpolation of dndt disabled" << std::endl;
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

        if (currentIncrement_cp*delT >=  timeBeforeC){
          // ***** Interpolation of twin energy ******
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
      if (currentIncrement_cp*delT >=  timeBeforeC){
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
          if (userInputs_pf.print_timing_with_output && currentIncrement_pf < userInputs_pf.totalIncrements_pf){
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
