//solve_pf() method for MatrixFreePDE class

#include "../../include/matrixFreePDE.h"

//solve BVP
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::solve_pf(){
    //log time
    computing_timer_pf.enter_subsection("matrixFreePDE: solve");
    pcout << "\nsolving...\n\n";

    //time dependent BVP
    if (isTimeDependentBVP){

        // If grain reassignment is activated, reassign grains
        if (userInputs_pf.grain_remapping_activated and (currentIncrement_pf%userInputs_pf.skip_grain_reassignment_steps == 0 or currentIncrement_pf == 0) ) {
            reassignGrains();
        }

        // For any nonlinear equation, set the initial guess as the solution to Laplace's equations
        generatingInitialGuess = true;
        setNonlinearEqInitialGuess();
        generatingInitialGuess = false;

        // Do an initial solve to set the elliptic fields
        solveIncrement(true);

        //output initial conditions for time dependent BVP
        if (userInputs_pf.outputTimeStepList[currentOutput] == currentIncrement_pf) {

            for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
                constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
                constraintsOtherSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
                solutionSet[fieldIndex]->update_ghost_values();
            }
            outputResults();
            currentOutput++;
        }

        if (userInputs_pf.checkpointTimeStepList[currentCheckpoint] == currentIncrement_pf) {
            save_checkpoint();
            currentCheckpoint++;
        }

        // Increase the current increment from 0 to 1 now that the initial conditions have been output
        currentIncrement_pf++;

        // Cycle up to the proper output and checkpoint counters
        while (userInputs_pf.outputTimeStepList.size() > 0 && userInputs_pf.outputTimeStepList[currentOutput] < currentIncrement_pf){
            currentOutput++;
        }
        while (userInputs_pf.checkpointTimeStepList.size() > 0 && userInputs_pf.checkpointTimeStepList[currentCheckpoint] < currentIncrement_pf){
            currentCheckpoint++;
        }

        //time stepping
        pcout << "\nTime stepping parameters: timeStep: " << userInputs_pf.dtValue << "  timeFinal: " << userInputs_pf.finalTime << "  timeIncrements: " << userInputs_pf.totalIncrements << "\n";

        // This is the main time-stepping loop
        for (; currentIncrement_pf<=userInputs_pf.totalIncrements; ++currentIncrement_pf){
            //increment current time
            currentTime_pf+=userInputs_pf.dtValue;
            if (currentIncrement_pf%userInputs_pf.skip_print_steps==0){
                pcout << "\ntime increment:" << currentIncrement_pf << "  time: " << currentTime_pf << "\n";
            }

            //check and perform adaptive mesh refinement
            adaptiveRefine(currentIncrement_pf);

            // Update the list of nuclei (if relevant)
            updateNucleiList();

            // If grain reassignment is activated, reassign grains
            if (userInputs_pf.grain_remapping_activated and (currentIncrement_pf%userInputs_pf.skip_grain_reassignment_steps == 0 or currentIncrement_pf == 0) ) {
                reassignGrains();
            }

            //solve time increment
            solveIncrement(false);

            // Output results to file (on the proper increments)
            if (userInputs_pf.outputTimeStepList[currentOutput] == currentIncrement_pf) {
                for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
                    constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
                    constraintsOtherSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
                    solutionSet[fieldIndex]->update_ghost_values();
                }
                outputResults();
                if (userInputs_pf.print_timing_with_output && currentIncrement_pf < userInputs_pf.totalIncrements){
                    computing_timer_pf.print_summary();
                }

                currentOutput++;
            }

            // Create a checkpoint (on the proper increments)
            if (userInputs_pf.checkpointTimeStepList[currentCheckpoint] == currentIncrement_pf) {
                save_checkpoint();
                currentCheckpoint++;
            }

        }
    }

    //time independent BVP
    else{
        generatingInitialGuess = false;

        //solve
        solveIncrement(false);

        //output results to file
        outputResults();
    }

    //log time
    computing_timer_pf.leave_subsection("matrixFreePDE: solve");
}

#include "../../include/matrixFreePDE_template_instantiations.h"
