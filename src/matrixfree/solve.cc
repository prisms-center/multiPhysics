//solve() method for MatrixFreePDE class

#include "../../include/matrixFreePDE.h"
#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/fe/fe_tools.h>

//solve BVP
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::solve(){
    //log time
    computing_timer.enter_subsection("matrixFreePDE: solve");
    pcout << "\nsolving...\n\n";

    //time dependent BVP
    if (isTimeDependentBVP){

        // If grain reassignment is activated, reassign grains
        if (userInputs.grain_remapping_activated and (currentIncrement%userInputs.skip_grain_reassignment_steps == 0 or currentIncrement == 0) ) {
            reassignGrains();
        }

        // For any nonlinear equation, set the initial guess as the solution to Laplace's equations
        generatingInitialGuess = true;
        setNonlinearEqInitialGuess();
        generatingInitialGuess = false;

        // Do an initial solve to set the elliptic fields
        solveIncrement(true);

        //For PRISMS-MP interpolate solution vector for Into PRISMS-Plasticity mesh

        //output initial conditions for time dependent BVP
        if (userInputs.outputTimeStepList[currentOutput] == currentIncrement) {

            for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
                constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
                constraintsOtherSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
                solutionSet[fieldIndex]->update_ghost_values();
            }

            outputResults();

            //For PRISMS-MP interpolate solution vector for Into PRISMS-Plasticity mesh
            Functions::FEFieldFunction<dim,vectorType> fe_function_1 (*dofHandlersSet[0], *solutionSet[0]);
            //Functions::FEFieldFunction<dim,vectorType> solution_function (*dof_handler_cp, *solution_cp);
            VectorTools::interpolate (*dof_handler_cp, fe_function_1, *solution_cp);
            QGauss<dim> quadrature_cp (degree+1);
            pcout << "n_dofs_per_cell: " << fe_cp->n_dofs_per_cell() << std::endl;
            FullMatrix<double> dof_to_qpoint_matrix (quadrature_cp.size(), fe_cp->n_dofs_per_cell());
            FETools::compute_interpolation_to_quadrature_points_matrix(*fe_cp,quadrature_cp,dof_to_qpoint_matrix);
            FEValues<dim> fe_values_cp (*fe_cp, quadrature_cp, update_values | update_gradients | update_JxW_values | update_quadrature_points);
            
            typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_cp->begin_active(),
                                               endc = dof_handler_cp->end(),
                                               dg_cell = dof_handler_cp->begin_active();
                                               
            unsigned int num_local_cells =triangulation_cp.n_locally_owned_active_cells();
            std::vector<std::vector<double>>   twinfraction_iter1;
            twinfraction_iter1.resize(num_local_cells,std::vector<double>(quadrature_cp.size(),0));
            std::vector<double> solution_values_cp(quadrature_cp.size());
            unsigned int cellID=0;
            pcout << "Number of Local Cells (CP mesh): " << num_local_cells << "\n";
            for (; cell != endc; ++cell){
                fe_values_cp.reinit(cell);
                fe_values_cp.get_function_values(*solution_cp, solution_values_cp);
                for (unsigned int q=0; q<quadrature_cp.size(); ++q){
                    twinfraction_iter1[cellID][q]= solution_values_cp[q];
                }
            if (cell->is_locally_owned()){cellID++;}
            }
            outputResults_cp();
            currentOutput++;
        }

        if (userInputs.checkpointTimeStepList[currentCheckpoint] == currentIncrement) {
            save_checkpoint();
            currentCheckpoint++;
        }

        // Increase the current increment from 0 to 1 now that the initial conditions have been output
        currentIncrement++;

        // Cycle up to the proper output and checkpoint counters
        while (userInputs.outputTimeStepList.size() > 0 && userInputs.outputTimeStepList[currentOutput] < currentIncrement){
            currentOutput++;
        }
        while (userInputs.checkpointTimeStepList.size() > 0 && userInputs.checkpointTimeStepList[currentCheckpoint] < currentIncrement){
            currentCheckpoint++;
        }

        //time stepping
        pcout << "\nTime stepping parameters: timeStep: " << userInputs.dtValue << "  timeFinal: " << userInputs.finalTime << "  timeIncrements: " << userInputs.totalIncrements << "\n";

        // This is the main time-stepping loop
        for (; currentIncrement<=userInputs.totalIncrements; ++currentIncrement){
            //increment current time
            currentTime+=userInputs.dtValue;
            if (currentIncrement%userInputs.skip_print_steps==0){
                pcout << "\ntime increment:" << currentIncrement << "  time: " << currentTime << "\n";
            }

            //check and perform adaptive mesh refinement
            adaptiveRefine(currentIncrement);

            // Update the list of nuclei (if relevant)
            updateNucleiList();

            // If grain reassignment is activated, reassign grains
            if (userInputs.grain_remapping_activated and (currentIncrement%userInputs.skip_grain_reassignment_steps == 0 or currentIncrement == 0) ) {
                reassignGrains();
            }

            //solve time increment
            solveIncrement(false);

            // Output results to file (on the proper increments)
            if (userInputs.outputTimeStepList[currentOutput] == currentIncrement) {
                for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
                    constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
                    constraintsOtherSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
                    solutionSet[fieldIndex]->update_ghost_values();
                }
                outputResults();
                
                //For PRISMS-MP interpolate solution vector for Into PRISMS-Plasticity mesh
                Functions::FEFieldFunction<dim,vectorType> fe_function_1 (*dofHandlersSet[0], *solutionSet[0]);
                //Functions::FEFieldFunction<dim,vectorType> solution_function (*dof_handler_cp, *solution_cp);
                VectorTools::interpolate (*dof_handler_cp, fe_function_1, *solution_cp);

            outputResults_cp();
            currentOutput++;
                if (userInputs.print_timing_with_output && currentIncrement < userInputs.totalIncrements){
                    computing_timer.print_summary();
                }

                currentOutput++;
            }

            // Create a checkpoint (on the proper increments)
            if (userInputs.checkpointTimeStepList[currentCheckpoint] == currentIncrement) {
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
    computing_timer.leave_subsection("matrixFreePDE: solve");
}

#include "../../include/matrixFreePDE_template_instantiations.h"
