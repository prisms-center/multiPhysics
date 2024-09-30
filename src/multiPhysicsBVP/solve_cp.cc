//solve method for multiPhysicsBVP class
#include "../../include/multiPhysicsBVP.h"
#include "customPDE.h"

#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/fe/fe_tools.h>

//loop over increments and solve each increment
template <int dim, int degree>
void MultiPhysicsBVP<dim,degree>::solve_cp(){
  pcout << "begin solve...\n\n";
  bool success;
  //load increments
  unsigned int successiveIncs=0;
//cp_twinfraction_addition
	 //local variables
  QGauss<dim>  quadrature(userInputs_cp.quadOrder);
  FEValues<dim> fe_values (FE, quadrature, update_values | update_gradients | update_JxW_values);
  const unsigned int   dofs_per_cell   = FE.dofs_per_cell;
  const unsigned int   num_quad_points = quadrature.size();
  unsigned int num_local_cells = triangulation_cp.n_locally_owned_active_cells();
   std::vector<double> twin_init(userInputs_cp.numTwinSystems1);
    for (unsigned int i=0;i<userInputs_cp.numTwinSystems1;i++){
    twin_init[i]=0.0;
      }
  twinfraction_iter1.resize(num_local_cells,std::vector<std::vector<double> >(num_quad_points,twin_init));
  dtwinfraction_iter1.resize(num_local_cells,std::vector<std::vector<double> >(num_quad_points,twin_init));
	//
  if(userInputs_cp.enableAdaptiveTimeStepping){
    for (;totalLoadFactor<totalIncrements;){
      ++currentIncrement_cp;
      loadFactorSetByModel=std::min(loadFactorSetByModel, totalIncrements-totalLoadFactor);
      pcout << "\nincrement: "  << currentIncrement_cp << std::endl;
      char buffer[100];
      sprintf(buffer, "current load factor: %12.6e\ntotal load factor:   %12.6e\n", loadFactorSetByModel, totalLoadFactor);
      pcout << buffer;

      //call updateBeforeIncrement, if any
      updateBeforeIncrement();

      if (!userInputs_cp.flagTaylorModel){
        //solve time increment
        success=solveNonLinearSystem();
      }

      //call updateAfterIncrement, if any
      if ((success)||(userInputs_cp.flagTaylorModel)){
        updateAfterIncrement();

        //update totalLoadFactor
        totalLoadFactor+=loadFactorSetByModel;

        //increase loadFactorSetByModel, if succesiveIncForIncreasingTimeStep satisfied.
        successiveIncs++;

        if (successiveIncs>=userInputs_cp.succesiveIncForIncreasingTimeStep){
          loadFactorSetByModel*=userInputs_cp.adaptiveLoadIncreaseFactor;
          char buffer1[100];
          sprintf(buffer1, "current increment increased. Restarting increment with loadFactorSetByModel: %12.6e\n", loadFactorSetByModel);
          pcout << buffer1;
        }
        #if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 3)&&(DEAL_II_VERSION_MAJOR==9)))
        computing_timer_cp.enter_section("postprocess");
        #else
	    computing_timer_cp.enter_subsection("postprocess");
        #endif

        if (currentIncrement_cp%userInputs_cp.skipOutputSteps==0)
        if (userInputs_cp.writeOutput) output();

        #if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 3)&&(DEAL_II_VERSION_MAJOR==9)))
        computing_timer_cp.exit_section("postprocess");
        #else
	    computing_timer_cp.leave_subsection("postprocess");
        #endif

      }
      else{
        successiveIncs=0;
      }
    }
    char buffer[100];
    sprintf(buffer, "\nfinal load factor  : %12.6e\n", totalLoadFactor);
    pcout << buffer;
  }
  else
    for (;currentIncrement_cp<totalIncrements; ++currentIncrement_cp){
      pcout << "\nincrement: "  << currentIncrement_cp << std::endl;

      //Interpolate twin fraction and twinfraction change from PF mesh into CPFE mesh 
      //Initial interpolation of twin seed (Order parameter "n", solutionSet[0]) in phase-field mesh into variable dtwinfraction_iter in CPFE mesh
      
      //Define object fe_function_1 for phase field data given dof handler and the given solution vector. . This object can evaluate the finite 
      //element solution at given points in the domain.
      // Accessing pf_object through the virtual function
      auto& pf_obj = this->get_pf_object();
      pcout << "\npf_obj access successful " << std::endl;
      Functions::FEFieldFunction<dim, vectorType_pf> fe_function_1(*pf_obj.getDofHandlersSet()[0], *pf_obj.getSolutionSet()[0]);
      pcout << "\nCreated fe_function_1 object " << std::endl;
      //Interpolate into the CPFE domain
      vectorType_cp solution_cp;  // No pointer, just a regular object
      IndexSet own_dofs = dofHandler.locally_owned_dofs();
      DoFTools::extract_locally_relevant_dofs (dofHandler, locally_relevant_dofs);
      solution_cp.reinit(own_dofs, locally_relevant_dofs ,mpi_communicator);
      VectorTools::interpolate (dofHandler, fe_function_1, solution_cp); //Works but not in parallel
      pcout << "\nInterpolated into solution_cp " << std::endl;
      //Define dof_to_qpoint_matrix 
      FullMatrix<double> dof_to_qpoint_matrix (quadrature.size(),FE.n_dofs_per_cell());
      pcout << "\nDeclared dof_to_qpoint_matrix " << std::endl;
      //Interpolate into CPFE mesh quadrature points
      FETools::compute_interpolation_to_quadrature_points_matrix(FE,quadrature,dof_to_qpoint_matrix);
      pcout << "\nCalculated dof_to_qpoint_matrix " << std::endl;
      //FEValues<dim> fe_values_cp (*FE, quadrature, update_values | update_gradients | update_JxW_values | update_quadrature_points);
      typename DoFHandler<dim>::active_cell_iterator cell = dofHandler.begin_active(),
                                                  endc = dofHandler.end(),
                                                  dg_cell = dofHandler.begin_active();
      pcout << "\nDefined cell iterator" << std::endl;
      //unsigned int num_local_cells = triangulation_cp.n_locally_owned_active_cells();
      //std::vector<std::vector<double>> twinfraction_iter1; //This should be dtwinfraction_iter
      //twinfraction_iter1.resize(num_local_cells,std::vector<double>(quadrature_cp.size(),0));
      std::vector<double> solution_values_cp(quadrature.size());
      pcout << "\nDefined solution_values_cp for the cell" << std::endl;
      unsigned int cellID=0;
      for (; cell != endc; ++cell){
        fe_values.reinit(cell);
        fe_values.get_function_values(solution_cp, solution_values_cp);
        for (unsigned int q=0; q<quadrature.size(); ++q){
          twinfraction_iter1[cellID][q][0]= solution_values_cp[q]; //Passing the solution for twin of index =0
        }
        if (cell->is_locally_owned()){cellID++;}
      } 
      //End of Initial Interpolation
      pcout << "\nInterpolation complete" << std::endl;

    if (userInputs_cp.enableIndentationBCs){
        MultiPhysicsBVP<dim,degree>::updateBeforeIncrement();
        if (!userInputs_cp.continuum_Isotropic)
            updateBeforeIncrement();
    }
    else{
        updateBeforeIncrement();
    }

    if (!userInputs_cp.flagTaylorModel){
      //solve time increment
      success=solveNonLinearSystem();
    }

    //call updateAfterIncrement, if any
    if ((success)||(userInputs_cp.flagTaylorModel)){
      updateAfterIncrement();

      //update totalLoadFactor
      totalLoadFactor+=loadFactorSetByModel;

      //increase loadFactorSetByModel, if succesiveIncForIncreasingTimeStep satisfied.
      successiveIncs++;
      //output results to file
      //
      #if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 3)&&(DEAL_II_VERSION_MAJOR==9)))
      computing_timer_cp.enter_section("postprocess");
      #else
      computing_timer_cp.enter_subsection("postprocess");
      #endif

      //////////////////////TabularOutput Start///////////////
      std::vector<unsigned int> tabularTimeInputIncInt;
      std::vector<double> tabularTimeInputInc;
      if (userInputs_cp.tabularOutput){

        tabularTimeInputInc=userInputs_cp.tabularTimeOutput;
        for(unsigned int i=0;i<userInputs_cp.tabularTimeOutput.size();i++){
          tabularTimeInputInc[i]=tabularTimeInputInc[i]/delT;
        }
        tabularTimeInputIncInt.resize(userInputs_cp.tabularTimeOutput.size(),0);
        ///Converting to an integer always rounds down, even if the fraction part is 0.99999999.
        //Hence, I add 0.1 to make sure we always get the correct integer.
        for(unsigned int i=0;i<userInputs_cp.tabularTimeOutput.size();i++){
          tabularTimeInputIncInt[i]=int(tabularTimeInputInc[i]+0.1);
        }
      }
      //////////////////////TabularOutput Finish///////////////
      if (((!userInputs_cp.tabularOutput)&&((currentIncrement_cp+1)%userInputs_cp.skipOutputSteps==0))||((userInputs_cp.tabularOutput)&& (std::count(tabularTimeInputIncInt.begin(), tabularTimeInputIncInt.end(), (currentIncrement_cp+1))==1))){
        if (userInputs_cp.writeOutput) output();
      }

      #if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 3)&&(DEAL_II_VERSION_MAJOR==9)))
      computing_timer_cp.exit_section("postprocess");
      #else
      computing_timer_cp.leave_subsection("postprocess");
      #endif
    }
    else{
      successiveIncs=0;
    }
  }
}
#include "../../include/multiPhysicsBVP_template_instantiations.h"
