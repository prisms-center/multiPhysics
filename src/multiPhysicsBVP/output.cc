//output method for multiPhysicsBVP class
#include "../../include/multiPhysicsBVP.h"
#include <fstream>

//output results
template <int dim, int degree>
void MultiPhysicsBVP<dim,degree>::output(){
  DataOut<dim> data_out, data_out_Scalar;
  data_out.attach_dof_handler (dofHandler);
  data_out_Scalar.attach_dof_handler (dofHandler_Scalar);

  //add displacement field
  data_out.add_data_vector (solutionWithGhosts,
			    nodal_solution_names,
			    DataOut<dim>::type_dof_data,
			    nodal_data_component_interpretation);

  //add postprocessing fields
  unsigned int numPostProcessedFieldsWritten=0;
  for (unsigned int field=0; field<numPostProcessedFields; field++){
    if(!userInputs_cp.output_Eqv_strain)
        if (postprocessed_solution_names[field].compare(std::string("Eqv_strain"))==0) continue;
    if(!userInputs_cp.output_Eqv_stress)
        if (postprocessed_solution_names[field].compare(std::string("Eqv_stress"))==0) continue;
    if (userInputs_cp.continuum_Isotropic){
      if(!userInputs_cp.output_alpha)
        if (postprocessed_solution_names[field].compare(std::string("alpha"))==0) continue;
    }
    else {
      if(!userInputs_cp.output_Twin)
        if (postprocessed_solution_names[field].compare(std::string("Twin"))==0) continue;
    }


//pcout<<"field="<<field<<"step1\n";

  if(!userInputs_cp.output_Var1)
    if (postprocessed_solution_names[field].compare(std::string("output_Var1"))==0) continue;
  if(!userInputs_cp.output_Var2)
    if (postprocessed_solution_names[field].compare(std::string("output_Var2"))==0) continue;
  if(!userInputs_cp.output_Var3)
    if (postprocessed_solution_names[field].compare(std::string("output_Var3"))==0) continue;
  if(!userInputs_cp.output_Var4)
    if (postprocessed_solution_names[field].compare(std::string("output_Var4"))==0) continue;
  if(!userInputs_cp.output_Var5)
    if (postprocessed_solution_names[field].compare(std::string("output_Var5"))==0) continue;
  if(!userInputs_cp.output_Var6)
    if (postprocessed_solution_names[field].compare(std::string("output_Var6"))==0) continue;
  if(!userInputs_cp.output_Var7)
    if (postprocessed_solution_names[field].compare(std::string("output_Var7"))==0) continue;
  if(!userInputs_cp.output_Var8)
    if (postprocessed_solution_names[field].compare(std::string("output_Var8"))==0) continue;
  if(!userInputs_cp.output_Var9)
    if (postprocessed_solution_names[field].compare(std::string("output_Var9"))==0) continue;
  if(!userInputs_cp.output_Var10)
    if (postprocessed_solution_names[field].compare(std::string("output_Var10"))==0) continue;
  if(!userInputs_cp.output_Var11)
    if (postprocessed_solution_names[field].compare(std::string("output_Var11"))==0) continue;
  if(!userInputs_cp.output_Var12)
    if (postprocessed_solution_names[field].compare(std::string("output_Var12"))==0) continue;
  if(!userInputs_cp.output_Var13)
    if (postprocessed_solution_names[field].compare(std::string("output_Var13"))==0) continue;
  if(!userInputs_cp.output_Var14)
    if (postprocessed_solution_names[field].compare(std::string("output_Var14"))==0) continue;
  if(!userInputs_cp.output_Var15)
    if (postprocessed_solution_names[field].compare(std::string("output_Var15"))==0) continue;
  if(!userInputs_cp.output_Var16)
    if (postprocessed_solution_names[field].compare(std::string("output_Var16"))==0) continue;
  if(!userInputs_cp.output_Var17)
    if (postprocessed_solution_names[field].compare(std::string("output_Var17"))==0) continue;
  if(!userInputs_cp.output_Var18)
    if (postprocessed_solution_names[field].compare(std::string("output_Var18"))==0) continue;
  if(!userInputs_cp.output_Var19)
    if (postprocessed_solution_names[field].compare(std::string("output_Var19"))==0) continue;
  if(!userInputs_cp.output_Var20)
    if (postprocessed_solution_names[field].compare(std::string("output_Var20"))==0) continue;
  if(!userInputs_cp.output_Var21)
    if (postprocessed_solution_names[field].compare(std::string("output_Var21"))==0) continue;
  if(!userInputs_cp.output_Var22)
    if (postprocessed_solution_names[field].compare(std::string("output_Var22"))==0) continue;
  if(!userInputs_cp.output_Var23)
    if (postprocessed_solution_names[field].compare(std::string("output_Var23"))==0) continue;
  if(!userInputs_cp.output_Var24)
    if (postprocessed_solution_names[field].compare(std::string("output_Var24"))==0) continue;




  //if(!userInputs_cp.output_alpha)
    //if (postprocessed_solution_names[field].compare(std::string("alpha"))==0) continue;
  //if(!userInputs_cp.output_tau_vm)
    //if (postprocessed_solution_names[field].compare(std::string("tau_vm"))==0) continue;
    //
    data_out_Scalar.add_data_vector (*postFieldsWithGhosts[field],
				     postprocessed_solution_names[field].c_str());
    numPostProcessedFieldsWritten++;
  }


  //add material id to output file
  Vector<float> material (triangulation_cp.n_active_cells());
  unsigned int matID=0;
  unsigned int cellID=0;
  typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = triangulation_cp.begin_active(), endc = triangulation_cp.end();
  if(userInputs_cp.readExternalMesh){
    for (; cell!=endc; ++cell){
      material(matID) = cell->material_id(); matID++;}
  }
  else{
    for (; cell!=endc; ++cell){
      if(cell->is_locally_owned()){
        material(cell->active_cell_index())=postprocessValuesAtCellCenters(cellID,0);
        cellID++;
      }
    }
  }

  data_out.add_data_vector (material, "meshGrain_ID");
//cp_quad_trial_2
//cp_trial2

Vector<double> distributed_values_per_cell;//(triangulation.n_active_cells());
Vector<double> global_values_per_cell;//(triangulation.n_active_cells());
//vectorType global_values_per_cell;
 QGauss<dim>  quadrature(userInputs_cp.quadOrder);
 FEValues<dim> fe_values_2 (FE_Scalar, quadrature, update_values | update_gradients | update_JxW_values);
  const unsigned int   dofs_per_cell_2   = FE_Scalar.dofs_per_cell;
 //   const unsigned int   dofs_per_cell   = FE_Scalar.n_dofs();
//  locally_owned_dofs_Scalar = dofHandler_Scalar.locally_owned_dofs ();
 // DoFTools::extract_locally_relevant_dofs (dofHandler_Scalar, locally_relevant_dofs_Scalar);
//FullMatrix<double> qpoint_to_dof_matrix (triangulation.n_active_cells(),
//                                         quadrature.size());
 FullMatrix<double> qpoint_to_dof_matrix (dofs_per_cell_2,
                                         quadrature.size());
//const unsigned int   num_quad_points = quadrature.size();
//const unsigned int n_local_cells =triangulation.n_locally_owned_active_cells();
//std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
FETools::compute_projection_from_quadrature_points_matrix
          (FE_Scalar,
           quadrature, quadrature,
           qpoint_to_dof_matrix);
        //    locally_owned_dofs = dofHandler_Scalar.locally_owned_dofs ();
        //    typename DoFHandler<dim>::active_cell_iterator cell = dofHandler.begin_active(), endc = dofHandler.end();
             Vector<double> local_values_at_qpoints(quadrature.size()); 
             local_values_at_qpoints.reinit(quadrature.size());
             //qpoint_to_dof_matrix.reinit(dofs_per_cell,quadrature.size());
            
             distributed_values_per_cell.reinit(FE_Scalar.n_dofs_per_cell());
             global_values_per_cell.reinit(dofHandler_Scalar.n_dofs());
             //for(unsigned int i=0;i<dofHandler_Scalar.n_dofs();i++){
             //global_values_per_cell(i)=0;
            // }
            // global_values_per_cell.reinit (locally_owned_dofs_Scalar, locally_relevant_dofs_Scalar, mpi_communicator); distributed_values_per_cell=0;
//PETScWrappers::MPI::Vector global_values_per_cell(mpi_communicator, triangulation.n_active_cells(), n_local_cells);
typename DoFHandler<dim>::active_cell_iterator cell_t = dofHandler_Scalar.begin_active(), endc_t = dofHandler_Scalar.end(),dg_cell = dofHandler_Scalar.begin_active();
     cellID=0;
    for (; cell_t!=endc_t; ++cell_t,++dg_cell) {
       if (cell_t->is_locally_owned()){
    fe_values_2.reinit (cell_t);
     //cell_t->get_dof_indices (local_dof_indices);
        for (unsigned int q=0; q<quadrature.size(); ++q)
          local_values_at_qpoints[q] =postprocessValues(cellID, q, 3, 0);// 0.5;//stateVar_conv_try[cellID][q];
         //  cell->set_user_index(cellID);
    //      fe_values.reinit (cell_t);
       // cell_t->get_dof_indices (local_dof_indices);
        
        qpoint_to_dof_matrix.vmult (distributed_values_per_cell,
                                    local_values_at_qpoints);
       dg_cell->set_dof_values (distributed_values_per_cell,global_values_per_cell);
    // 	constraintsMassMatrix.distribute_local_to_global(distributed_values_per_cell, local_dof_indices, global_values_per_cell);
  //  postFieldsWithGhosts[4].push_back(distributed_values_per_cell);
     cellID++;
 
      }
      }
         //      postFieldsWithGhosts[4]->compress(VectorOperation::add);
   //   qpoint_to_dof_matrix.vmult (distributed_values_per_cell,
//                                    local_values_at_qpoints);
     //for (unsigned int i = 0; i < distributed_values_per_cell.size(); ++i){
      // if (distributed_values_per_cell(i) != 0)          
     //  global_values_per_cell(i) = distributed_values_per_cell(i);}
      //global_values_per_cell.compress(VectorOperation::add);
      // data_out_Scalar.add_data_vector ( distributed_values_per_cell, "SSD");
     //  data_out_Scalar.add_data_vector ( global_values_per_cell, "SSD");
  //  data_out_Scalar.add_data_vector (*postFieldsWithGhosts[4],
//				     "SSD");
     //   unsigned int n_twin_systems_Size = userInputs.numTwinSystems1;
	//call base class project() function to project post processed fields
//	for (unsigned int field=8; field<8+n_twin_systems_Size;field++){
//	*solution_cp=*postFieldsWithGhosts[field];
//	  Functions::FEFieldFunction<dim,vectorType> fe_function_1 (*dofHandler, *postFieldsWithGhosts[field]);
//	}
//cp_trial2

//cp_quad_trial_2
  Vector<float> FieldsAtCellCenters (triangulation_cp.n_active_cells());
  for (unsigned int field=1; field<numPostProcessedFieldsAtCellCenters; field++){

    if(!userInputs_cp.outputCellCenters_Var1)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var1"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var2)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var2"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var3)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var3"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var4)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var4"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var5)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var5"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var6)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var6"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var7)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var7"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var8)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var8"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var9)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var9"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var10)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var10"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var11)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var11"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var12)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var12"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var13)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var13"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var14)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var14"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var15)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var15"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var16)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var16"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var17)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var17"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var18)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var18"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var19)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var19"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var20)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var20"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var21)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var21"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var22)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var22"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var23)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var23"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var24)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var24"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var25)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var25"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var26)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var26"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var27)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var27"))==0) continue;
    if(!userInputs_cp.outputCellCenters_Var28)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var28"))==0) continue;

    cellID=0;
    cell = triangulation_cp.begin_active(), endc = triangulation_cp.end();
    for (; cell!=endc; ++cell){
      if(cell->is_locally_owned()){
        FieldsAtCellCenters(cell->active_cell_index())=postprocessValuesAtCellCenters(cellID,field);
        cellID++;
      }
    }

    data_out.add_data_vector (FieldsAtCellCenters,
				     postprocessedFieldsAtCellCenters_solution_names[field].c_str());
  }

  //add subdomain id to output file
  Vector<float> subdomain (triangulation_cp.n_active_cells());
  for (unsigned int i=0; i<subdomain.size(); ++i)
    subdomain(i) = triangulation_cp.locally_owned_subdomain();
  data_out.add_data_vector (subdomain, "subdomain");
  if (numPostProcessedFieldsWritten>0){
    data_out_Scalar.add_data_vector (subdomain, "subdomain");
  }

   data_out.build_patches ();

   if (numPostProcessedFieldsWritten>0){
     data_out_Scalar.add_data_vector (material, "meshGrain_ID");
     data_out_Scalar.build_patches ();
   }

  //write to results file
  std::string dir(userInputs_cp.outputDirectory);
  dir+="/";

  //
  unsigned int incrementDigits= (totalIncrements_cp<10000 ? 4 : std::ceil(std::log10(totalIncrements_cp))+1);
  unsigned int domainDigits   = (Utilities::MPI::n_mpi_processes(mpi_communicator)<10000 ? 4 : std::ceil(std::log10(Utilities::MPI::n_mpi_processes(mpi_communicator)))+1);

  const std::string filename = (dir+"solution-" +
				Utilities::int_to_string (currentIncrement_cp,incrementDigits) +
				"." +
				Utilities::int_to_string (triangulation_cp.locally_owned_subdomain(),
							  domainDigits));
  std::ofstream outputFile ((filename + ".vtu").c_str());
  data_out.write_vtu (outputFile);
  //write projected fields, if any
  if (numPostProcessedFieldsWritten>0){
    const std::string filenameForProjectedFields = (dir+"projectedFields-" +
						    Utilities::int_to_string (currentIncrement_cp,incrementDigits) +
						    "." +
						    Utilities::int_to_string (triangulation_cp.locally_owned_subdomain(),
									      domainDigits));
    std::ofstream outputFileForProjectedFields ((filenameForProjectedFields + ".vtu").c_str());
    data_out_Scalar.write_vtu (outputFileForProjectedFields);
  }


  //create pvtu record
  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
    std::vector<std::string> filenames, filenamesForProjectedFields;
    for (unsigned int i=0;i<Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD); ++i){
      filenames.push_back ("solution-" +
			   Utilities::int_to_string (currentIncrement_cp, incrementDigits) +
			   "." +
			   Utilities::int_to_string (i, domainDigits) +
			   + ".vtu");
      if (numPostProcessedFieldsWritten>0){
	filenamesForProjectedFields.push_back ("projectedFields-" +
					       Utilities::int_to_string (currentIncrement_cp, incrementDigits) +
					       "." +
					       Utilities::int_to_string (i, domainDigits) +
					       + ".vtu");
      }
    }
    const std::string filenamepvtu = (dir+"solution-" +
				      Utilities::int_to_string (currentIncrement_cp,incrementDigits) +
				      ".pvtu");
    std::ofstream master_output (filenamepvtu.c_str());
    data_out.write_pvtu_record (master_output, filenames);
    pcout << "CPFE output written to: " << filenamepvtu.c_str();
    //
    if (numPostProcessedFieldsWritten>0){
      const std::string filenamepvtuForProjectedFields = (dir+"projectedFields-" +
							  Utilities::int_to_string (currentIncrement_cp,incrementDigits) +
							  ".pvtu");
      std::ofstream master_outputForProjectedFields (filenamepvtuForProjectedFields.c_str());
      data_out_Scalar.write_pvtu_record (master_outputForProjectedFields, filenamesForProjectedFields);
      pcout << " and " << filenamepvtuForProjectedFields.c_str() ;
    }
    pcout << " \n\n";
  }
}
#include "../../include/multiPhysicsBVP_template_instantiations.h"
