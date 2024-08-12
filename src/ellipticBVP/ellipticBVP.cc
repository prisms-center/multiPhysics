//constructor and destructor for ellipticBVP class
#include "../../include/multiPhysicsBVP.h"

//constructor
template <int dim>
//ellipticBVP<dim>::ellipticBVP (userInputParameters _userInputs_cp):multiPhysicsBVP_templamultiPhysicsBVP_template_instatiations.h
MultiPhysicsBVP<dim>:: MultiPhysicsBVP(userInputParameters_pf<dim>, userInputParameters_cp _userInputs_cp);
  Subscriptor(),
  mpi_communicator (MPI_COMM_WORLD),
  userInputs_cp(_userInputs_cp),
  triangulation_cp (mpi_communicator,
		 typename Triangulation<dim>::MeshSmoothing
		 (Triangulation<dim>::smoothing_on_refinement |
		  Triangulation<dim>::smoothing_on_coarsening)),
  triangulation2_cp (mpi_communicator,
		 typename Triangulation<dim>::MeshSmoothing
		 (Triangulation<dim>::smoothing_on_refinement |
		  Triangulation<dim>::smoothing_on_coarsening)),
  FE (FE_Q<dim>(_userInputs_cp.feOrder), dim),
  FE_Scalar (FE_Q<dim>(_userInputs_cp.feOrder), 1),
  dofHandler (triangulation_cp),
  dofHandler_Scalar (triangulation_cp),
  delT(_userInputs_cp.delT),
  totalT(_userInputs_cp.totalTime),
  currentIteration(0),
  currentIncrement_cp(0),
  resetIncrement(false),
  loadFactorSetByModel(1.0),
  totalLoadFactor(0.0),
  pcout (std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0),
  computing_timer_cp (pcout, TimerOutput::summary, TimerOutput::wall_times),
  numPostProcessedFields(0)
{
  //Nodal Solution names - this is for writing the output file
  for (unsigned int i=0; i<dim; ++i){
    nodal_solution_names.push_back("u");
    nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  }
  if(userInputs_cp.enableCyclicLoading)
    cycleTime=4*userInputs_cp.quarterCycleTime;
  totalIncrements=std::round(totalT/delT);
  if(userInputs_cp.enableTabularPeriodicBCs)
    periodicTotalIncrements=userInputs_cp.periodicTabularTime/delT;
}

//destructor
template <int dim>
MultiPhysicsBVP<dim>::~MultiPhysicsBVP ()
{
}

#include "../../include/multiPhysicsBVP_template_instatiations.h"
