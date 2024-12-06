//constructor and destructor for MultiPhysics class
#include "../../include/multiPhysicsBVP.h"

//constructor
template <int dim, int degree>
MultiPhysicsBVP<dim,degree>:: MultiPhysicsBVP(userInputParameters_pf<dim> _userInputs_pf, userInputParameters_cp _userInputs_cp):
  Subscriptor(),
  pcout (std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0),
  mpi_communicator (MPI_COMM_WORLD),
  //Plasticity variables
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
  timeBeforeS(_userInputs_cp.timeBeforeSeeding),
  currentIteration(0),
  currentIncrement_cp(0),
  resetIncrement(false),
  loadFactorSetByModel(1.0),
  totalLoadFactor(0.0),
  computing_timer_cp (pcout, TimerOutput::summary, TimerOutput::wall_times),
  numPostProcessedFields(0),
  //Phase-Field variables
  userInputs_pf(_userInputs_pf),
  triangulation_pf (MPI_COMM_WORLD),
  currentFieldIndex(0),
  isTimeDependentBVP(false),
  isEllipticBVP(false),
  hasExplicitEquation(false),
  hasNonExplicitEquation(false),
  parabolicFieldIndex(0),
  ellipticFieldIndex(0),
  currentTime_pf(0.0),
  currentIncrement_pf(0),
  currentOutput(0),
  currentCheckpoint(0),
  current_grain_reassignment(0),
  computing_timer_pf (pcout, TimerOutput::summary, TimerOutput::wall_times),
  first_integrated_var_output_complete(false)
{
  //Nodal Solution names - this is for writing the output file
  for (unsigned int i=0; i<dim; ++i){
    nodal_solution_names.push_back("u");
    nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  }
  if(userInputs_cp.enableCyclicLoading)
    cycleTime=4*userInputs_cp.quarterCycleTime;
  totalIncrements_cp=std::round(totalT/delT);
  if(userInputs_cp.enableTabularPeriodicBCs)
    periodicTotalIncrements=userInputs_cp.periodicTabularTime/delT;
}

//destructor
template <int dim, int degree>
MultiPhysicsBVP<dim,degree>::~MultiPhysicsBVP ()
{
     matrixFreeObject.clear();

   // Delete the pointers contained in several member variable vectors
   // The size of each of these must be checked individually in case an exception is thrown
   // as they are being initialized.
   for(unsigned int iter=0; iter<locally_relevant_dofsSet.size(); iter++){
       delete locally_relevant_dofsSet[iter];
   }
   for(unsigned int iter=0; iter<constraintsDirichletSet.size(); iter++){
       delete constraintsDirichletSet[iter];
   }
   for(unsigned int iter=0; iter<soltransSet.size(); iter++){
       delete soltransSet[iter];
   }
   for(unsigned int iter=0; iter<dofHandlersSet.size(); iter++){
       delete dofHandlersSet[iter];
   }
   for(unsigned int iter=0; iter<FESet.size(); iter++){
       delete FESet[iter];
   }
   for(unsigned int iter=0; iter<solutionSet.size(); iter++){
       delete solutionSet[iter];
   }
   for(unsigned int iter=0; iter<residualSet.size(); iter++){
       delete residualSet[iter];
   }
}

#include "../../include/multiPhysicsBVP_template_instantiations.h"