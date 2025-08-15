// constructor and destructor for MultiPhysics class
#include "../../include/multiPhysicsBVP.h"

// constructor
template <int dim, int degree>
MultiPhysicsBVP<dim, degree>::MultiPhysicsBVP(userInputParameters_pf<dim> _userInputs_pf,
                                              userInputParameters_cp      _userInputs_cp)
  : Subscriptor()
  , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  , mpi_communicator(MPI_COMM_WORLD)
  , userInputs_cp(_userInputs_cp)
  , userInputs_pf(_userInputs_pf)
  , triangulation_cp(mpi_communicator,
                     typename Triangulation<dim>::MeshSmoothing(
                       Triangulation<dim>::smoothing_on_refinement |
                       Triangulation<dim>::smoothing_on_coarsening))
  , triangulation2_cp(mpi_communicator,
                      typename Triangulation<dim>::MeshSmoothing(
                        Triangulation<dim>::smoothing_on_refinement |
                        Triangulation<dim>::smoothing_on_coarsening))
  , FE(FE_Q<dim>(_userInputs_cp.feOrder), dim)
  , FE_Scalar(FE_Q<dim>(_userInputs_cp.feOrder), 1)
  , dofHandler(triangulation_cp)
  , dofHandler_Scalar(triangulation_cp)
  , delT(_userInputs_cp.delT)
  , totalT(_userInputs_cp.totalTime)
  , seedingT(_userInputs_cp.seedingTime)
  , timeBeforeC(_userInputs_cp.timeBeforeCoupling)
  , delT_pf_adjust(_userInputs_cp.delT_pf_adjust)
  , currentIteration(0)
  , currentIncrement_cp(0)
  , resetIncrement(false)
  , loadFactorSetByModel(1.0)
  , totalLoadFactor(0.0)
  , computing_timer_cp(pcout, TimerOutput::summary, TimerOutput::wall_times)
  , numPostProcessedFields(0)
{
  // Nodal Solution names - this is for writing the output file
  for (unsigned int i = 0; i < dim; ++i)
    {
      nodal_solution_names.push_back("u");
      nodal_data_component_interpretation.push_back(
        DataComponentInterpretation::component_is_part_of_vector);
    }
  if (userInputs_cp.enableCyclicLoading)
    {
      cycleTime = 4 * userInputs_cp.quarterCycleTime;
    }
  totalIncrements_cp = std::round(totalT / delT);
  if (userInputs_cp.enableTabularPeriodicBCs)
    {
      periodicTotalIncrements = userInputs_cp.periodicTabularTime / delT;
    }
}

// destructor
template <int dim, int degree>
MultiPhysicsBVP<dim, degree>::~MultiPhysicsBVP()
{}

#include "../../include/multiPhysicsBVP_template_instantiations.h"