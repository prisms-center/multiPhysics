// constructor and destructor for matrixFreePDE class

#include "matrixFreePDE.h"

// constructor
template <int dim, int degree>
MatrixFreePDE<dim, degree>::MatrixFreePDE(userInputParameters_pf<dim> _userInputs,
                                          userInputParameters_cp & _userInputs_cp)
  : Subscriptor()
  , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  , userInputs(_userInputs)
  , triangulation(MPI_COMM_WORLD)
  , currentFieldIndex(0)
  , isTimeDependentBVP(false)
  , isEllipticBVP(false)
  , hasExplicitEquation(false)
  , hasNonExplicitEquation(false)
  , currentTime(0.0)
  , currentIncrement(0)
  , currentOutput(0)
  , currentCheckpoint(0)
  , current_grain_reassignment(0)
  , computing_timer(pcout, TimerOutput::summary, TimerOutput::wall_times)
  , first_integrated_var_output_complete(false)
{}

// destructor
template <int dim, int degree>
MatrixFreePDE<dim, degree>::~MatrixFreePDE()
{
  matrixFreeObject.clear();

  // Delete the pointers contained in several member variable vectors
  // The size of each of these must be checked individually in case an exception
  // is thrown as they are being initialized.
  for (const auto &locally_relevant_dofs : locally_relevant_dofsSet)
    {
      delete locally_relevant_dofs;
    }
  for (const auto &constraintsDirichlet : constraintsDirichletSet)
    {
      delete constraintsDirichlet;
    }
  for (const auto &soltrans : soltransSet)
    {
      delete soltrans;
    }
  for (const auto &dofHandlers : dofHandlersSet)
    {
      delete dofHandlers;
    }
  for (const auto &FE : FESet)
    {
      delete FE;
    }
  for (const auto &solution : solutionSet)
    {
      delete solution;
    }
  for (const auto &residual : residualSet)
    {
      delete residual;
    }
}

template<int dim, int degree>
std::vector<const dealii::DoFHandler<dim>*>& MatrixFreePDE<dim,degree>::getDofHandlersSet()
{
  return dofHandlersSet;
}

template<int dim, int degree>
std::vector<vectorType_pf*>& MatrixFreePDE<dim,degree>::getSolutionSet()
{
    return solutionSet;
}

template<int dim, int degree>
std::vector<const dealii::AffineConstraints<double>*>& MatrixFreePDE<dim,degree>::getConstraintsDirichletSet()
{
    return constraintsDirichletSet;
}

template<int dim, int degree>
std::vector<const dealii::AffineConstraints<double>*>& MatrixFreePDE<dim,degree>::getConstraintsOtherSet()
{
    return constraintsOtherSet;
}

template<int dim, int degree>
void MatrixFreePDE<dim,degree>::getOutputResults()
{
    return outputResults();
}

template<int dim, int degree>
double& MatrixFreePDE<dim,degree>::getCurrentTime()
{
    return currentTime;
}

template<int dim, int degree>
unsigned int& MatrixFreePDE<dim,degree>::getCurrentIncrement()
{
    return currentIncrement;
}

template<int dim, int degree>
unsigned int& MatrixFreePDE<dim,degree>::getCurrentOutput()
{
    return currentOutput;
}

#include "../../include/matrixFreePDE_template_instantiations.h"