/*
 * initialConditions.h
 *
 *  Created on: Feb 27, 2017
 *      Author: stephendewitt
 */

#ifndef INCLUDE_INITIALCONDITIONS_H_
#define INCLUDE_INITIALCONDITIONS_H_

#include "matrixFreePDE.h"
#include "userInputParameters_pf.h"

template <int dim, int degree>
class InitialCondition : public dealii::Function<dim>
{
public:
  const unsigned int                index;
  const userInputParameters_pf<dim> userInputs;
  dealii::Vector<double>            values;

  InitialCondition(const unsigned int                _index,
                   const userInputParameters_pf<dim> _userInputs,
                   MatrixFreePDE<dim, degree>     *_matrixfreepde)
    : dealii::Function<dim>(1)
    , index(_index)
    , userInputs(_userInputs)
    , matrixfree_pde(_matrixfreepde)
  {
    std::srand(dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) + 1);
  }

  // IC for scalar values
  double
  value(const dealii::Point<dim> &p, const unsigned int component = 0) const
  {
    double                 scalar_IC = 0.0;
    dealii::Vector<double> vector_IC(dim);

    matrixfree_pde->setInitialCondition(p, index, scalar_IC, vector_IC);
    return scalar_IC;
  };

private:
  MatrixFreePDE<dim, degree> *matrixfree_pde;
};

template <int dim, int degree>
class InitialConditionVector : public dealii::Function<dim>
{
public:
  const unsigned int                index;
  const userInputParameters_pf<dim> userInputs;
  dealii::Vector<double>            values;

  InitialConditionVector(const unsigned int                _index,
                         const userInputParameters_pf<dim> _userInputs,
                         MatrixFreePDE<dim, degree>     *_matrixfreepde)
    : dealii::Function<dim>(dim)
    , index(_index)
    , userInputs(_userInputs)
    , matrixfree_pde(_matrixfreepde)
  {
    std::srand(dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) + 1);
  }

  // IC for vector values
  void
  vector_value(const dealii::Point<dim> &p, dealii::Vector<double> &vector_IC) const
  {
    double scalar_IC = 0.0;
    vector_IC.reinit(dim);
    matrixfree_pde->setInitialCondition(p, index, scalar_IC, vector_IC);
  };

private:
  MatrixFreePDE<dim, degree> *matrixfree_pde;
};

#endif /* INCLUDE_INITIALCONDITIONS_H_ */
