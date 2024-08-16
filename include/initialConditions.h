/*
 * initialConditions.h
 *
 *  Created on: Feb 27, 2017
 *      Author: stephendewitt
 */

#ifndef INCLUDE_INITIALCONDITIONS_H_
#define INCLUDE_INITIALCONDITIONS_H_

#include "userInputParameters_pf.h"
#include "multiPhysicsBVP.h"

template <int dim, int degree>
class InitialCondition : public dealii::Function<dim>
{
public:
  const unsigned int index;
  const userInputParameters_pf<dim> userInputs_pf;
  dealii::Vector<double> values;
  InitialCondition (const unsigned int _index, const userInputParameters_pf<dim> _userInputs_pf, MultiPhysicsBVP<dim,degree>* _multiphysics_pde) : dealii::Function<dim>(1), index(_index), userInputs_pf(_userInputs_pf),multiphysics_pde(_multiphysics_pde) {
    std::srand(dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)+1);
  }
  // IC for scalar values
  double value (const dealii::Point<dim> &p, const unsigned int component=0) const {
      double scalar_IC = 0.0;
      dealii::Vector<double> vector_IC(dim);

      multiphysics_pde->setInitialCondition(p, index, scalar_IC, vector_IC);
      return scalar_IC;
  };

private:
    MultiPhysicsBVP<dim,degree>* multiphysics_pde;
};

template <int dim, int degree>
class InitialConditionVector : public dealii::Function<dim>
{
public:
  const unsigned int index;
  const userInputParameters_pf<dim> userInputs_pf;
  dealii::Vector<double> values;
  InitialConditionVector (const unsigned int _index, const userInputParameters_pf<dim> _userInputs_pf, MultiPhysicsBVP<dim,degree>* _multiphysics_pde) : dealii::Function<dim>(dim), index(_index), userInputs_pf(_userInputs_pf),multiphysics_pde(_multiphysics_pde) {
    std::srand(dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)+1);
  }

  // IC for vector values
  void vector_value (const dealii::Point<dim> &p,dealii::Vector<double> &vector_IC) const {
      double scalar_IC = 0.0;
      vector_IC.reinit(dim);
      multiphysics_pde->setInitialCondition(p, index, scalar_IC, vector_IC);
  };

private:
    MultiPhysicsBVP<dim,degree>* multiphysics_pde;
};

#endif /* INCLUDE_INITIALCONDITIONS_H_ */
