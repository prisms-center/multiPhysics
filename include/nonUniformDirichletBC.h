#ifndef INCLUDE_NONUNIFORMDIRICHLETBCS_H_
#define INCLUDE_NONUNIFORMDIRICHLETBCS_H_

#include "userInputParameters_pf.h"
#include "multiPhysicsBVP.h"

template <int dim, int degree>
class NonUniformDirichletBC : public dealii::Function<dim>
{
public:

  dealii::Vector<double> values;

  NonUniformDirichletBC (const unsigned int _index, const unsigned int _direction, const double _time, MultiPhysicsBVP<dim,degree>* _multiphysics_pde) : dealii::Function<dim>(1), index(_index), direction(_direction), time(_time), multiphysics_pde(_multiphysics_pde) {
    std::srand(dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)+1);
  }
  // IC for scalar values
  double value (const dealii::Point<dim> &p, const unsigned int component=0) const {
      double scalar_BC = 0.0;
      dealii::Vector<double> vector_BC(dim);

      multiphysics_pde->setNonUniformDirichletBCs(p, index, direction, time, scalar_BC, vector_BC);

      return scalar_BC;
  };

private:
    const unsigned int index;
    const unsigned int direction;
    const double time;
    MultiPhysicsBVP<dim,degree>* multiphysics_pde;
};


template <int dim, int degree>
class NonUniformDirichletBCVector : public dealii::Function<dim>
{
public:

  dealii::Vector<double> values;

  NonUniformDirichletBCVector (const unsigned int _index, const unsigned int _direction, const double _time, MultiPhysicsBVP<dim,degree>* _multiphysics_pde) : dealii::Function<dim>(dim), index(_index), direction(_direction), time(_time), multiphysics_pde(_multiphysics_pde) {
    std::srand(dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)+1);
  }

  // IC for vector values
  void vector_value (const dealii::Point<dim> &p,dealii::Vector<double> &vector_BC) const {
      double scalar_BC = 0.0;

      vector_BC.reinit(dim);
      multiphysics_pde->setNonUniformDirichletBCs(p, index, direction, time, scalar_BC, vector_BC);
  };

private:
    const unsigned int index;
    const unsigned int direction;
    const double time;
    MultiPhysicsBVP<dim,degree>* multiphysics_pde;
};


#endif // INCLUDE_NONUNIFORMDIRICHLETBCS_H_
