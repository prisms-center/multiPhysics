/*
 * typeDefs.h
 *
 *  Created on: Feb 24, 2017
 *      Author: stephendewitt
 */

//#ifndef INCLUDE_TYPEDEFS_H_
//#define INCLUDE_TYPEDEFS_H_

//#include <deal.II/base/quadrature.h>
//#include <deal.II/base/timer.h>
//#include <deal.II/lac/vector.h>
//#include <deal.II/lac/affine_constraints.h>
//#include <deal.II/fe/fe_system.h>
//#include <deal.II/fe/fe_q.h>
//#include <deal.II/fe/fe_values.h>
//#include <deal.II/grid/tria.h>
//#include <deal.II/grid/tria_accessor.h>
//#include <deal.II/grid/tria_iterator.h>
//#include <deal.II/grid/grid_tools.h>
//#include <deal.II/dofs/dof_tools.h>
//#include <deal.II/dofs/dof_handler.h>
// #include <deal.II/numerics/vector_tools.h>
// #include <deal.II/lac/la_parallel_vector.h>
// #include <deal.II/matrix_free/matrix_free.h>
// #include <deal.II/matrix_free/fe_evaluation.h>
// #include <deal.II/base/config.h>
// #include <deal.II/base/exceptions.h>
// #include <deal.II/distributed/tria.h>
//#include <deal.II/distributed/solution_transfer.h>
//#include <deal.II/grid/manifold_lib.h>

//define PF data types
#ifndef scalarType_pf
typedef dealii::VectorizedArray<double> scalarType_pf;
#endif
#ifndef vectorType_pf
typedef dealii::LinearAlgebra::distributed::Vector<double> vectorType_pf;
#endif
//define FE system types
#ifndef typeScalar
typedef dealii::FEEvaluation<dim,degree,degree+1,1,double>  typeScalar;
#endif
#ifndef typeVector
typedef dealii::FEEvaluation<dim,degree,degree+1,dim,double>  typeVector;
#endif
//define data value types
#ifndef scalarvalueType_pf
typedef dealii::VectorizedArray<double> scalarvalueType_pf;
#endif
#ifndef vectorvalueType_pf
typedef dealii::Tensor<1, dim, dealii::VectorizedArray<double> > vectorvalueType_pf;
#endif
#if problemDIM==1
#ifndef scalargradType_pf
typedef dealii::VectorizedArray<double> scalargradType_pf;
#endif
#ifndef vectorgradType_pf
typedef dealii::VectorizedArray<double> vectorgradType_pf;
#endif
#ifndef vectorhessType_pf
typedef dealii::VectorizedArray<double> vectorhessType_pf;
#endif
#else
#ifndef scalargradType_pf
typedef dealii::Tensor<1, dim, dealii::VectorizedArray<double> > scalargradType_pf;
#endif
#ifndef scalarhessType_pf
typedef dealii::Tensor<2,dim,dealii::VectorizedArray<double> > scalarhessType_pf;
#endif
#ifndef vectorgradType_pf
typedef dealii::Tensor<2, dim, dealii::VectorizedArray<double> > vectorgradType_pf;
#endif
#ifndef vectorhessType_pf
typedef dealii::Tensor<3, dim, dealii::VectorizedArray<double> > vectorhessType_pf;
#endif
#endif

//#endif /* INCLUDE_TYPEDEFS_H_ */
