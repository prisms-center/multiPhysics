//List of deal.II headers needed for the plasticity codes

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/timer.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/grid/grid_refinement.h>

#if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 1)&&(DEAL_II_VERSION_MAJOR==9)))
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/constraint_matrix.h>
#else
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/affine_constraints.h>
#endif

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

//List of deal.II headers needed for the prisms-pf codes
#include <deal.II/base/quadrature.h>
//#include <deal.II/base/timer.h>
//#include <deal.II/lac/vector.h>
#include <deal.II/lac/affine_constraints.h>
//#include <deal.II/fe/fe_system.h>
//#include <deal.II/fe/fe_q.h>
//#include <deal.II/fe/fe_values.h>
#if (DEAL_II_VERSION_MAJOR == 9 && DEAL_II_VERSION_MINOR > 3)
#include <deal.II/fe/mapping_fe.h>
#endif
//#include <deal.II/grid/tria.h>
//#include <deal.II/grid/tria_accessor.h>
//#include <deal.II/grid/tria_iterator.h>
//#include <deal.II/grid/grid_tools.h>
//#include <deal.II/dofs/dof_tools.h>
//#include <deal.II/dofs/dof_handler.h>
//#include <deal.II/numerics/vector_tools.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
//#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/grid/manifold_lib.h>