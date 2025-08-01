//Methods related to the user model functionality
#include "../../include/multiPhysicsBVP.h"

#ifdef enableUserModel
template <int dim, int degree>
void MultiPhysicsBVP<dim,degree>::initQuadHistory(){
  QGauss<dim>  quadrature(quadOrder);
  const unsigned int   num_quad_points = quadrature.size();
  //initialize the quadHistory table
  quadHistory.reinit(TableIndices<3> (triangulation.n_locally_owned_active_cells(), num_quad_points, numQuadHistoryVariables));
}
#endif
#include "../../include/multiPhysicsBVP_template_instantiations.h"
