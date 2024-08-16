//methods to apply initial conditions
#include "../../include/multiPhysicsBVP.h"

//methods to apply initial conditions
template <int dim, int degree>
void MultiPhysicsBVP<dim,degree>::applyInitialConditions_cp(){
  //pcout << "applying the default zero initial condition\n";
  //default method to apply zero initial conditions on all fields
  VectorTools::interpolate (dofHandler,
			    Functions::ZeroFunction<dim>(dim),
			    solution);
}
#include "../../include/multiPhysicsBVP_template_instantiations.h"
