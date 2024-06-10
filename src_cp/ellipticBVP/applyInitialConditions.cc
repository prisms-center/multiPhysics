//methods to apply initial conditions
#include "../../include_cp/ellipticBVP.h"

//methods to apply initial conditions
template <int dim>
void ellipticBVP<dim>::applyInitialConditions(){
  //pcout << "applying the default zero initial condition\n";
  //default method to apply zero initial conditions on all fields
  VectorTools::interpolate (dofHandler,
			    Functions::ZeroFunction<dim>(dim),
			    solution);
}
#include "../../include_cp/ellipticBVP_template_instantiations.h"
