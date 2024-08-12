//methods to allow for pre/post iteration level updates
#include "../../include/multiPhysicsBVP.h"

//method called before each iteration
template <int dim>
void multiPhysicsBVP<dim>::updateBeforeIteration(){
  //default method does nothing
}

//method called after each iteration
template <int dim>
void multiPhysicsBVP<dim>::updateAfterIteration(){
  //default method does nothing
}

//method called after each iteration
template <int dim>
bool multiPhysicsBVP<dim>::testConvergenceAfterIteration(){
  //default method resets solution to previously converged solution if resetIncrement flagis true
  if (resetIncrement){
    solution=oldSolution;
    solutionWithGhosts=oldSolution;

    resetIncrement=false;
    char buffer[100];
    sprintf(buffer,
	    "current increment reset by model. Restarting increment with loadFactorSetByModel: %12.6e\n",
	    loadFactorSetByModel);
    pcout << buffer;
    return false;
  }
  return true;
}
#include "../../include/multiPhysicsBVP_template_instantiations.h"
