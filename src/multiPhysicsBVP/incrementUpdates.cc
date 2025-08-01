//methods to allow for pre/post increment level updates
#include "../../include/multiPhysicsBVP.h"

//method called before each increment
template <int dim, int degree>
void MultiPhysicsBVP<dim,degree>::updateBeforeIncrement(){
  //default method does nothing
  //Overwritten in crystal plasticity
  //pcout << "eBVP::updateBeforeIncrement\n";
  if (userInputs_cp.enableIndentationBCs){
  //    pcout << "enableIndentationBCs\n";
      updateIndentPos();
      //newton_rhs_uncondensed += n
  }
}

//method called after each increment
template <int dim, int degree>
void MultiPhysicsBVP<dim,degree>::updateAfterIncrementBase(){
    if (userInputs_cp.enableIndentationBCs){
        measureIndentationLoad();
        indenterLoad = Utilities::MPI::sum(indenterLoad, mpi_communicator);
        pcout << "         Indenter Load: "
              << indenterLoad
              << std::endl;
    }
  //default method does nothing
}

template <int dim, int degree>
void MultiPhysicsBVP<dim,degree>::updateAfterIncrement(){

    //default method does nothing
}


#include "../../include/multiPhysicsBVP_template_instantiations.h"
