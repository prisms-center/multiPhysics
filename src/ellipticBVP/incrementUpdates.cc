//methods to allow for pre/post increment level updates
#include "../../include/multiPhysicsBVP.h"

//method called before each increment
template <int dim>
void multiPhysicsBVP<dim>::updateBeforeIncrement(){
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
template <int dim>
void multiPhysicsBVP<dim>::updateAfterIncrementBase(){
    if (userInputs_cp.enableIndentationBCs){
        measureIndentationLoad();
        indenterLoad = Utilities::MPI::sum(indenterLoad, mpi_communicator);
        pcout << "         Indenter Load: "
              << indenterLoad
              << std::endl;
    }
  //default method does nothing
}

template <int dim>
void multiPhysicsBVP<dim>::updateAfterIncrement(){

    //default method does nothing
}


#include "../../include/multiPhysicsBVP_template_instatiations.h"
