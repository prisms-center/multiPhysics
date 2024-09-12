//run method for multiPhysicsBVP class
#include "../../../include/crystalPlasticity.h"
#include <sys/stat.h>

//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

template <int dim>
void crystalPlasticity<dim>::run(){

  const int dir_err = mkdir(this->userInputs_cp.outputDirectory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  //initialization
  #if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 3)&&(DEAL_II_VERSION_MAJOR==9)))
  this->computing_timer_cp.enter_section("mesh and initialization");
  #else
  this->computing_timer_cp.enter_subsection("mesh and initialization");
  #endif

  //READING AND INITIALIZATION (PRISMS-Plasticity)
  //read mesh;
  this->mesh();
  //initialize FE objects and global data structures
  this->init_cp();
  this->initProjection();
  #if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 3)&&(DEAL_II_VERSION_MAJOR==9)))
  this->computing_timer_cp.exit_section("mesh and initialization");
  #else
  this->computing_timer_cp.leave_subsection("mesh and initialization");
  #endif
  //this->solve_cp();

  //READING AND INITIALIZATION (PRISMS-PF)
  pf_object.buildFields();
  pf_object.init_pf();
  pf_object.solve_pf();

}
#include "../../../include/crystalPlasticity_template_instantiations.h"
