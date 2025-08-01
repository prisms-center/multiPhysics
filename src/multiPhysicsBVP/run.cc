//run method for multiPhysicsBVP class
#include "../../include/multiPhysicsBVP.h"
#include <sys/stat.h>

//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

template <int dim, int degree>
void MultiPhysicsBVP<dim,degree>::run(){

  const int dir_err = mkdir(userInputs_cp.outputDirectory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  //initialization
  #if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 3)&&(DEAL_II_VERSION_MAJOR==9)))
  computing_timer_cp.enter_section("mesh and initialization");
  #else
  computing_timer_cp.enter_subsection("mesh and initialization");
  #endif

  //READING AND INITIALIZATION (PRISMS-Plasticity)
  //read mesh;
  mesh();
  //initialize FE objects and global data structures
  init_cp();
  initProjection();
  #if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 3)&&(DEAL_II_VERSION_MAJOR==9)))
  computing_timer_cp.exit_section("mesh and initialization");
  #else
  computing_timer_cp.leave_subsection("mesh and initialization");
  #endif

  //READING AND INITIALIZATION (PRISMS-PF)
  //pf_object.buildFields();
  //init_pf();

  //solve_cp();
  //solve_pf();
}
#include "../../include/multiPhysicsBVP_template_instantiations.h"
