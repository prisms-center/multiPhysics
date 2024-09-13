#include "../../include/customPDE.h"

// ===========================================================================
// FUNCTION FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim, int degree>
void customPDE<dim,degree>::setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC){
    // ---------------------------------------------------------------------
    // ENTER THE INITIAL CONDITIONS HERE 
    // ---------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the initial condition for each variable
    // according to its variable index

    // The initial condition is a set of overlapping circles/spheres defined
    // by a hyperbolic tangent function. The center of each circle/sphere is
    // given by "center" and its radius is given by "radius".

  double center[2] = {0.5,0.5};
  double dist, edist;
  double b0=a0*std::sqrt((1.0-ecc*ecc));
  double cth, sth, nX, nY;
  double pi=3.141592;
  scalar_IC = 0;
  
  //Elliptical seed of semimajor axis
 
  if (index==0){
    //Calculating distance from center of the system
    dist = 0.0;
    for (unsigned int dir = 0; dir < 2; dir++){
      dist += (p[dir]-center[dir]*userInputs_pf.domain_size[dir])*(p[dir]-center[dir]*userInputs_pf.domain_size[dir]);
    }
    dist = std::sqrt(dist);
    
    //Calculating distance from center to perimeter of the ellipse
    cth=(p[0]-center[0]*userInputs_pf.domain_size[0])/(dist + 1.0e-7);
    sth=(p[1]-center[1]*userInputs_pf.domain_size[1])/(dist + 1.0e-7);
    
    //Rotated unit vector with respect to the twin plane
    
    nX=cth*std::cos(0.5*pi-th)-sth*std::sin(0.5*pi-th);
    nY=cth*std::sin(0.5*pi-th)+sth*std::cos(0.5*pi-th);

    //Distance to the center of the ellipse
    edist = a0*b0/std::sqrt(b0*b0*(1.0-nX*nX) + a0*a0*nX*nX);
                     
    scalar_IC =  0.5*(1.0-std::tanh((dist-edist)/(1.0*del0)));

    if (scalar_IC > 1.0) scalar_IC = 1.0;
                     
  } else {
    scalar_IC = 0.0;
  }
  // ---------------------------------------------------------------------
}

// ===========================================================================
// FUNCTION FOR NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS
// ===========================================================================

template <int dim, int degree>
void customPDE<dim,degree>::setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC)
{
    // --------------------------------------------------------------------------
    // ENTER THE NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS HERE
    // --------------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the boundary condition for each variable
    // according to its variable index. This function can be left blank if there
    // are no non-uniform Dirichlet boundary conditions. For BCs that change in
    // time, you can access the current time through the variable "time". The
    // boundary index can be accessed via the variable "direction", which starts
    // at zero and uses the same order as the BC specification in parameters.in
    // (i.e. left = 0, right = 1, bottom = 2, top = 3, front = 4, back = 5).


    // -------------------------------------------------------------------------

}

#include "../../include/customPDE_template_instantiations.h"