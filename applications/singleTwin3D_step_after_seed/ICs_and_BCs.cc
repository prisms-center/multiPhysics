#include "customPDE.h"

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

  double center[3] = {0.5,0.5,0.5};
  dealii::Tensor<1, dim> nc;
  dealii::Tensor<1, dim> n;
  dealii::Tensor<1, dim> n_int;
  dealii::Tensor<1, dim> td_int;
  dealii::Tensor<1, dim> tn_int;
  nc.clear(); n.clear(); n_int.clear(); td_int.clear(); tn_int.clear();
  double dist, edist;
  double b0=a0*std::sqrt((1.0-ecc*ecc));
  double nXc, nYc, nZc, nX, nY, nZ, nZ_reg;
  double pi = 3.14159265358979323846;
  scalar_IC = 0;

  if (index==0){
    //Calculating distance from center of the system
    dist = 0.0;
    for (unsigned int dir = 0; dir < dim; dir++){
      dist += (p[dir]-center[dir]*userInputs_pf.domain_size[dir])*(p[dir]-center[dir]*userInputs_pf.domain_size[dir]);
    }
    dist = std::sqrt(dist);
    
    //Elliptical seed
    //Calculating distance from center to perimeter of the ellipse
    //Components of the normal vector from the center to point p
    nc[0] = (p[0]-center[0]*userInputs_pf.domain_size[0])/(dist + 1.0e-7);
    nc[1] = (p[1]-center[1]*userInputs_pf.domain_size[1])/(dist + 1.0e-7);
    nc[2] = (p[2]-center[2]*userInputs_pf.domain_size[2])/(dist + 1.0e-7);

    //Rotating td and tn according to Euler angles
    // Euler angles in radians (ZXZ convention)
    double phi1 = pi*euler_angs[0]/180.0; // e.g., 0.785398 for 45 degrees
    double Phi  = pi*euler_angs[1]/180.0; // e.g., 1.0472   for 60 degrees
    double phi2 = pi*euler_angs[2]/180.0; // e.g., 0.523599 for 30 degrees

    // Step 1: Build the rotation matrices Rz(phi1), Rx(Phi), Rz(phi2)
    dealii::Tensor<2, dim> Rz1, Rx, Rz2;
    Rz1.clear(); Rx.clear(); Rz2.clear();

    // Initialize identity diagonals first
    for (unsigned int i = 0; i < dim; ++i) {
    Rz1[i][i] = 1.0;
    Rx[i][i] = 1.0;
    Rz2[i][i] = 1.0;
    }

    // Rz(phi1)
    Rz1[0][0] = std::cos(phi1); Rz1[0][1] = -std::sin(phi1);
    Rz1[1][0] = std::sin(phi1); Rz1[1][1] =  std::cos(phi1);

    // Rx(Phi)
    Rx[1][1] = std::cos(Phi); Rx[1][2] = -std::sin(Phi);
    Rx[2][1] = std::sin(Phi); Rx[2][2] =  std::cos(Phi);

    // Rz(phi2)
    Rz2[0][0] = std::cos(phi2); Rz2[0][1] = -std::sin(phi2);
    Rz2[1][0] = std::sin(phi2); Rz2[1][1] =  std::cos(phi2);

    // Step 2: Compute Q = Rz1 * Rx * Rz2
    dealii::Tensor<2, dim> Q_temp, Q_1;
    Q_temp.clear(); Q_1.clear();

    // Q_temp = Rz1 * Rx
    for (unsigned int i = 0; i < dim; ++i)
        for (unsigned int j = 0; j < dim; ++j)
            for (unsigned int k = 0; k < dim; ++k)
                Q_temp[i][j] += Rz1[i][k] * Rx[k][j];

    // Q1 = Q_temp * Rz2
    for (unsigned int i = 0; i < dim; ++i)
        for (unsigned int j = 0; j < dim; ++j)
            for (unsigned int k = 0; k < dim; ++k)
                Q_1[i][j] += Q_temp[i][k] * Rz2[k][j];
    

    //Rotated normal vector, tn_int
    for (unsigned int i = 0; i < dim; ++i)
        for (unsigned int j = 0; j < dim; ++j)
            tn_int[i] += Q_1[i][j] * tn[j];

    //Rotated direction vector, td_int
    for (unsigned int i = 0; i < dim; ++i)
        for (unsigned int j = 0; j < dim; ++j)
            td_int[i] += Q_1[i][j] * td[j];


    //Rotated twin direction and twin normal vectors
    dealii::Tensor<1, dim> e_X = td_int;
    dealii::Tensor<1, dim> e_Y = tn_int;
    dealii::Tensor<1, dim> e_Z;

    // Compute e_Z = e_X cross e_Y
    e_Z[0] = e_X[1]*e_Y[2] - e_X[2]*e_Y[1];
    e_Z[1] = e_X[2]*e_Y[0] - e_X[0]*e_Y[2];
    e_Z[2] = e_X[0]*e_Y[1] - e_X[1]*e_Y[0];

    // Normalize e_Z
    double norm_eZ = std::sqrt(e_Z*e_Z);
    for (unsigned int i = 0; i < dim; ++i)
        e_Z[i] /= norm_eZ;

    // Construct rotation matrix Q as 3 column vectors
    dealii::Tensor<2, dim> Q;
    for (unsigned int i = 0; i < dim; ++i) {
        Q[i][0] = e_X[i];
        Q[i][1] = e_Y[i];
        Q[i][2] = e_Z[i];
    }
    //Intermediate vector, n_int
    for (unsigned int i = 0; i < 3; ++i)
        for (unsigned int j = 0; j < 3; ++j)
            n[i] += Q[j][i] * nc[j];

    //Rotated unit vector with respect to the twin plane    
    nX = n[0];
    nY = n[1];
    nZ = n[2];

    //Distance from center of the ellipse to a point in the ellipse that intersects the line defined by (nX,nY,nZ)
    //Regularized nZ such that nX^2 + nY^2 + nZ^2 = 1
    nZ_reg = std::sqrt(1.0-nX*nX-nY*nY);
    edist = 1.0/std::sqrt((nX/a0)*(nX/a0) + (nY/b0)*(nY/b0) + (nZ_reg/a0)*(nZ_reg/a0));
                     
    scalar_IC =  0.5*(1.0-std::tanh((dist-edist)/(1.0*del0*std::sqrt(edist/a0))));

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
