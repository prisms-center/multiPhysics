//base class for matrix Free implementation of PDE's
#ifndef MULTIPHYSICSBVP_H
#define MULTIPHYSICSBVP_H

// general headers
#include <fstream>
#include <sstream>
#include <iterator>

//dealii headers
#include "dealIIheaders.h"

using namespace dealii;

// PRISMS-PF headers
#include "fields.h"
#include "userInputParameters_pf.h"
#include "nucleus.h"
#include "variableValueContainer.h"
#include "variableContainer.h"
#include "SimplifiedGrainRepresentation.h"

// PRISMS-Plasticity headers
#include "userInputParameters_cp.h"
#include "crystalOrientationsIO.h"

//compiler directives to handle warnings
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic pop

//define data types
typedef PETScWrappers::MPI::Vector vectorType_cp;
typedef PETScWrappers::MPI::SparseMatrix matrixType_cp;

// define data types
#ifndef scalarType_pf
typedef dealii::VectorizedArray<double> scalarType_pf;
#endif
#ifndef vectorType_pf
typedef dealii::LinearAlgebra::distributed::Vector<double> vectorType_pf;
#endif

//macro for constants
#define constV(a) make_vectorized_array(a)

template <int dim, int degree>
class MultiPhysicsBVP:public Subscriptor
{
public:
  /**
  * Class contructor
  */
  MultiPhysicsBVP(userInputParameters_pf<dim>, userInputParameters_cp _userInputs_cp);
  ~MultiPhysicsBVP();

  //PRISMS-PF functions
  /**
  * Initializes the mesh, degrees of freedom, constraints and data structures using the user provided
  * inputs in the application parameters file.
  */
  
  virtual void init_pf();

  virtual void makeTriangulation(parallel::distributed::Triangulation<dim> &) const;

  /**
  * Initializes the data structures for enabling unit tests.
  *
  * This method initializes the MatrixFreePDE object with a fixed geometry, discretization and
  * other custom selected options specifically to help with unit tests, and should not be called
  * in any of the physical models.
  */
  void initForTests();

  /**
  * This method implements the time stepping algorithm and invokes the solveIncrement() method.
  */
  void solve_pf();
  /**
  * This method essentially converts the MatrixFreePDE object into a matrix object which can be
  * used with matrix free iterative solvers. Provides the A*x functionality for solving the system of
  * equations AX=b.
  */
  void vmult (vectorType_pf &dst, const vectorType_pf &src) const;
  /**
  * Vector of all the physical fields in the problem. Fields are identified by dimentionality (SCALAR/VECTOR),
  * the kind of PDE (ELLIPTIC/PARABOLIC) used to compute them and a character identifier  (e.g.: "c" for composition)
  * which is used to write the fields to the output files.
  */
  std::vector<Field<dim> >                  fields;

  void buildFields();

  // Parallel message stream
  ConditionalOStream  pcout;

  // Initial conditions function
  virtual void setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC) = 0;

  // Non-uniform boundary conditions function
  virtual void setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC) = 0;
  
  //PRISMS-Plasticity functions
  void run   ();
  crystalOrientationsIO<dim> orientations_Mesh;

protected:
  //PRISMS-PF functions
  /**
  * Initializes the mesh, degrees of freedom, constraints and data structures using the user provided
  * inputs in the application parameters file.
  */
  userInputParameters_pf<dim> userInputs_pf;

  unsigned int totalDOFs;

  //Total cpfe dofs
  unsigned int totalDOFs_cp;

  // Virtual methods to set the attributes of the primary field variables and the postprocessing field variables
  //virtual void setVariableAttriubutes() = 0;
  //virtual void setPostProcessingVariableAttriubutes(){};
  variableAttributeLoader var_attributes;

  // Elasticity matrix variables
  const static unsigned int CIJ_tensor_size = 2*dim-1+dim/3;

  // Method to reinitialize the mesh, degrees of freedom, constraints and data structures when the mesh is adapted
  void reinit  ();

  /**
  * Method to reassign grains when multiple grains are stored in a single order parameter.
  */
  void reassignGrains();

  std::vector<SimplifiedGrainRepresentation<dim>> simplified_grain_representations;

  /**
   * Method to solve each time increment of a time-dependent problem. For time-independent problems
   * this method is called only once. This method solves for all the fields in a staggered manner (one after another)
   * and also invokes the corresponding solvers: Explicit solver for Parabolic problems, Implicit (matrix-free) solver for Elliptic problems.
   */
  virtual void solveIncrement (bool skip_time_dependent);
  /* Method to write solution fields to vtu and pvtu (parallel) files.
  *
  * This method can be enabled/disabled by setting the flag writeOutput to true/false. Also,
  * the user can select how often the solution files are written by setting the flag
  * skipOutputSteps in the parameters file.
  */
  void outputResults();

  /*Parallel mesh object which holds information about the FE nodes, elements and parallel domain decomposition
  */
  parallel::distributed::Triangulation<dim> triangulation_pf;

  /*A vector of finite element objects used in a model. For problems with only one primal field,
  *the size of this vector is one,otherwise the size is the number of primal fields in the problem.
  */
  std::vector<FESystem<dim>*>          FESet;
  
  /*A vector of all the constraint sets in the problem. A constraint set is a map which holds the mapping between the degrees
  *of freedom and the corresponding degree of freedom constraints. Currently the type of constraints stored are either
  *Dirichlet boundary conditions or hanging node constraints for adaptive meshes.
  */
  std::vector<const AffineConstraints<double>*> constraintsDirichletSet, constraintsOtherSet;
  /*A vector of all the degree of freedom objects is the problem. A degree of freedom object handles the serial/parallel distribution
  *of the degrees of freedom for all the primal fields in the problem.*/
  std::vector<const DoFHandler<dim>*>  dofHandlersSet;

  /*A vector of the locally relevant degrees of freedom. Locally relevant degrees of freedom in a parallel implementation is a collection of the
  *degrees of freedom owned by the current processor and the surrounding ghost nodes which are required for the field computations in this processor.
  */
  std::vector<const IndexSet*>         locally_relevant_dofsSet;
  /*Copies of constraintSet elements, but stored as non-const to enable application of constraints.*/
  std::vector<AffineConstraints<double>*>       constraintsDirichletSet_nonconst, constraintsOtherSet_nonconst;
  /*Copies of dofHandlerSet elements, but stored as non-const.*/
  std::vector<DoFHandler<dim>*>        dofHandlersSet_nonconst;
  /*Copies of locally_relevant_dofsSet elements, but stored as non-const.*/
  std::vector<IndexSet*>               locally_relevant_dofsSet_nonconst;
  /*Vector all the solution vectors in the problem. In a multi-field problem, each primal field has a solution vector associated with it.*/
  std::vector<vectorType_pf*>             solutionSet;
  /*Vector all the residual (RHS) vectors in the problem. In a multi-field problem, each primal field has a residual vector associated with it.*/
  std::vector<vectorType_pf*>             residualSet;
  /*Vector of parallel solution transfer objects. This is used only when adaptive meshing is enabled.*/
  std::vector<parallel::distributed::SolutionTransfer<dim, vectorType_pf>*> soltransSet;

  // Objects for vectors
  DoFHandler<dim>* vector_dofHandler;
  FESystem<dim>* vector_fe;
  MatrixFree<dim,double> vector_matrixFreeObject;

  //matrix free objects
  /*Object of class MatrixFree<dim>. This is primarily responsible for all the base matrix free functionality of this MatrixFreePDE<dim> class.
  *Refer to deal.ii documentation of MatrixFree<dim> class for details.
  */
  MatrixFree<dim,double>               matrixFreeObject;
  /*Vector to store the inverse of the mass matrix diagonal. Due to the choice of spectral elements with Guass-Lobatto quadrature, the mass matrix is diagonal.*/
  vectorType_pf                           invM;
  /*Vector to store the solution increment. This is a temporary vector used during implicit solves of the Elliptic fields.*/
  vectorType_pf                           dU_vector, dU_scalar;

  //matrix free methods
  /*Current field index*/
  unsigned int currentFieldIndex;
  /*Method to compute the inverse of the mass matrix*/
  void computeInvM();

  /*AMR methods*/
  /**
  * Method that actually changes the triangulation based on refine/coarsen flags set previously.
  */
  void refineGrid();

  /**
  * Method to control the overall flow of adaptive mesh refinement.
  */
  void adaptiveRefine(unsigned int _currentIncrement_pf);

  /**
  * Virtual method to define the the criterion for refining or coarsening the mesh. This method sets refine/coarsen flags that are read by refineGrid.
  */
  virtual void adaptiveRefineCriterion();

  /*Method to compute the right hand side (RHS) residual vectors*/
  void computeExplicitRHS();
  void computeNonexplicitRHS();

  //virtual methods to be implemented in the derived class
  /*Method to calculate LHS(implicit solve)*/
  void getLHS(const MatrixFree<dim,double> &data,
          vectorType_pf &dst,
          const vectorType_pf &src,
          const std::pair<unsigned int,unsigned int> &cell_range) const;


  bool generatingInitialGuess;
  void getLaplaceLHS(const MatrixFree<dim,double> &data,
              vectorType_pf &dst,
              const vectorType_pf &src,
              const std::pair<unsigned int,unsigned int> &cell_range) const;


  void setNonlinearEqInitialGuess();
  void computeLaplaceRHS(unsigned int fieldIndex);
  void getLaplaceRHS (const MatrixFree<dim,double> &data,
          vectorType_pf &dst,
          const vectorType_pf &src,
          const std::pair<unsigned int,unsigned int> &cell_range) const;


  /*Method to calculate RHS (implicit/explicit). This is an abstract method, so every model which inherits MatrixFreePDE<dim> has to implement this method.*/
  void getExplicitRHS (const MatrixFree<dim,double> &data,
          std::vector<vectorType_pf*> &dst,
          const std::vector<vectorType_pf*> &src,
          const std::pair<unsigned int,unsigned int> &cell_range) const;

  void getNonexplicitRHS (const MatrixFree<dim,double> &data,
            std::vector<vectorType_pf*> &dst,
            const std::vector<vectorType_pf*> &src,
            const std::pair<unsigned int,unsigned int> &cell_range) const;

  virtual void explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
                                                        dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const=0;

  virtual void nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
                                                        dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const=0;

  virtual void equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
                                dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const=0;

  virtual void postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
                                                              variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
                                                              const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const=0;
  void computePostProcessedFields(std::vector<vectorType_pf*> &postProcessedSet);

  void getPostProcessedFields(const dealii::MatrixFree<dim,double> &data,
                                                                                      std::vector<vectorType_pf*> &dst,
                                                                                      const std::vector<vectorType_pf*> &src,
                                                                                      const std::pair<unsigned int,unsigned int> &cell_range);

  //methods to apply dirichlet BC's
  /*Map of degrees of freedom to the corresponding Dirichlet boundary conditions, if any.*/
  std::vector<std::map<dealii::types::global_dof_index, double>*> valuesDirichletSet;
  /*Virtual method to mark the boundaries for applying Dirichlet boundary conditions.  This is usually expected to be provided by the user.*/
  void markBoundaries_pf(parallel::distributed::Triangulation<dim> &) const;
  /** Method for applying Dirichlet boundary conditions.*/
  void applyDirichletBCs_pf();

  /** Method for applying Neumann boundary conditions.*/
  void applyNeumannBCs();

  // Methods to apply periodic BCs
  void setPeriodicity_pf();
  void setPeriodicityConstraints_pf(AffineConstraints<double>*, const DoFHandler<dim>*) const;
  void getComponentsWithRigidBodyModes(std::vector<int> &) const;
  void setRigidBodyModeConstraints(const std::vector<int>, AffineConstraints<double>*, const DoFHandler<dim>*) const;

  //methods to apply initial conditions
  /*Virtual method to apply initial conditions.  This is usually expected to be provided by the user in IBVP (Initial Boundary Value Problems).*/

  void applyInitialConditions_pf();

  // --------------------------------------------------------------------------
  // Methods for saving and loading checkpoints
  // --------------------------------------------------------------------------

  void save_checkpoint();

  void load_checkpoint_triangulation();
  void load_checkpoint_fields();
  void load_checkpoint_time_info();

  void move_file(const std::string&, const std::string&);

  void verify_checkpoint_file_exists(const std::string filename);

  // --------------------------------------------------------------------------
  // Nucleation methods and variables
  // --------------------------------------------------------------------------
  // Vector of all the nuclei seeded in the problem
  std::vector<nucleus<dim> > nuclei;

  // Method to get a list of new nuclei to be seeded
  void updateNucleiList();
  std::vector<nucleus<dim> > getNewNuclei();
  void getLocalNucleiList(std::vector<nucleus<dim> > & newnuclei) const;
  void safetyCheckNewNuclei(std::vector<nucleus<dim> > newnuclei, std::vector<unsigned int> & conflict_ids);
  void refineMeshNearNuclei(std::vector<nucleus<dim> > newnuclei);
  double weightedDistanceFromNucleusCenter(const dealii::Point<dim,double> center, const std::vector<double> semiaxes, const dealii::Point<dim,double> q_point_loc, const unsigned int var_index) const;
  dealii::VectorizedArray<double> weightedDistanceFromNucleusCenter(const dealii::Point<dim,double> center, const std::vector<double> semiaxes, const dealii::Point<dim,dealii::VectorizedArray<double> > q_point_loc, const unsigned int var_index) const;


  // Method to obtain the nucleation probability for an element, nontrival case must be implemented in the subclass
  virtual double getNucleationProbability(variableValueContainer, double, dealii::Point<dim>, unsigned int) const {return 0.0;};
  //utility functions
  /*Returns index of given field name if exists, else throw error.*/
  unsigned int getFieldIndex(std::string _name);

  std::vector<double> freeEnergyValues;
  void outputFreeEnergy(const std::vector<double>& freeEnergyValues) const;

  /*Method to compute the integral of a field.*/
  void computeIntegral(double& integratedField, int index, std::vector<vectorType_pf*> postProcessedSet);

  //variables for time dependent problems
  /*Flag used to see if invM, time stepping in run(), etc are necessary*/
  bool isTimeDependentBVP;
  /*Flag used to mark problems with Elliptic fields.*/
  bool isEllipticBVP;

  bool hasExplicitEquation;
  bool hasNonExplicitEquation;
  bool has_Dirichlet_BCs;
  //
  unsigned int parabolicFieldIndex, ellipticFieldIndex;
  double currentTime_pf;
  unsigned int currentIncrement_pf, currentOutput, currentCheckpoint, current_grain_reassignment;

  /*Timer and logging object*/
  mutable TimerOutput computing_timer_pf;

  std::vector<double> integrated_postprocessed_fields;

  bool first_integrated_var_output_complete;

  // Methods and variables for integration
  double integrated_var;
  unsigned int integral_index;
  std::mutex assembler_lock;

  void computeIntegralMF(double& integratedField, int index, const std::vector<vectorType_pf*> postProcessedSet);

  void getIntegralMF (const MatrixFree<dim,double> &data,
          std::vector<vectorType_pf*> &dst,
          const std::vector<vectorType_pf*> &src,
          const std::pair<unsigned int,unsigned int> &cell_range);

  //PRISMS-Plasticity functions

  //parallel objects
  MPI_Comm   mpi_communicator;
  IndexSet   locally_owned_dofs;
  IndexSet   locally_owned_dofs_Scalar;
  IndexSet   locally_relevant_dofs,locally_relevant_dofs_Mod, locally_relevant_ghost_dofs;
  IndexSet   locally_relevant_dofs_Scalar;
  IndexSet   dof_FN,dof_FP;
  IndexSet   dof_FXN,dof_FXP;
  IndexSet   dof_FYN,dof_FYP;
  IndexSet   dof_FZN,dof_FZP;
  IndexSet   dof_Boundary_Layer2,vertices_DOFs,global_dof_Boundary_Layer2,Edges_DOFs_1,Edges_DOFs_2,Edges_DOFs_3;
  //INDENTATION
  IndexSet      active_set, debug_set, frozen_set;
  unsigned int active_set_size, old_active_set_size, freeze_out_iterations;
  bool  active_set_size_changed;

  unsigned int globalDOF_V000_1,globalDOF_V000_2,globalDOF_V000_3,globalDOF_V001_1,globalDOF_V001_2,globalDOF_V001_3,globalDOF_V010_1,globalDOF_V010_2,globalDOF_V010_3,globalDOF_V100_1,globalDOF_V100_2,globalDOF_V100_3;
  unsigned int globalDOF_V011_1,globalDOF_V011_2,globalDOF_V011_3,globalDOF_V110_1,globalDOF_V110_2,globalDOF_V110_3,globalDOF_V101_1,globalDOF_V101_2,globalDOF_V101_3,globalDOF_V111_1,globalDOF_V111_2,globalDOF_V111_3;

  unsigned int local_globalDOF_V000_1,local_globalDOF_V000_2,local_globalDOF_V000_3,local_globalDOF_V001_1,local_globalDOF_V001_2,local_globalDOF_V001_3,local_globalDOF_V010_1,local_globalDOF_V010_2,local_globalDOF_V010_3,local_globalDOF_V100_1,local_globalDOF_V100_2,local_globalDOF_V100_3;
  unsigned int local_globalDOF_V011_1,local_globalDOF_V011_2,local_globalDOF_V011_3,local_globalDOF_V110_1,local_globalDOF_V110_2,local_globalDOF_V110_3,local_globalDOF_V101_1,local_globalDOF_V101_2,local_globalDOF_V101_3,local_globalDOF_V111_1,local_globalDOF_V111_2,local_globalDOF_V111_3;

  unsigned int global_size_dof_Boundary_Layer2,global_size_dof_Edge;

  std::vector<unsigned int> global_vector_dof_Boundary_Layer2,dofNodalDisplacement;
  std::vector<double> deluNodalDisplacement,initPosIndenter,finalPosIndenter, prevPosIndenter, currentPosIndenter;
  std::vector<Point<dim>> KeyPosIndenter;
  double indenterSize, indenterTolerance, indenterLoad;// depthRefinementMultiple;
  unsigned int indenterShape, indenterFace, loadFace, indentDof, indentationKeyFrames;//, refinementFactor;
  bool roughIndenter;
  std::vector<std::vector<unsigned int> >  vertices_Constraints_Matrix,edges_Constraints_Matrix,faces_Constraints_Matrix, global_Edges_DOFs_Vector_Array;
  std::vector<std::vector<double>>  vertices_Constraints_Coef,edges_Constraints_Coef,faces_Constraints_Coef, global_Edges_DOFs_Coord_Vector_Array,nodalDisplacement;
  std::vector<std::vector<int>> periodicBCsInput; // 	Periodic BCs Input
  std::vector<std::vector<double>> periodicBCsInput2,periodicBCsInput2_Orig; // 	Periodic BCs Input
  std::vector<unsigned int> vertices_Constraint_Known, vertices_DOFs_vector;
  std::vector<std::vector<unsigned int>> edges_DOFs;
  std::vector<IndexSet> faces_dof_Index_vector;

  unsigned int numberVerticesConstraint,numberEdgesConstraint,numberFacesConstraint,totalNumVerticesDOFs,totalNumEdgesDOFs,totalNumEachEdgesNodes,totalNumFacesDOFs;

  Point<dim> vertexNode,node, cellCenter, nodeU, nodedU, nodeU2;
  std::vector<double> displace_local;

  //User input parameters object
  userInputParameters_cp userInputs_cp;

  //FE data structres
  parallel::distributed::Triangulation<dim> triangulation_cp,triangulation2_cp;
  FESystem<dim>      FE;
  FESystem<dim>      FE_Scalar;
  DoFHandler<dim>    dofHandler;
  DoFHandler<dim>    dofHandler_Scalar;
  std::vector<unsigned int> cellOrientationMap_Mesh;

  //methods
  virtual void mesh();
  void init_cp();
  void assemble();
  #if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 1)&&(DEAL_II_VERSION_MAJOR==9)))
  ConstraintMatrix   constraints, constraints_PBCs_Inc0, constraints_PBCs_IncNot0, constraints_PBCs_Inc0Neg;
  ConstraintMatrix   constraintsMassMatrix;
  void solveLinearSystem(ConstraintMatrix& constraintmatrix, matrixType_cp& A, vectorType_cp& b, vectorType_cp& x, vectorType_cp& xGhosts, vectorType_cp& dxGhosts);
  void solveLinearSystem2(ConstraintMatrix& constraintmatrix, matrixType_cp& A, vectorType_cp& b, vectorType_cp& x, vectorType_cp& xGhosts, vectorType_cp& dxGhosts);
  ///Periodic BCs Implementation
  void setFaceConstraints(ConstraintMatrix& constraintmatrix);
  void setEdgeConstraints(ConstraintMatrix& constraintmatrix);
  void setNodeConstraints(ConstraintMatrix& constraintmatrix);
  #else
  AffineConstraints<double>   constraints, constraints_PBCs_Inc0, constraints_PBCs_IncNot0, constraints_PBCs_Inc0Neg;
  AffineConstraints<double>   constraintsMassMatrix;
  //INDENTATION
  AffineConstraints<double>   indentation_constraints, hanging_constraints;
  void solveLinearSystem(AffineConstraints<double>& constraintmatrix, matrixType_cp& A, vectorType_cp& b, vectorType_cp& x, vectorType_cp& xGhosts, vectorType_cp& dxGhosts);
  void solveLinearSystem2(AffineConstraints<double>& constraintmatrix, matrixType_cp& A, vectorType_cp& b, vectorType_cp& x, vectorType_cp& xGhosts, vectorType_cp& dxGhosts);
  ///Periodic BCs Implementation
  void setNodeConstraints(AffineConstraints<double>& constraintmatrix);
  void setEdgeConstraints(AffineConstraints<double>& constraintmatrix);
  void setFaceConstraints(AffineConstraints<double>& constraintmatrix);
  #endif

  bool solveNonLinearSystem();
  void solve_cp();
  void output();
  void initProjection();
  void projection();
  void markBoundaries_cp();

  //virtual methods to be implemented in derived class
  //method to calculate elemental Jacobian and Residual,
  //which should be implemented in the derived material model class

  virtual void getElementalValues(FEValues<dim>& fe_values, unsigned int dofs_per_cell, unsigned int num_quad_points, FullMatrix<double>& elementalJacobian, Vector<double>&elementalResidual) = 0;

  //methods to allow for pre/post iteration updates
  virtual void updateBeforeIteration();
  virtual void updateAfterIteration();
  virtual bool testConvergenceAfterIteration();
  //methods to allow for pre/post increment updates
  virtual void updateBeforeIncrement();
  virtual void updateAfterIncrement();
  void updateAfterIncrementBase();
  //methods to apply dirichlet BC's and initial conditions
  void applyDirichletBCs_cp();
  void applyInitialConditions_cp();
  void setBoundaryValues(const Point<dim>& node, const unsigned int dof, bool& flag, double& value);
  ///////These functions are for Periodic BCs Implementation
  void setPeriodicity_cp();
  void setPeriodicityConstraintsInit();
  void setPeriodicityConstraints_cp();
  void setPeriodicityConstraintsInc0();
  void setPeriodicityConstraintsIncNot0();
  void setPeriodicityConstraintsInc0Neg();
  ///////These functions are for Indentation BCs Implementation
  void updateIndentPos();
  void setIndentation(const Point<dim>& node, const unsigned int dof, bool& flag, double& value);
  void setIndentation2(const Point<dim>& node, const unsigned int dof, bool& flag, double& value,
                      double& criterion);
  void setIndentationConstraints();
  void displaceNode(const Point<dim> & p, const Point<dim> & u);
  void meshRefineIndentation();
  bool flagActiveSet(const Point<dim> & p);
  bool flagActiveSetLambda(const Point<dim> & p, double& criterion);
  bool flagActiveSetLambda2(const Point<dim> & p, double& criterion);
  //void updateActiveSet();
  void measureIndentationLoad();
  void setActiveSet();
  void setActiveSet2();
  void setFrozenSet();
  void assemble_mass_matrix_diagonal();
  double Obstacle(const Point<dim> & p, const unsigned int & component, const std::vector<double> & ind);
  ///////These functions are for DIC BCs evaluation
  void bcFunction1(double _yval, double &value_x, double &value_y, double _currentIncr);
  void bcFunction2(double _yval, double &value_x, double &value_y, double _currentIncr);
  void bcFunction3(double _yval, double &value_x, double &value_y, double _currentIncr);
  void bcFunction4(double _yval, double &value_x, double &value_y, double _currentIncr);
  std::map<types::global_dof_index, Point<dim> > supportPoints;
  //parallel data structures
  vectorType_cp solution, oldSolution, residual;
  vectorType_cp solutionWithGhosts, solutionIncWithGhosts;
  //INDENTATION
  vectorType_cp  newton_rhs_uncondensed, newton_rhs_uncondensed_inc, diag_mass_matrix_vector;
  vectorType_cp  lambda;
  matrixType_cp jacobian;
  // Boundary condition variables
  std::vector<std::vector<bool>> faceDOFConstrained;
  std::vector<std::vector<double>> deluConstraint;
  FullMatrix<double> tabularDisplacements,tabularInputNeumannBCs;
  unsigned int timeCounter,timeCounter_Neumann;
  double currentTime_cp;
  /////DIC bc names
  FullMatrix<double> bc_new1,bc_new2,bc_new3,bc_new4;
  FullMatrix<double> Fprev=IdentityMatrix(dim);
  FullMatrix<double> F,deltaF;
  FullMatrix<double> targetVelGrad;
  //misc variables
  double delT,totalT,cycleTime;
  unsigned int currentIteration, currentIncrement_cp;
  unsigned int totalIncrements,periodicTotalIncrements;
  bool resetIncrement;
  double loadFactorSetByModel;
  double totalLoadFactor;
  //INDENTATION
  double currentIndentDisp;
  //parallel message stream (already defined in the PRISMS-PF section)
  //ConditionalOStream  pcout;
  //compute-time logger
  TimerOutput computing_timer_cp;
  //output variables
  //solution name array
  std::vector<std::string> nodal_solution_names;
  std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation;
  //post processing
  unsigned int numPostProcessedFields;
  unsigned int numPostProcessedFieldsAtCellCenters;
  //postprocessed scalar variable name array (only scalar variables supported currently, will be extended later to vectors and tensors, if required.)
  std::vector<std::string> postprocessed_solution_names, postprocessedFieldsAtCellCenters_solution_names;
  //postprocessing data structures
  std::vector<vectorType_cp*> postFields, postFieldsWithGhosts, postResidual;
  matrixType_cp massMatrix;
  Table<4,double> postprocessValues;
  Table<2,double> postprocessValuesAtCellCenters;
  //user model related variables and methods
  #ifdef enableUserModel
  unsigned int numQuadHistoryVariables;
  Table<3, double> quadHistory;
  virtual void initQuadHistory();
  #endif
};

#endif