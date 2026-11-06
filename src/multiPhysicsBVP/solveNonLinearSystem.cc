//solve non linear system of equations for multiPhysicsBVP class
#include "../../include/multiPhysicsBVP.h"

// Newton iteration with
//  (1) a convergence exit (Absolute/Relative nonLinear solver tolerance, if set > 0), and
//  (2) a backtracking (damped) update: if the residual grows after a full Newton step,
//      the step is halved and the residual re-evaluated, up to maxBacktracks times.
//      Backtracking is only applied from the second correction on, so the boundary-
//      condition increment applied in the first iteration is never scaled.
template <int dim, int degree>
bool MultiPhysicsBVP<dim,degree>::solveNonLinearSystem(){
  //residuals
  double relNorm=1.0, initialNorm=1.0e-16, currentNorm=0.0;
  //INDENTATION
  unsigned int extra_non_linear_iterations=0;
  //non linear iterations
  char buffer[250];
  currentIteration=0;

  // ---- damping state (added) ----
  const unsigned int maxBacktracks = 6;      // step can shrink to 1/64
  unsigned int       nBacktracks   = 0;
  double             lambda        = 1.0;
  double             acceptedNorm  = std::numeric_limits<double>::max();
  bool               haveStep      = false;  // a Newton step that may be backtracked
  vectorType_cp      solutionBeforeStep(solution);
  vectorType_cp      stepIncrement(solution);
  const bool useAbsTol = userInputs_cp.absNonLinearTolerance > 0.0;
  const bool useRelTol = userInputs_cp.relNonLinearTolerance > 0.0;
  bool converged = false;

  while (currentIteration < userInputs_cp.maxNonLinearIterations + extra_non_linear_iterations){
    //call updateBeforeIteration, if any
    updateBeforeIteration();

    //Calling assemble
    #if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 3)&&(DEAL_II_VERSION_MAJOR==9)))
    computing_timer_cp.enter_section("assembly");
    assemble();
    computing_timer_cp.exit_section("assembly");
    #else
    computing_timer_cp.enter_subsection("assembly");
    assemble();
    computing_timer_cp.leave_subsection("assembly");
    #endif

    if (!resetIncrement){
      //Calculate residual norms
      currentNorm=residual.l2_norm();

      // ---- backtracking (added): residual grew (or is not finite) after the last step ----
      if (haveStep && currentIteration >= 2 && nBacktracks < maxBacktracks &&
          (!std::isfinite(currentNorm) || currentNorm > acceptedNorm))
        {
          lambda *= 0.5;
          nBacktracks++;
          sprintf(buffer,
            "  backtrack %u: residual %8.2e > %8.2e, retrying with step %6.4f\n",
            nBacktracks, currentNorm, acceptedNorm, lambda);
          pcout << buffer;
          solution = solutionBeforeStep;
          solution.add(lambda, stepIncrement);
          solutionWithGhosts = solution;
          continue; // re-assemble at the shortened step (does not count as an iteration)
        }
      if (nBacktracks == maxBacktracks && haveStep && currentNorm > acceptedNorm)
        pcout << "  backtracking exhausted, accepting step anyway\n";
      nBacktracks  = 0;
      lambda       = 1.0;
      acceptedNorm = currentNorm;

      initialNorm=std::max(initialNorm, currentNorm);
      relNorm=currentNorm/initialNorm;
      //print iteration information
      sprintf(buffer,
        "nonlinear iteration %3u [current residual: %8.2e, initial residual: %8.2e, relative residual: %8.2e]\n",
        currentIteration,
        currentNorm,
        initialNorm,
        relNorm);
      pcout << buffer;

      // ---- convergence exit (added); never before the BC increment has been applied ----
      if (currentIteration >= 1 &&
          ((useAbsTol && currentNorm < userInputs_cp.absNonLinearTolerance) ||
           (useRelTol && relNorm < userInputs_cp.relNonLinearTolerance)))
        {
          converged = true;
          pcout << "nonlinear solver converged\n";
          break;
        }

      //if not converged, solveLinearSystem Ax=b
      solutionBeforeStep = solution;
      #if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 3)&&(DEAL_II_VERSION_MAJOR==9)))
      computing_timer_cp.enter_section("solve");
      solveLinearSystem(constraints, jacobian, residual, solution, solutionWithGhosts, solutionIncWithGhosts);
      computing_timer_cp.exit_section("solve");
      #else
      computing_timer_cp.enter_subsection("solve");
      solveLinearSystem(constraints, jacobian, residual, solution, solutionWithGhosts, solutionIncWithGhosts);
      computing_timer_cp.leave_subsection("solve");
      #endif
      stepIncrement = solution;
      stepIncrement -= solutionBeforeStep;
      haveStep = true;

      currentIteration++;
    }
    //call updateAfterIteration, if any
    if (userInputs_cp.enableIndentationBCs){

        if ((old_active_set_size!=active_set_size)&&(currentIteration>1))
        {
            freeze_out_iterations = userInputs_cp.freezeActiveSetSteps;
            if (extra_non_linear_iterations > 12)
                pcout << "active set not converging??? continuing...";
            else
            {
              extra_non_linear_iterations += userInputs_cp.freezeActiveSetSteps + 1;
              pcout << "Allowing extra non-linear iteration for Active Set size change "<<old_active_set_size<<"-->"<<active_set_size<<"\n";
            }
        }
        old_active_set_size=active_set_size;
    }

    updateAfterIteration();
  }

  if (!converged && (useAbsTol || useRelTol))
    pcout << "WARNING: nonlinear solver reached the iteration limit without converging\n";

  //update old solution to new converged solution
  oldSolution=solution;
  return true;
}
#include "../../include/multiPhysicsBVP_template_instantiations.h"
