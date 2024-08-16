//solve method for multiPhysicsBVP class
#include "../../include/multiPhysicsBVP.h"

//loop over increments and solve each increment
template <int dim, int degree>
void MultiPhysicsBVP<dim,degree>::solve_cp(){
  pcout << "begin solve...\n\n";
  bool success;
  //load increments
  unsigned int successiveIncs=0;

  if(userInputs_cp.enableAdaptiveTimeStepping){
    for (;totalLoadFactor<totalIncrements;){
      ++currentIncrement_cp;
      loadFactorSetByModel=std::min(loadFactorSetByModel, totalIncrements-totalLoadFactor);
      pcout << "\nincrement: "  << currentIncrement_cp << std::endl;
      char buffer[100];
      sprintf(buffer, "current load factor: %12.6e\ntotal load factor:   %12.6e\n", loadFactorSetByModel, totalLoadFactor);
      pcout << buffer;

      //call updateBeforeIncrement, if any
      updateBeforeIncrement();

      if (!userInputs_cp.flagTaylorModel){
        //solve time increment
        success=solveNonLinearSystem();
      }

      //call updateAfterIncrement, if any
      if ((success)||(userInputs_cp.flagTaylorModel)){
        updateAfterIncrement();

        //update totalLoadFactor
        totalLoadFactor+=loadFactorSetByModel;

        //increase loadFactorSetByModel, if succesiveIncForIncreasingTimeStep satisfied.
        successiveIncs++;

        if (successiveIncs>=userInputs_cp.succesiveIncForIncreasingTimeStep){
          loadFactorSetByModel*=userInputs_cp.adaptiveLoadIncreaseFactor;
          char buffer1[100];
          sprintf(buffer1, "current increment increased. Restarting increment with loadFactorSetByModel: %12.6e\n", loadFactorSetByModel);
          pcout << buffer1;
        }
        #if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 3)&&(DEAL_II_VERSION_MAJOR==9)))
        computing_timer_cp.enter_section("postprocess");
        #else
	    computing_timer_cp.enter_subsection("postprocess");
        #endif

        if (currentIncrement_cp%userInputs_cp.skipOutputSteps==0)
        if (userInputs_cp.writeOutput) output();

        #if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 3)&&(DEAL_II_VERSION_MAJOR==9)))
        computing_timer_cp.exit_section("postprocess");
        #else
	    computing_timer_cp.leave_subsection("postprocess");
        #endif

      }
      else{
        successiveIncs=0;
      }
    }
    char buffer[100];
    sprintf(buffer, "\nfinal load factor  : %12.6e\n", totalLoadFactor);
    pcout << buffer;
  }
  else
    for (;currentIncrement_cp<totalIncrements; ++currentIncrement_cp){
    pcout << "\nincrement: "  << currentIncrement_cp << std::endl;
    if (userInputs_cp.enableIndentationBCs){
        MultiPhysicsBVP<dim,degree>::updateBeforeIncrement();
        if (!userInputs_cp.continuum_Isotropic)
            updateBeforeIncrement();
    }
    else{
        updateBeforeIncrement();
    }
    //call updateBeforeIncrement, if any


    if (!userInputs_cp.flagTaylorModel){
      //solve time increment
      success=solveNonLinearSystem();
    }

    //call updateAfterIncrement, if any
    if ((success)||(userInputs_cp.flagTaylorModel)){
      updateAfterIncrement();

      //update totalLoadFactor
      totalLoadFactor+=loadFactorSetByModel;

      //increase loadFactorSetByModel, if succesiveIncForIncreasingTimeStep satisfied.
      successiveIncs++;
      //output results to file
      //
      #if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 3)&&(DEAL_II_VERSION_MAJOR==9)))
      computing_timer_cp.enter_section("postprocess");
      #else
      computing_timer_cp.enter_subsection("postprocess");
      #endif

      //////////////////////TabularOutput Start///////////////
      std::vector<unsigned int> tabularTimeInputIncInt;
      std::vector<double> tabularTimeInputInc;
      if (userInputs_cp.tabularOutput){

        tabularTimeInputInc=userInputs_cp.tabularTimeOutput;
        for(unsigned int i=0;i<userInputs_cp.tabularTimeOutput.size();i++){
          tabularTimeInputInc[i]=tabularTimeInputInc[i]/delT;
        }
        tabularTimeInputIncInt.resize(userInputs_cp.tabularTimeOutput.size(),0);
        ///Converting to an integer always rounds down, even if the fraction part is 0.99999999.
        //Hence, I add 0.1 to make sure we always get the correct integer.
        for(unsigned int i=0;i<userInputs_cp.tabularTimeOutput.size();i++){
          tabularTimeInputIncInt[i]=int(tabularTimeInputInc[i]+0.1);
        }
      }
      //////////////////////TabularOutput Finish///////////////
      if (((!userInputs_cp.tabularOutput)&&((currentIncrement_cp+1)%userInputs_cp.skipOutputSteps==0))||((userInputs_cp.tabularOutput)&& (std::count(tabularTimeInputIncInt.begin(), tabularTimeInputIncInt.end(), (currentIncrement_cp+1))==1))){
        if (userInputs_cp.writeOutput) output();
      }

      #if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 3)&&(DEAL_II_VERSION_MAJOR==9)))
      computing_timer_cp.exit_section("postprocess");
      #else
      computing_timer_cp.leave_subsection("postprocess");
      #endif
    }
    else{
      successiveIncs=0;
    }
  }
}
#include "../../include/multiPhysicsBVP_template_instantiations.h"
