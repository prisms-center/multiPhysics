//#include "../../include/ParseCommandLineOpts.h"
#include "../../include/inputFileReader.h"
#include "../../src/variableAttributeLoader/variableAttributeLoader.cc"

//tension BVP
//general headers
#include <fstream>
#include <sstream>
#include <iostream>
using namespace std;

#include "../../include/crystalPlasticity.h"

// Header file for postprocessing that may or may not exist
//#ifdef POSTPROCESS_FILE_EXISTS
#include "../../src/customPDE/postprocess.cc"
//#else
//void variableAttributeLoader::loadPostProcessorVariableAttributes(){}
//#endif

//main
int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  try
    {
      deallog.depth_console(0);

      ParameterHandler parameter_handler;

      std::list<std::string> args;
      for (int i=1; i<argc; ++i) args.push_back (argv[i]);

      if (args.size() == 0){
        std::cerr<<"Provide name of a parameter file."<<std::endl;
        exit (1);
      }

      const std::string parameter_file_cp = args.front ();
      const std::string parameter_file_pf = "parameters_pf.prm";

      //Read Crystal Plasticity parameters
      userInputParameters_cp userInputs_cp(parameter_file_cp,parameter_handler);

      //Read PhaseField Parameters
      variableAttributeLoader variable_attributes;
      inputFileReader input_file_reader(parameter_file_pf,variable_attributes);
      userInputParameters_pf<3> userInputs_pf(input_file_reader,input_file_reader.parameter_handler,variable_attributes);

      crystalPlasticity<3> problem(userInputs_pf, userInputs_cp);

      //reading materials atlas files
      
      if(!userInputs_cp.readExternalMesh)
      {
        problem.orientations.loadOrientations(userInputs_cp.grainIDFile,
  					    userInputs_cp.headerLinesGrainIDFile,
  					    userInputs_cp.grainOrientationsFile,
  					    userInputs_cp.numPts,
  					    userInputs_cp.span);
      }
      
      problem.orientations.loadOrientationVector(userInputs_cp.grainOrientationsFile, userInputs_cp.enableMultiphase, userInputs_cp.additionalVoxelInfo);

      problem.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    }

  return 0;
}
