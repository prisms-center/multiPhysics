//#include "../../include/ParseCommandLineOpts.h"
#include "../../include/inputFileReader.h"
#include "../../src/variableAttributeLoader/variableAttributeLoader.cc"

//tension BVP
//general headers
#include <fstream>
#include <sstream>
#include <iostream>
//using namespace std;

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

    if (argc != 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <CPFE parameters file> <PF parameters file>" << std::endl;
        return 1;
    }

    const std::string parameter_file_cp = argv[1];
    const std::string parameter_file_pf = argv[2];

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
