//#include "../../include/ParseCommandLineOpts.h"
#include "../../include/inputFileReader.h"
#include "../../src/variableAttributeLoader/variableAttributeLoader.cc"

//tension BVP
//general headers
#include <fstream>
#include <sstream>
#include <iostream>

#include "../../include/crystalPlasticity.h"

// Header files for user-defined PF application functions
#include "customPDE.h"
#include "ICs_and_BCs.cc"
#include "equations.cc"
#include "postprocess.cc"

//main
int main (int argc, char **argv)
{
    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
    try
      {
        deallog.depth_console(0);

        ParameterHandler parameter_handler;

        if (argc != 2)
          {
            std::cerr << "Usage: " << argv[0]
                      << " <parameters file>" << std::endl;
            return 1;
          }

        const std::string parameter_file = argv[1];

        // Read PhaseField Parameters (declares PF parameters)
        variableAttributeLoader variable_attributes;
        inputFileReader input_file_reader(parameter_file,variable_attributes,parameter_handler);

        // Read Crystal Plasticity parameters (declares CP parameters and parses file)
        userInputParameters_cp userInputs_cp(parameter_file,parameter_handler);

        userInputParameters_pf<3> userInputs_pf(input_file_reader,input_file_reader.parameter_handler,variable_attributes);

        customPDE<3,1> pf_problem(userInputs_pf, userInputs_cp);
        crystalPlasticity<3> problem(userInputs_pf, userInputs_cp, pf_problem);

        // reading materials atlas files
      
        if(!userInputs_cp.readExternalMesh)
          {
            problem.orientations.loadOrientations(userInputs_cp.grainIDFile,
                              userInputs_cp.headerLinesGrainIDFile,
                              userInputs_cp.grainOrientationsFile,
                              userInputs_cp.numPts,
                              userInputs_cp.span);
          }
      
        problem.orientations.loadOrientationVector(userInputs_cp.grainOrientationsFile, userInputs_cp.enableMultiphase, userInputs_cp.additionalVoxelInfo);

        // Store a reference to the CPFE orientations to customPDE so they can be used in the phase-field model
        pf_problem.cpfe_orientations = &problem.orientations;

        problem.run();
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
