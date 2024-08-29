##
#  CMake script for the PRISMS-MP applications:
##

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

# Find deal.II installation
FIND_PACKAGE(deal.II 8.3.0 REQUIRED
	HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})

DEAL_II_INITIALIZE_CACHED_VARIABLES()

# Set up the debug, release, and run targets
ADD_CUSTOM_TARGET(debug
  COMMAND +env ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND +env ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
  )

ADD_CUSTOM_TARGET(release
  COMMAND +env ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND +env ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )

ADD_CUSTOM_TARGET(run COMMAND main
  COMMENT "Run with ${CMAKE_BUILD_TYPE} configuration"
  )

# Set up the paths to the library files
INCLUDE_DIRECTORIES(include)

# Multiphysics source files
FILE(GLOB multiphysics_sources src/multiPhysicsBVP/*.cc)
FILE(GLOB headers include/*.h)

# Phase-field source files
FILE(GLOB userinputparameters_pf_sources src/userInputParameters_pf/*.cc)
FILE(GLOB solverparameters_fp_sources src/SolverParameters/*.cc)
FILE(GLOB equationdependencyparser_sources src/EquationDependencyParser/*.cc)
FILE(GLOB floodfiller_sources src/FloodFiller/*.cc)
FILE(GLOB orderparameterremapper_sources src/OrderParameterRemapper/*.cc)
FILE(GLOB simplifiedgrainrepresentation_sources src/SimplifiedGrainRepresentation/*.cc)
FILE(GLOB customPDE_sources src/customPDE/*.cc)

# Plasticity source files
FILE(GLOB userInputParameters_cp src/userInputParameters_cp/*.cc)
FILE(GLOB EBVPfunctions src/ellipticBVP/*.cc)
FILE(GLOB CP_masterclass src/materialModels/crystalPlasticity/*.cc)
FILE(GLOB Continuum_masterclass src/materialModels/continuumPlasticity/*.cc)
##FILE(GLOB ContinuumModels_masterclass src/materialModels/continuumPlasticity/models/*.cc)
FILE(GLOB crystalOrientationsIO src/utilityObjects/*.cc)

PROJECT(prisms_mp)
SET(CMAKE_BUILD_TYPE Release)

# Append extra flags for the GNU compiler to suppress some warnings
#if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#        set(DEAL_II_CXX_FLAGS_DEBUG "${DEAL_II_CXX_FLAGS_DEBUG} -Wno-maybe-uninitialized -Wno-deprecated-declarations -Wno-comment -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable")
#        set(DEAL_II_CXX_FLAGS_RELEASE "${DEAL_II_CXX_FLAGS_RELEASE} -Wno-maybe-uninitialized -Wno-deprecated-declarations -Wno-comment -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable")
#endif()

set(DEAL_II_CXX_FLAGS_RELEASE "${DEAL_II_CXX_FLAGS_RELEASE} -Wno-maybe-uninitialized -Wno-deprecated-declarations -Wno-comment -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable")

ADD_LIBRARY(${PROJECT_NAME} ${multiphysics_sources} ${userinputparameters_pf_sources} ${solverparameters_pf_sources} ${equationdependencyparser_sources} ${floodfiller_sources} ${orderparameterremapper_sources} ${simplifiedgrainrepresentation_sources} ${customPDE_sources} ${userInputParameters_cp} ${CP_masterclass} ${Continuum_masterclass} ${crystalOrientationsIO} src/utilities/sortIndexEntryPairList.cc src/variableAttributeLoader/variableAttributeLoader.cc src/utilities/vectorBCFunction.cc src/inputFileReader/inputFileReader.cc src/parallelNucleationList/parallelNucleationList.cc src/variableContainer/variableContainer.cc )
DEAL_II_SETUP_TARGET(${PROJECT_NAME})

PROJECT(prisms_mp_debug)
SET(CMAKE_BUILD_TYPE Debug)

# Append extra flags for the GNU compiler to suppress some warnings
#if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#        set(DEAL_II_CXX_FLAGS_DEBUG "${DEAL_II_CXX_FLAGS_DEBUG} -Wno-maybe-uninitialized -Wno-deprecated-declarations -Wno-comment -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable")
#        set(DEAL_II_CXX_FLAGS_RELEASE "${DEAL_II_CXX_FLAGS_RELEASE} -Wno-maybe-uninitialized -Wno-deprecated-declarations -Wno-comment -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable")
#endif()

set(DEAL_II_CXX_FLAGS_DEBUG "${DEAL_II_CXX_FLAGS_DEBUG} -Wno-maybe-uninitialized -Wno-deprecated-declarations -Wno-comment -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable")

ADD_LIBRARY(${PROJECT_NAME} ${multiphysics_sources} ${userinputparameters_pf_sources} ${solverparameters_pf_sources} ${equationdependencyparser_sources} ${floodfiller_sources} ${orderparameterremapper_sources} ${simplifiedgrainrepresentation_sources} ${customPDE_sources} ${userInputParameters_cp} ${CP_masterclass} ${Continuum_masterclass} ${crystalOrientationsIO} src/utilities/sortIndexEntryPairList.cc src/variableAttributeLoader/variableAttributeLoader.cc src/utilities/vectorBCFunction.cc src/inputFileReader/inputFileReader.cc src/parallelNucleationList/parallelNucleationList.cc src/variableContainer/variableContainer.cc )
DEAL_II_SETUP_TARGET(${PROJECT_NAME})

