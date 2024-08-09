##
#  CMake script for the PRISMS-MP applications:
##

set(TARGET "main")

set(TARGET_SRC
  ${TARGET}.cc
  )

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

# Find deal.II installation
FIND_PACKAGE(deal.II 8.3.0 REQUIRED
	HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})

DEAL_II_INITIALIZE_CACHED_VARIABLES()

# Set up the paths to the library files
INCLUDE_DIRECTORIES(include)

PROJECT(${TARGET})

deal_ii_invoke_autopilot()

# Add an executable or library
# add_executable(test_program main.cpp)

# Optionally, link to any necessary libraries if available
# target_link_libraries(test_program ${LIBRARIES})

