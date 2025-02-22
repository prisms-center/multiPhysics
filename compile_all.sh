#!/bin/bash

# Exit script on error
set -e

# Define the root directory (assuming script is in the root folder)
ROOT_DIR="$(pwd)"
CORE_DIR="$ROOT_DIR"
APP_DIR="$ROOT_DIR/applications/crystalPlasticity"

# Build Core Library
echo "Compiling Core Library..."
pushd "$CORE_DIR"

echo "Removing old compilation files"
rm -rf CMakeFiles Makefile libprisms_mp*  cmake_install.cmake CMakeCache.txt
cmake .
make -j 4  # Use all available CPU cores
popd

# Build Application
echo "Compiling Application..."
pushd "$APP_DIR"

echo "Removing old compilation files"
rm -rf CMakeFiles Makefile results_cp* solution* integratedFields.txt cmake_install.cmake CMakeCache.txt
cmake .
make -j 4
popd

echo "Compilation completed successfully!"

