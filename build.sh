#!/bin/bash
# build.sh - A simple build script for the pymatops Python module

# Get the pybind11 include flags (this calls the pybind11 CLI)
PYBIND11_FLAGS=$(python3 -m pybind11 --includes)

# Get the Python extension suffix (e.g., .so, .pyd)
EXT_SUFFIX=$(python3-config --extension-suffix)

LIB_NAME="pyMatOps"
SRC_FILE="bindings/bindings.cpp"
OUT_DIR="pyMatOps"

# Ensure the output directory exists
mkdir -p $OUT_DIR

# Compile the bindings into a shared library
g++ -O3 -Wall -shared -std=c++11 -fPIC ${PYBIND11_FLAGS} -Iinclude -o $OUT_DIR/${LIB_NAME}${EXT_SUFFIX} $SRC_FILE

echo "Build complete. You can now import 'pyMatOps' in Python."