#!/bin/bash
# build.sh - A simple build script for the pymatops Python module

set -xe

# Build this version
VERSION=$1

# Check if user provided version number
if [ -z "$VERSION" ]; then
    echo "Error: VERSION argument is missing!"
    echo "Usage: $0 <version>"
    exit 1
fi

# Make include dir if not there
mkdir -p include/

# Ensure wget is installed
if ! command -v wget &>/dev/null; then
    echo "Error: wget is not installed. Please install it and try again."
    exit 1
fi

wget "https://github.com/Theorist-Git/matOps/releases/download/${VERSION}/matOps.hpp" -P include/ || {
    echo "Error: Failed to download matOps.hpp. Check version: ${VERSION}"
    exit 1
}

echo "Building matOps.hpp ${VERSION}"

# Ensure pybind11 and python3-config are available
if ! python3 -m pybind11 --includes &>/dev/null; then
    echo "Error: pybind11 is not installed. Run 'pip install pybind11'."
    exit 1
fi

if ! command -v python3-config &>/dev/null; then
    echo "Error: python3-config is missing. Ensure Python development headers are installed."
    exit 1
fi

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
g++ -O3 -Wall -shared -std=c++11 -fPIC -fvisibility=hidden ${PYBIND11_FLAGS} -Iinclude -o $OUT_DIR/${LIB_NAME}${EXT_SUFFIX} $SRC_FILE

echo "Build complete. You can now import 'pyMatOps' in Python."