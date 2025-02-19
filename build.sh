#!/bin/bash
# build.sh - A simple build script for the matOpsPy Python module
# sudo docker run --rm -it -v "$(pwd)":/io quay.io/pypa/manylinux2014_x86_64 /bin/bash

set -e

# Default: do not install new matOps.hpp
INSTALL_NEW=0

# Check for the flag --install-new as the first argument
if [ "$1" == "--install-new" ]; then
    INSTALL_NEW=1
    shift  # Remove the flag from the arguments so that $1 becomes the version
fi

if [ $INSTALL_NEW -eq 1 ]; then
    # Check if user provided a version number
    if [ -z "$1" ]; then
        echo "Error: VERSION argument is missing!"
        echo "Usage: $0 --install-new <version>"
        exit 1
    fi

    VERSION=$1

    # Create include directory if not present
    mkdir -p include/

    # Remove existing matOps.hpp if it exists
    if [ -f "include/matOps.hpp" ]; then
        echo "Removing existing matOps.hpp"
        rm -f include/matOps.hpp
    fi

    # Ensure wget is installed
    if ! command -v wget &>/dev/null; then
        echo "Error: wget is not installed. Please install it and try again."
        exit 1
    fi

    # Download the new matOps.hpp from the specified version
    wget "https://github.com/Theorist-Git/matOps/releases/download/${VERSION}/matOps.hpp" -P include/ || {
        echo "Error: Failed to download matOps.hpp. Check version: ${VERSION}"
        exit 1
    }

    echo "Downloaded matOps.hpp version ${VERSION}"
fi

echo "Building matOpsPy"

# Ensure pybind11 is installed and python3-config is available
if ! python3 -m pybind11 --includes &>/dev/null; then
    echo "Error: pybind11 is not installed. Run 'pip install pybind11'."
    exit 1
fi

if ! command -v python3-config &>/dev/null; then
    echo "Error: python3-config is missing. Ensure Python development headers are installed."
    exit 1
fi

# Get pybind11 include flags and the Python extension suffix
PYBIND11_FLAGS=$(python3 -m pybind11 --includes)
EXT_SUFFIX=$(python3-config --extension-suffix)

LIB_NAME="matOpsPy"
SRC_FILE="bindings/bindings.cpp"
OUT_DIR="matOpsPy"

# Ensure the output directory exists
mkdir -p "$OUT_DIR"

# Compile the bindings into a shared library
g++ -O3 -Wall -shared -std=c++11 -fPIC -fvisibility=hidden ${PYBIND11_FLAGS} -Iinclude -o "$OUT_DIR/${LIB_NAME}${EXT_SUFFIX}" "$SRC_FILE"

echo "Build complete. Run 'pip install -e .' to install matOpsPy into your environment."
