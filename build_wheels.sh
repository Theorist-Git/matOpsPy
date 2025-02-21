#!/bin/bash
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
        echo "Error: wget is not installed."
        echo "Installing wget"
        yum install -y wget
    fi

    # Download the new matOps.hpp from the specified version
    wget "https://github.com/Theorist-Git/matOps/releases/download/${VERSION}/matOps.hpp" -P include/ || {
        echo "Error: Failed to download matOps.hpp. Check version: ${VERSION}"
        exit 1
    }

    echo "Downloaded matOps.hpp version ${VERSION}"
fi

echo "Building matOpsPy-${VERSION}"

rm -rf build/ dist/ matOpsPy.egg-info/

# Ensure latest packaging tools are installed
for PYBIN in /opt/python/cp38-cp38/bin /opt/python/cp39-cp39/bin /opt/python/cp310-cp310/bin /opt/python/cp311-cp311/bin /opt/python/cp312-cp312/bin /opt/python/cp313-cp313/bin; do
    "${PYBIN}/python" -m pip install --upgrade pip setuptools wheel auditwheel
    "${PYBIN}/python" -m pip install pybind11
    "${PYBIN}/python" setup.py bdist_wheel
done

# Repair wheels to ensure manylinux compliance
auditwheel repair dist/*.whl -w dist/
rm -f dist/*linux_x86_64.whl

echo "matOpsPy-${VERSION} built successfully"
echo "Run twine upload dist/* to upload matOpsPy-${VERSION} to PyPI"
