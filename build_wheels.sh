#!/bin/bash
set -e

# Ensure latest packaging tools are installed
for PYBIN in /opt/python/*/bin; do
    "${PYBIN}/python" -m pip install --upgrade pip setuptools wheel auditwheel
    "${PYBIN}/python" -m pip install pybind11
    "${PYBIN}/python" setup.py bdist_wheel
done

# Repair wheels to ensure manylinux compliance
auditwheel repair dist/*.whl -w dist/

