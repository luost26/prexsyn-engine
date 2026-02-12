#!/bin/bash

if [[ -n "$CONDA_PREFIX" ]]; then
    echo "Installing Boost to $CONDA_PREFIX"
else
    echo "Error: CONDA_PREFIX is not set. Please activate a conda environment."
    exit 1
fi

if [[ -d "rdkit" ]]; then
    echo "RDKit directory already exists. Skipping clone."
else
    git clone https://github.com/rdkit/rdkit.git -b Release_2025_09_5 rdkit --depth 1
fi

cd rdkit
mkdir -p build
cd build

cmake -DCMAKE_BUILD_TYPE=Release \
    -DRDK_INSTALL_INTREE=OFF \
    -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX" \
    -DRDK_BUILD_CPP_TESTS=OFF \
    -DRDK_BUILD_PYTHON_WRAPPERS=OFF \
    -G Ninja ..

cmake --build .
cmake --install . --prefix $CONDA_PREFIX
