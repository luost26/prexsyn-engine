#!/bin/bash

if [[ -n "$CONDA_PREFIX" ]]; then
    echo "Installing Boost to $CONDA_PREFIX"
else
    echo "Error: CONDA_PREFIX is not set. Please activate a conda environment."
    exit 1
fi

if [[ -d "boost" ]]; then
    echo "Boost directory already exists. Skipping clone."
else
    wget https://archives.boost.io/release/1.86.0/source/boost_1_86_0.tar.gz
    tar -xzf boost_1_86_0.tar.gz
    mv boost_1_86_0 boost
fi
cd boost


./bootstrap.sh --prefix=$CONDA_PREFIX --with-toolset=clang
./b2
./b2 install
