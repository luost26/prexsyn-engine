#!/bin/bash

PREFIX="$1"

if [[ -z "$PREFIX" ]]; then
    echo "Usage: $0 <install_prefix>"
    exit 1
fi

if [[ -d "rdkit" ]]; then
    echo "RDKit directory already exists. Skipping clone."
else
    git clone https://github.com/rdkit/rdkit.git -b Release_2025_09_5 rdkit --depth 1
fi

cd rdkit
rm -rf build
mkdir -p build
cd build

cmake -DCMAKE_BUILD_TYPE=Release \
    -DRDK_INSTALL_INTREE=OFF \
    -DCMAKE_PREFIX_PATH="$PREFIX" \
    -DCMAKE_INSTALL_PREFIX="$PREFIX" \
    -DRDK_BUILD_CPP_TESTS=OFF \
    -DRDK_BUILD_PYTHON_WRAPPERS=OFF \
    -DRDK_BUILD_CHEMDRAW_SUPPORT=OFF \
    -DRDK_BUILD_CAIRO_SUPPORT=OFF \
    -DRDK_BUILD_FREETYPE_SUPPORT=OFF \
    -DRDK_INSTALL_COMIC_FONTS=OFF \
    -DRDK_BUILD_DESCRIPTORS3D=OFF \
    -G Ninja ..

cmake --build .
cmake --install . --prefix $PREFIX
