#!/bin/bash

PREFIX="$1"

if [[ -z "$PREFIX" ]]; then
    echo "Usage: $0 <install_prefix>"
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

if [[ -n "$CC" ]]; then
    ARGS="--with-toolset=$(basename $CC)"
else
    ARGS=""
fi

./bootstrap.sh --prefix=$PREFIX --with-libraries=serialization,iostreams,program_options $ARGS
./b2 --clean-all
./b2 link=static cxxflags=-fPIC cflags=-fPIC
./b2 install
