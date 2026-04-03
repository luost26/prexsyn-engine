#!/bin/bash

PREFIX="$1"

if [[ -z "$PREFIX" ]]; then
    echo "Usage: $0 <install_prefix>"
    exit 1
fi

mkdir -p "$PREFIX"

LIBRARIES=(serialization iostreams program_options json stacktrace format graph rational flyweight math property_tree crc multiprecision assign multi_array)

if [[ -d "boost" ]]; then
    echo "Boost directory already exists. Skipping clone."
else
    git clone https://github.com/boostorg/boost.git -b boost-1.86.0 boost --depth 1
fi
cd boost

git submodule update --depth 1 -q --init tools/boostdep
for lib in "${LIBRARIES[@]}"; do
    git submodule update --depth 1 -q --init "libs/$lib"
    python tools/boostdep/depinst/depinst.py -X test -g "--depth 1" $lib
done

if [[ -n "$BOOST_TOOLSET" ]]; then
    ARGS="--with-toolset=$BOOST_TOOLSET"
else
    ARGS=""
fi

./bootstrap.sh --prefix=$PREFIX $ARGS #--with-libraries=$(IFS=,; echo "${LIBRARIES[*]}")
# ./b2 --clean-all
./b2 link=static cxxflags=-fPIC cflags=-fPIC
./b2 install
