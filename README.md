# PrexSyn Engine

PrexSyn Engine is the C++ backend library for [PrexSyn](https://github.com/luost26/prexsyn). It provides a high-throughput data pipeline that generates synthetic pathways annotated with molecular properties to train PrexSyn models. It also includes a fast detokenizer for reconstructing synthetic pathways and product molecules from model outputs.

## Installation & Usage ([v0](https://github.com/luost26/prexsyn-engine/tree/dev-v0))

[![version](https://anaconda.org/conda-forge/prexsyn-engine/badges/version.svg)](https://anaconda.org/conda-forge/prexsyn-engine)
[![platforms](https://anaconda.org/conda-forge/prexsyn-engine/badges/platforms.svg)](https://anaconda.org/conda-forge/prexsyn-engine)

PrexSyn Engine v0 is automatically installed from conda-forge as a dependency when you install PrexSyn. Please refer to the [documentation](https://prexsyn.readthedocs.io/en/latest/prexsyn-engine/) and  [PrexSyn repository](https://github.com/luost26/prexsyn) for usage instructions.

## 🚧 PrexSyn Engine v1 Work in Progress

The current version ([v0](https://github.com/luost26/prexsyn-engine/tree/dev-v0)) of PrexSyn Engine is dynamically linked to RDKit and Boost to enable interoperability with RDKit's Mol objects on the Python side. 
However, we have found that the restrictions and risks associated with this design choice outweigh the benefits. For example, dynamic linking limits installation to Conda only and prevents distribution via PyPI. Moreover, future changes in RDKit API/ABI (which are very likely) could break the compatibility and cause runtime errors that are hard to debug.

Therefore, we are working on a new version of PrexSyn Engine that is statically linked to RDKit. This new version will be distributed as a PyPI wheel package. It will be free from the risk of breaking changes in upstream libraries. The development branch for PrexSyn Engine v1 is [here](https://github.com/luost26/prexsyn-engine/tree/dev-v1).
