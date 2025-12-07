# PrexSyn Engine

[![test status](https://img.shields.io/github/actions/workflow/status/luost26/prexsyn-engine/test.yml?style=flat-square&label=test)](https://github.com/luost26/prexsyn-engine/actions/workflows/test.yml)
[![build status](https://img.shields.io/github/actions/workflow/status/luost26/prexsyn-engine/build-conda.yml?style=flat-square)](https://github.com/luost26/prexsyn-engine/actions/workflows/build-conda.yml)
[![version](https://anaconda.org/luost26/prexsyn-engine/badges/version.svg)](https://anaconda.org/luost26/prexsyn-engine)
[![platforms](https://anaconda.org/luost26/prexsyn-engine/badges/platforms.svg)](https://anaconda.org/luost26/prexsyn-engine)
![python](https://img.shields.io/badge/dynamic/yaml?url=https%3A%2F%2Fraw.githubusercontent.com%2Fluost26%2Fprexsyn-engine%2Frefs%2Fheads%2Fmain%2Frattler-recipe%2Fvariants.yaml&query=%24.python&style=flat-square&label=python)

PrexSyn Engine is the C++ backend library for [PrexSyn](https://github.com/luost26/prexsyn). It provides a high-throughput data pipeline that generates synthetic pathways annotated with molecular properties to train PrexSyn models. It also includes a fast detokenizer for reconstructing synthetic pathways and product molecules from model outputs.


## Installation & Usage

PrexSyn Engine can only be installed via Conda. To install, run the following command:

```bash
conda install luost26::prexsyn-engine
```

Please refer to the [documentation](https://prexsyn.readthedocs.io/en/latest/prexsyn-engine/) and  [PrexSyn repository](https://github.com/luost26/prexsyn) for usage instructions.
