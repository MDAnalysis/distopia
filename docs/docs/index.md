# Distopia

Distopia is a package to rapidly calculate distances, angles and dihedrals under periodic boundary conditions in single and double precision. Explicit SIMD vectorisation allows awesome speedups over autovectorised code (up to 10x). The distopia package consists of consists of the python layer (distopia) and a C++ library (libdistopia) that does the heavy lifting.

Distopia can be used out of the box by building the library using the instructions in [Building and Testing](building_and_testing.md). Examples of use are given in [Examples](examples.md).

## Current Status

Distopia is currently under active development and should be considered alpha This means the API is liable to change rapidly without warning.

## OS and compiler support

We currently support x86_64 and ARM with the Clang and GCC family of compilers. Windows and MSVC probably work, but are not tested in CI, so use at your own risk. 


## Participating

Ask questions on the mdnalysis-discussion mailing list and be part of the conversation. You can also join the MDAnalysis Discord Server to talk with other users and developers. In order to join our Discord server, use the following [invitation](https://discord.com/invite/fXTSfDJyxE). Please report bugs or enhancement requests through the [Issue Tracker](https://github.com/MDAnalysis/distopia/issues). Distopia is open source and welcomes your contributions. Fork the repository on GitHub and submit a pull request!

