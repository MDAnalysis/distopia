distopia
--------

## Faster distance calculations using SIMD intrinsics


Distopia is a package to rapidly calculate distances, angles and dihedrals under periodic boundary conditions in single and double precision. Explicit SIMD vectorisation allows awesome speedups over autovectorised code (up to 10x). The distopia package consists of consists of the python layer (distopia) and a C++ library (libdistopia) that does the heavy lifting.

Distopia can be used out of the box by building the library using the instructions in Building and Testing. Examples of use are given in Examples.


## Getting started

Get started by reading the documentation here: https://www.mdanalysis.org/distopia

## Current Status
Distopia is currently under active development and should be considered alpha This means the API is liable to change rapidly without warning.

## OS and compiler support

We currently support x86_64 and ARM with the Clang and GCC family of compilers. Windows and MSVC may work, but are not tested in CI, so use at your own risk.

## Dependencies

We use several frameworks developed by Google including [Highway](https://github.com/google/highway), [GoogleTest](https://github.com/google/googletest) and [GoogleBench](https://github.com/google/benchmark). Thank you google!

## Participating
Ask questions on the mdanalysis-discussion mailing list and be part of the conversation. You can also join the MDAnalysis Discord Server to talk with other users and developers. In order to join our Discord server, use the following [invitation](https://discord.com/invite/fXTSfDJyxE). Please report bugs or enhancement requests through the Issue Tracker. Distopia is open source and welcomes your contributions. Fork the repository on GitHub and submit a pull request!
