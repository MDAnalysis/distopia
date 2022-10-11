.. distopia documentation master file, created by
   sphinx-quickstart on Wed Aug 25 23:17:34 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to distopia's documentation!
====================================

:Release: |release|
:Date: |today|

**Distopia** is a package to rapidly calculate distances
under periodic boundary conditions. Explicit SIMD vectorisation allows awesome
speedups over autovectorised code (up to 10x). The distopia package consists of
consists of the python layer (**distopia**) and a C++ library (**libdistopia**)
that does the heavy lifting. 

Distopia can be used out of the box by building the library using the
instructions in :ref:`Building and testing distopia`.

**Current Status:**

Distopia is currently under active development and should be considered
**pre 0.1.0**. This means the API is liable to change rapidly without warning.

**OS and compiler support:**

We currently support x86 linux and mac-os machines, and the Clang and GCC
family of compilers. Windows and MSVC support along with explicitly vectorised
ARM code is on the TODO list. x86_64 SIMD widths up to AVX512 are supported.

**Participating:**

Ask questions on the mdnalysis-discussion mailing list and be part of the
conversation. You can also join the MDAnalysis Discord Server to talk with
other users and developers. (In order to join our Discord server, use the
following invitation_. Please report bugs or
enhancement requests through the Issue Tracker. Distopia is open source and
welcomes your contributions. Fork the repository on GitHub and submit a pull
request!

.. _invitation: https://discord.gg/fXTSfDJyxE


.. toctree::
   :maxdepth: 4
   :caption: Contents:
   
   ./building_distopia.rst
   ./api/distopia.rst
   ./api/libdistopia.rst
   





Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
