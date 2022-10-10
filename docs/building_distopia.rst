Building and testing distopia
#############################

Building distopia requires scikit-build, CMake, numpy and a modern C++ compiler
(preferably GCC or Clang).

the most basic build options is to the scikit-build setup as follows

.. code-block:: bash

  python setup.py build 

Or to install 

.. code-block:: bash

  python setup.py install 


Make sure that the resulting shared library is in your LD_LIBRARY_PATH for the library to
work correctly. 

Distopia
--------

The distopia python layer will be built and bind to the shared library regardless
of the options specified for its compilation. See below. 

Libdistopia
-----------

The C++ library component (**libdistopia**) can be built using several possible configurations
using various CMake flags.  There are two main modes, building distopia for use with
a **single instruction set**, or for **dispatch**.

Single instruction set builds
-----------------------------

This is the default mode for distopia, in which a single version of the distopia
functions are built for a single set of SIMD flags. Distopia can attempt to select 
the highest level of SIMD supported on the current computer by setting 
`DISTOPIA_AUTO_SELECT_SIMD=ON`.

Otherwise you can manually select a single instruction set by setting
`DISTOPIA_MANUAL_SELECT_SIMD=ON` and specifying ONE of the following
instruction set flags

* `-DDISTOPIA_USE_SSE1` for SSE
* `-DDISTOPIA_USE_SSE2` for SSE2
* `-DDISTOPIA_USE_SSE3` for SSE3
* `-DDISTOPIA_USE_SSSE3` for SSSE3
* `-DDISTOPIA_USE_SSE4_1` for SSE4.1
* `-DDISTOPIA_USE_SSE4_2` for SSE4.2
* `-DDISTOPIA_USE_AVX` for AVX
* `-DDISTOPIA_USE_AVX2` for AVX2
* `-DDISTOPIA_USE_AVX512` for AVX512


Additionally you can enable aggressive optimisations for the current CPU 
(march=native, mtune=native) by specifying `DISTOPIA_AGGRESSIVE_MARCH=ON`.

An example of an automatically selected single instruction set build is

.. code-block:: bash

  python setup.py install -- -DDISTOPIA_AUTO_SELECT_SIMD=ON -DISTOPIA_AGGRESSIVE_MARCH=ON

An example of a manually selected single instruction set build is

.. code-block:: bash

  python setup.py install -- -DDISTOPIA_MANUAL_SELECT_SIMD=ON -DDISTOPIA_USE_AVX2=ON

Tests and benchmarks
--------------------

**libdistopia** comes with a set of tests and benchmarks that can be enabled with
`DISTOPIA_BUILD_TEST=ON` and `DISTOPIA_BUILD_TIMINGS=ON`.


Dispatch builds
---------------

**libdistopia** can also be built for **dispatch**. This means that the distopia functions 
are compiled multiple times for different instruction sets and packaged into one
shared library. A version of the function is then selected at runtime using runtime
dispatch on available CPU features.  

This is an advanced option if you are cloning the repo and building yourself
but the default version if you download a precompiled binary from conda.
Building for dispatch is incompatible with using `DISTOPIA_AUTO_SELECT_SIMD` or
`DISTOPIA_MANUAL_SELECT_SIMD` or `DISTOPIA_AGGRESSIVE_MARCH=ON`.

Libdistopia can be built for dispatch either by building using
`-DDISTOPIA_DISPATCH_MAX=ON` which builds for **every instruction set** or by 
using `-DDISTOPIA_DISPATCH_MANUAL=ON` and then specifying instruction sets 
using something like `-DDISTOPIA_USE_AVX=ON -DDISTOPIA_USE_AVX2=ON`.

Note that if building using manual dispatch, building for the SSE1 instruction
set is mandatory and is enabled automatically.

An example of building for maximum dispatch

.. code-block:: bash

  python setup.py install -- -DDISTOPIA_DISPATCH_MAX=ON

An example of selecting specific instruction sets for dispatch is

.. code-block:: bash

  python setup.py install -- -DDISTOPIA_DISPATCH_MANUAL=ON -DDISTOPIA_USE_AVX=ON -DDISTOPIA_USE_AVX2=ON