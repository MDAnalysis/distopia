Building and testing distopia
#############################

Building distopia requires CMake and a modern C++ compiler
(preferably GCC or Clang). Optionally Ninja can be used to build distopia much
faster than with `make`

To build:

.. code-block:: bash

  mkdir build
  cd build
  cmake ..
  make

or with Ninja:

.. code-block:: bash

  mkdir build 
  cd build
  cmake .. -GNinja
  ninja

To run the tests use `make test` or `ninja test`.


To control the instruction set use **one** the following CMake flags

* `-DDISTOPIA_USE_SSE1` for SSE
* `-DDISTOPIA_USE_SSE2` for SSE2
* `-DDISTOPIA_USE_SSE3` for SSE3
* `-DDISTOPIA_USE_SSSE3` for SSSE3
* `-DDISTOPIA_USE_SSE4_1` for SSE4.1
* `-DDISTOPIA_USE_SSE4_2` for SSE4.2
* `-DDISTOPIA_USE_AVX` for AVX
* `-DDISTOPIA_USE_AVX2` for AVX2

**Or you can let distopia choose for you (default)**


To benchmark methods (from `./build`) run the benchmarks binary using
`./benchmarks`


To assess code coverage build with `cmake -DDISTOPIA_COVERAGE=ON` and use
either make or ninja to build the coverage targets `tests_coverage` or
`test_kernels_coverage`.  You can then view the resulting HTML gcovr coverage
reports in your favourite browser.
