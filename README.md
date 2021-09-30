distopia
--------

Faster distance calculations for the year 2020 🚀 🚀 🚀 

**Documentation can be found at https://www.mdanalysis.org/distopia/**


To build:
```
  mkdir build
  cd build
  cmake ..
  make
 ```
or with Ninja:

```
  mkdir build 
  cd build
  cmake .. -GNinja
  ninja
```
 
To run the tests:
```
make test
```
or  with Ninja:
```
ninja test
```

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

To benchmark methods (from `./build`) run the benchmarks binary:
```
  ./benchmarks
```
Distopia's documentation is made with breathe using a combination of Doxygen and Sphinx
To build the docs, go to the doc directory and:

```
doxygen
make html
```
then point your browser to the generated documentation.

To assess code coverage build with `cmake -DDISTOPIA_COVERAGE=ON` and use either make or ninja to build the coverage targets `tests_coverage` or `test_kernels_coverage`.  You can then view the resulting HTML gcovr coverage reports in your favourite browser.

