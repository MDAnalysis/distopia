distopia
--------

Faster distance calculations for the year 2020 🚀 🚀 🚀 

To build:
```
  mkdir build
  cd build
  cmake ..
  make
 ```
 
To generate test coordinate data:
```python
  python generate_coords.py 1000000
```
Will generate 1000000 random coordinates.

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




To benchmark methods (from `./build`):
```
  ./timings ../data.txt 1000
```
Will run through the timings of available functions (depends on the instruction set used), averaging over 1000 runs.
 
