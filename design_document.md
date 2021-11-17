distopia design
---------------

This is a document to outline the rough design of distopia, including file by file.
The general design philosophy is to abstract away the use of intrinsics as much as possible and use generic algorithms.

The following is a general description of the way the API works.

distopia.h is the public header, include this to use the library. The functions in this file calls dispatchers in the respective file prefixed with **calc_** eg calc_bonds.cpp calc_angles.cpp. Each of these files also defines an **Inner** function that does the loop and memory access pattern.

 As the kernels are much faster in SOA format, memory access and AOS*SOA transforms are handled by a class that contains 3*SIMD_WIDTH data packets called **vector_triple**. The files responsible for this and its corresponding operator overloads are prefixed **vector_triple**.  The actual kernels used in **calc_** files are in **kernels.h** . Additional helper functions for the kernels are in **ops.h**. The setup to handle PBCs is in **box.h**. 
 
The operations for these key parts are then defined in terms of a series of required operations and operator overloads in **basemath.h**. For new SIMD types these operations must all be defined along with **macros** and **swizzles** that abstract away the basic data operations. For x86, this is all handled in **/src/lib/x86/**. An implementation for ARM would have to overload all of the operations provided by the low level API and the correct inclusion handled by an intermediary header. 

Other files are mostly for detecting the compiled-upon architecture and providing compiler hints and are hopefully self explanatory. 


File by file
------------

* src/
  * lib/   the distopia library
  * compare/   kernels from MDAnalysis and MDtraj to compare against

diving deeper into lib:

* src/lib
  * include/
    * distopia.h the distopia public header
  * src/ the actual source code
    * arch_config.h detects ISA and SIMD features available on current system
    * atan.cpp Jakub's implementation of an x86 atan based on https://opensource.apple.com/source/Libm/Libm-287.1/Source/Intel/atanf.s
    * atan.h header for atan2
    * basemath.h defines essential operations for floating point types
    * box.h defines class for simulation boxes
    * calc_angles.h dispatcher and inner loop for calculation of angles 
    * calc_bonds.h dispatcher and inner loop for calculation of distances 
    * compiler_hints.h compiler hints
    * distopia_type_traits.h custom type traits for vector types 
    * kernels.h The core operations for calculating distances, angles, dihedrals
    * ops.h helper routines for kernels
    * simd_config.h Allows setting the SIMD level with CMake, currently not functional
    * vector_triple.h defines a class for easy processing of coordinate data and AOS-SOA acces to transforms
    * vector_triple_basemath.h operator overloads for vector triple types
    * x86/ specific code to implement basemath and AOS-SOA for x86 vector types (currently up to 256 wide)
      * x86_basemath.h vector overloads for basemath.h
      * x86_swizzle.h AOS-SOA transforms and related shuffles
      * x86_tgintrin.h Macros that map the x86 intrinsics to a unified interface
      * x86_vector_operators.h operator overloads for x86 vectors

    * tests/
      * googletest/ the googletest source code
      * googlebench/ the googlebench source code
      * benchmark.cpp the benchmarks for distopia and compare kernels
      * test_kernels.cpp the tests for distopia kernels, tests that distopia matches the vanilla     implementation
      * test_mda_match.cpp tests that distopia matches the mda implementation
      * tests.cpp unit tests for individual components of distopia 



