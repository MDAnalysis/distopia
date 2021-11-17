distopia design
---------------

This is a document to outline the rough design of distopia file by file.
The general design philosophy is to abstract away the use of intrinsics as much as possible and use generic algorithms.

The following is a general description of the way the API works.

distopia.h is the public header, include this to use the library. The functions in this file calls dispatchers in the respective file prefixed with **calc_** eg calc_bonds.cpp calc_angles.cpp. Each of these files also defines an **Inner** function that does the loop and memory access pattern.

 As the kernels are much faster in SOA format memory access and AOS->SOA transforms handled by a class that handles 3*SIMD_WIDTH data packets called **vector_triple**. The files responsible for this and its corresponding operator overloads are prefixed **vector_triple**.  The actual kernels used in **calc_** files are in **kernels.h** . Additional helper functions for the kernels are in **ops.h**. The setup to handle PBCs is in **box.h**. 
 
The operations for these key parts are then defined in terms of a series of required operations and operator overloads in **basemath.h**. For new SIMD types these operations must all be defined along with macros and swizzles that abstract away the basic data operations. For x86, this is all handled in **/src/lib/x86/**. An implementation for ARM would have to overload all of the operations provided by the low level API and the correct inclusion handled by an intermediary header. 

Other files are mostly for detecting the compiled-upon architecture and providing compiler hints. 


File by file
------------

src/
    -> lib/   the distopia library
    -> compare/   kernels from MDAnalysis and MDtraj to compare against

diving deeper into lib:

src/lib
    -> include/
        -> distopia.h the distopia public header
    -> src/ the actual source code
        -> arch_config.h
        -> atan.cpp
        -> atan.h
        -> basemath.h
        -> box.h
        -> calc_angles.h
        -> calc_bonds.h
        -> compiler_hints.h
        -> distopia_type_traits.h
        -> kernels.h
        -> ops.h
        -> simd_config.h
        -> vector_triple.h
        -> vector_triple_basemath.h
    -> tests
        -> googletest the googletest source code
        -> googlebench the googlebench source code
        -> benchmark.cpp the benchmarks for distopia and compare kernels
        -> test_kernels.cpp the tests for distopia kernels, tests that distopia matches the vanilla     implementation
        -> test_mda_match.cpp tests that distopia matches the mda implementation
        -> tests.cpp unit tests for individual components of distopia 


