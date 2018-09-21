# dS_instanton_solver
Library of functions and classes for solving for bounce solutions in de Sitter space, including back-reaction, at arbitrary precision.

Currently this software is in an early stage of development, although the code is useable. Documentation is sparse, but will be added to in the future.

----------------------------------------------------------------------------------------------------------------------------------
COMPILATION

A makefile is included to compile the code into a library with gcc by default. Simply navigate to the folder where the code is stored and run 'make'.

The code has also been successfully compiled on Windows - more support and instructions on this will be provided in the future.

Dependencies:

BOOST library. The code uses a modified version of the boost odeint library, included in the folder odeint_with_event_detection. When compiling, it is vital that this folder be included BEFORE the main boost libraries, so that it overides the standard routine.

NB - this modified version does not work with the most recent version of the boost libraries. It has been tested and runs with version 1.60.0, which can be downloaded at https://www.boost.org/ or https://sourceforge.net/projects/boost/files/boost/1.60.0/. The makefile should be edited to refer to the directory where the boost libraries are installed. Future releases will aim to address this incompatibility.

MPFR - If arbitrary precision numbers are required, then the code makes use of the MPFR library available at https://www.mpfr.org/ by default. It has been tested with version 4.0.1, but should work with later versions. Note that other implementations of arbitrary precision arithmetic can also be used, although these may require editing multi_precision_definitions.h and/or multi_precision_definitions.cpp (especially the convert_type functions for converting between double and arbitrary precision numbers).

If only double precision is needed, the library can be compiled without MPFR.

----------------------------------------------------------------------------------------------------------------------------------
USAGE

The code provides several functions for solving the O(4) symmetric euclidean bounce equations for a scalar field both in flat space and in curved space, including back-reaction. Currently this compiles to a library of functions that can be used in other code. A python and MATLAB wrapper are aimed to be released in the future.

--------------------------------------------------------------------------------------------------------------------------------
ACKNOWLEDGEMENTS

The author would like to thank the developers of the boost libraries, a modified version of which forms the core of the ode-solver used in the code.
