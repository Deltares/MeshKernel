# Sepran library
## Information
A Two-Dimensional mesh triangular generator for an area which is defined by a closed polygon.

Library: https://repos.deltares.nl/repos/ds/trunk/guis/plugins_qt/packages/sepran/

## Usage in MeshKernel
1. The source is downloaded from Netlib
2. Only the following files are included in the current repository:
   - triangle.c
   - triangle.h

## Modifications
1. All original Fortran source code was completely translated to standard C++.
2. Although arrays on the C++ side can be indexed with a 0-based indexing, any index values assigned,
   must be in Fortran's 1-based indexing
