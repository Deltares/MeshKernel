# Triangle library
## Information
A Two-Dimensional Quality Mesh Generator and Delaunay Triangulator.

Source: www.cs.cmu.edu/~quake/triangle.html  
Library: http://www.netlib.org/voronoi/triangle.zip 

## Usage in MeshKernel
1. The source is downloaded from Netlib
2. Only the following files are included in the current repository:
   - triangle.c
   - triangle.h

## Modifications
The following modifications were necessary:
1. All variables declared with data type `unsigned long` were changed to `unsigned long long`
2. Literals `0l`, `1l`, `1l`, ... were changed to `0ull`, `1ull`, `1ull`
3. `printf` format specifier `%lx` were changed to `%llx`