# Introduction

The MeshKernel library is a C++ dynamic library that performs generation and manipulations of 2D grids.
The code was re-written from a FORTRAN code base to C++ for the following reasons:

-   Split the monolithic code base.   

-   Separating concerns.

-   Introduce encapsulation: from global data and global methods to object oriented design.

-   Introduce unit testing.

-   Enabling the visualization of the meshes during creation
    (interactivity).

-   Simplify the connection with C like languages, such C#. C++ and
    C# have the same array memory layout and pointers can be passed
    seamlessly, without the need for pointer conversion.
