# MeshKernel

Deltares C++ library for creating and editing 2D unstructured and curvilinear meshes, suitable for the D-Flow Flexible Mesh simulator.

The library is separated in an API namespace (MeshKernelApi), used for communication with the client and a backend namespace (MeshKernel), where the algorithms are implemented. 
The API namespace contains several structures used as parameters for the API methods (see API usage section). 
These structures must be mirrored in the client application and filled with appropriate values.

## Build

### Windows

MeshKernel depends on boost. On windows precompiled boost binaries (with MSVC compiler) can be downloaded here:

https://sourceforge.net/projects/boost/files/boost-binaries/ 

Once installed, modify boost environmental variables accordingly. For example:

BOOST_INCLUDEDIR=C:\Apps\boost_1_70_0

BOOST_LIBRARYDIR=C:\Apps\boost_1_70_0\lib64-msvc-14.1

Install CMAKE (https://cmake.org/download/) and execute bootstrap.bat. The visual studio solution is produced inside the build folder.


## Examples

1. Creating a triangular mesh inside a polygon

In this example a mesh is created by discretizing the polygon perimeter with the desired edge length

![alt tag](doc/figures/TriangularMeshInPolygon.jpg)

2. Mesh orthogonalization

Finite volume staggered flow solvers require the mesh to be as much orthogonal as possible. 
MeshKernel provides an algorithm to adapt the mesh and achieve a good balance between mesh orthogonality and smothness.

![alt tag](doc/figures/MeshOrthogonalization.jpg)

3. Curvilinear mesh generation

Curvilinear meshes for rivers can be generated using splines.

![alt tag](doc/figures/OrthogonalCurvilinearGrid.jpg)

4. Mesh refinement

A mesh can be refined in areas based on samples or polygon selections 

![alt tag](doc/figures/GridRefinement.jpg)


## API usage

Setting a triangular mesh and moving its 2nd node to position 1.0, 3.0:
```c++
// Create a new mesh entry into MeshKernel library
int meshKernelId;
int state = mkernel_new_mesh(meshKernelId);

// Mesh nodes and edges
std::vector<double> nodex{ 0.0, 3.0, 1.5};
std::vector<double> nodey{ 0.0, 0.0, 3.0};
std::vector<int> edge_nodes{0, 1, 1, 2, 2, 0};

// The MeshGeometry communication structures
MeshGeometryDimensions meshGeometryDimensions;
meshGeometryDimensions.numnode = 3;
meshGeometryDimensions.numedge = 3;

MeshGeometry meshGeometry;
meshGeometry.nodex = &nodex[0];
meshGeometry.nodey = &nodey[0];
meshGeometry.edge_nodes = &edge_nodes[0];

// if lat, lon coordinate isGeographic true, otherwise false
bool isGeographic = false;

// Set the mesh into the mesh entry created before
state = mkernel_set_state(meshKernelId, meshGeometryDimensions, meshGeometry, bool isGeographic);

// The new position
std::vector<double> newPositionX{ 1.0};
std::vector<double> newPositionY{ 3.0};

GeometryListNative geometryListIn;
geometryListIn.xCoordinates = &newPositionX[0];
geometryListIn.yCoordinates = &newPositionY[0];

// Move the second node to the new position  
int nodeIndex = 2; 
state = mkernel_move_node(meshKernelId, geometryListIn, nodeIndex);

// Deallocate the mesh entry
state = mkernel_deallocate_state(int meshKernelId);
```
