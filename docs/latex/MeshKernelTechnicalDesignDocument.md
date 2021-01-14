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

This document describes the design of MeshKernel
(chapter [2](#chap:design){reference-type="ref"
reference="chap:design"}), and then the responsibilities of the main
classes.

# MeshKernel design {#chap:design}

The choices made for the representation of the geometric entities have
a strong influence on library design and performance. In MeshKernel an
unstructured mesh is uniquely defined by two entities:

-   The nodes vector: represented using `std::vector<Point>`, where
    `Point` is a structure containing the 2D coordinates of a point (cartesian or
    spherical).

-   The edges vector: represented using an `std::vector<std::
    pair<size_t,size_t>>`, containing the start and the end indices of the
    edges in the node vector described above.

All other mesh properties are computed from these two entities, such as
the face nodes, the face edges, the faces mass centers, and the faces
circumcenters. See [mesh](#chap:mesh){reference-type="ref"
reference="chap:mesh"} for some more details.

MeshKernel has an API namespace (`meshkernelapi`)
and a back-end namespace (`meshkernel`), where the classes implementing the algorithms 
are included, as shown in Figure [2.1](#fig:classDiagram){reference-type="ref"
reference="fig:classDiagram"}. `meshkernelapi` contains the library
API methods and several structures are used for communicating with the clients. 
These structures must be replicated in the client applications 
and filled with appropriate values.

An example of library execution is shown in figure
[2.2](#fig:sequenceDiagram){reference-type="ref"
reference="fig:sequenceDiagram"}. When the client application creates a
new mesh two API calls are required: in the first call (`mkernel_new_mesh`) a
new entry is created in the `meshInstances` vector, and in the second call
(`mkernel_set_state`) the entry is set using the nodes and edges information
obtained from the client. After this call, the mesh with all computed mappings 
is stored in memory and is ready to be used by the algorithms (the library keeps the mesh state).

The client now calls the `mkernel_refine_mesh_based_on_samples` function. In
the local scope of the `mkernel_refine_mesh_based_on_samples` function 
an instance of the `MeshRefinement` class is created, the Refine method is executed and the 
resulting mesh is saved in the `meshInstances` vector. 
The client retrieves the last state of the mesh using the `mkernel_get_mesh function`, where all information
required for rendering the new mesh (the nodes, the edges, and the faces) is copied 
from the library state to flat arrays. This extra copy is required for clients that
cannot communicate an array of structures at the API level (C# and C++ clients could potentially do this, FORTRAN clients not ).

By using this design only the mesh instances are saved throughout the
API calls and all other algorithm classes act as mesh modifiers, 
and are destroyed automatically after the API call is completed. 
Exceptions to this rule are the algorithms supporting interactivity.
In these cases, the algorithms are divided into several API methods, 
and their instances survive until an explicit "delete" method is invoked (e.g. `mkernel_orthogonalize_delete`).

![MeshKernel library simplified class diagram. Rectangles with right
corners represent structures, rectangles with round corners represent
classes and yellow rectangles represent namespaces.](figures/MeshKernelClassDiagram_1.jpg){#fig:classDiagram
width="100%"}

![Sequence diagram for creating a new grid and performing mesh
refinement.](figures/sequence_diagram_refinement.jpg){#fig:sequenceDiagram
width="100%"}

# The Operations file{#chap:operations}

In the current implementation, several geometrical methods acting on
arrays or simple types (e.g. Points) are collected in the
Operation.cpp file. This choice was made made because some geometric operations are often re-used in several classes.
Operations include:

-   Vector dot product.

-   Find the indices in a vector equal to a certain value.

-   Sort a vector returning its permutation array.

-   Finding the root of a function using the golden section search algorithm.

-   Performing coordinate transformations (e.g. from cartesian to spherical).

-   Inquire if a point is in a polygon.

-   Outer/inner products of two segments.

-   Compute the normal vector of a segment originating from a specific point.

-   Compute squared distances the distances between two points.

-   Compute the circumcenter of a triangle.

-   Compute if two lines are crossing and the propreties of the intersection.

-   Interpolating values using the averaging algorithm.

All operations reported above supports cartesian, spherical, and
spherical accurate coordinate systems.

# The Mesh class{#chap:mesh}

MeshKernel can handle 2d meshes and 1d meshes. Algorithms require cartain mappings to be available for both mesh1d and mesh2d, such as a mapping listing all edge indices connected to a particular node. 
The methods computing these mappings are shared between Mesh2D and Mesh1D, and implemented in the Mesh base class.
The mesh base class also contains other common data members, such as the node coordinate, the edges definitions, the face definitions and the mesh projection.
The mesh base class has the following responsibilities:

-   Construct the mesh faces from the nodes and edges and other mesh
    mappings required by all algorithms (Mesh::FindFaces).
    Mesh::FindFaces is using recursion to find faces with up to 6 edges (meshkernel::maximumNumberOfEdgesPerFace).

-   Supporting mesh editing, namely:

    -   Node merging

    -   Node insertion

    -   Moving a node

    -   Inserting edges

    -   Deleting edges

    -   Merging nodes (merging two nodes at meshkernel::mergingDistance). This algorithm use an r-tree for inquiring
adjacent nodes, see later.

-   Converting a curvilinear grid to an unstructured mesh (converting
    constructor).

-   Holding the mesh projection (cartesian, spherical, or spherical
    accurate).

-   Making a quad mesh from a polygon or from parameters.

-   Making a triangular mesh from a polygon. This algorithm introduces a dependency on 
the Richard Shewchuk Triangle.c library, added as an external componet in extern\triangle folder.

The public interface of the mesh class contains several algorithms modifying the mesh class members.
Most of them are trivial, and we refer to the doxygen api documentation.
Others are documented below:

###### RemoveSmallFlowEdges {#removesmallflowedges .unnumbered}

An unstructured mesh can be used to calculate water flow. This involves
a pressure gradient between the circumcenters of neighbouring faces.
That procedure is numerically unreliable when the distance between face
circumcenters (flow edges) becomes too small. Let's consider figure
[\[fig:coincide_circumcenters\]](#fig:coincide_circumcenters){reference-type="ref"
reference="fig:coincide_circumcenters"}:

::: {.center}
:::

The algorithm works as follow:

-   Any degenerated triangle (e.g. those having a coinciding node) is
    removed by collapsing the second and third node into the first one.

-   The edges crossing small flow edges are found. The flow edge length
    is computed from the face circumcenters and compared to an estimated
    cut off distance. The cutoff distance is computed using the face
    areas as follow:
    $$\textrm{cutOffDistance} = \textrm{threshold} \cdot 0.5 \cdot (\sqrt{\textrm{Area}_I}+\sqrt{\textrm{Area}_{II}})$$

-   All small flow edges are flagged with invalid indices and removed
    from the mesh. Removal occors in Administrate function.

###### RemoveSmallTrianglesAtBoundaries {#removesmalltrianglesatboundaries .unnumbered}

This algorithm removes triangles having the following properties:

-   The are at mesh boundary.

-   One or more neighboring faces are non-triangles.

-   The ratio of the face area to the average area of neighboring non
    triangles is less than a minimum ratio (defaults to 0.2).

-   The absolute cosine of one internal angle is less than 0.2.

These triangles having the above properties are merged by collapsing the
face nodes to the node having the minimum absolute cosine (e.g. the node
where the internal angle is closer to 90 degrees).

# The Mesh1D class{#chap:mesh1d}

Mesh1D is a base class derived from Mesh. 
A mesh 1d is composed of a series of connected edges representing 1d real word features, such as pipes or a sewage network.

# The Mesh2D class{#chap:mesh2d}

The mesh class represents an unstructured mesh. When communicating with
the client only unstructured meshes are used. Some algorithms generate
curvilinear grids (see section 9), but these are converted to a mesh
instance when communicating with the client. 

# The RTree class{#chap:rtree}

The mesh class stores two RTree class instances, used for inquiring the closest mesh nodes and edge to a point.
RTree is a class wrapping the boost::geometry::index::rtree code, adding an interface for performing common queries such as 
inquiring the nearest neighbors inside a specified distance(`meshkernel::RTree::NearestNeighborsOnSquaredDistance`) or a vector of the nearest neighbors (`meshkernel::RTree::NearestNeighbors`).
RTee has a `m_queryCache`, a vector used for collecting all query results and avoid frequent re-allocations when the number of results changes.

# The Contacts class{#chap:contacts}

The responsibility of the contacts class is connecting a 1d mesh to a 2d mesh. The class has a reference to the 1d and 2d mesh that will be connected. 
A connection is defined by the indices of the connected 1d node and 2d face. The following algorithms are available for 1d-2d connections:

-   `ComputeSingleConnections`: each non-boundary 1d node is connected to single 2d face. Figure 3 shows two 2d meshes, a 1d mesh between them, and the 1d-2d connections (in red).
The boundary nodes of the 1d mesh (those sharing only one 1d edge) are not connected to any 2d face. For the 1d nodes not overlapping a 2d mesh, a ray starting from the current node n is computed (dashed blue ray).
This ray is normal to the segment connecting the previous (n-1) and next one 1d node (n+1, the connecting segment is shown with a green dashed line). 
The ray is extended for 5 times the length of the connecting segment. 
The current 1d node is connected to the first boundary 2d face crossing the ray, first in the left direction and then in the right direction. 
By doing so a 1d mesh can be connected on the left and right sides of a mesh 2d boundary, for example when the 1d part represents a river and the 2d part the river banks.
The 1d nodes overlapping the 2d mesh are directly connected to the face including them.

![1d mesh connecting to 2d mesh using the ComputeSingleConnections algorithm. Connections are shown in red.](figures/ComputeSingleConnections.jpg){#fig:ComputeSingleConnections="100%"}

-   `ComputeMultipleConnections`: each internal 1d node is connected to multiple 2d faces. This type of connections should be used when the lengths of the 1d mesh edges are
considerably larger than the 2d mesh edges and generating a single connection for each 1d node is not representative.
In this algorithm, only the internal 1d nodes are connected. Figure 4 shows a 1d mesh overlapping a 2d mesh. 
For the node n, the closest 2d faces within a search radius are found and it is determined if those faces cross are crossed by the current 1d edge starting at node n and ending at and n+1. 
If the answer is positive, a connection is generated between the face and the closest 1d node composing the current 1d edge (i.e. n or n+1). The procedure is repeated for each 1d node.

![1d mesh connecting to 2d mesh using the ComputeSingleConnections algorithm. Connections are shown in red.](figures/ComputeMultipleConnections.jpg){#fig:ComputeMultipleConnections="100%"}

# The mesh OrthogonalizationAndSmoothing class{#chap:orthoandsmoothing}

This class implements the mesh orthogonalization and smoothing algorithm
as described in D-Flow FM technical manual (consult this manual for the
mathematical details on the equations). The algorithm is composed of two
differential equations: the first equation maximizes orthogonalization
between edges and flow links and the second equation reduces the
differences of the internal mesh angles (mesh smoothness). For this
reason, the OrthogonalizationAndSmoothing class is composed of a
smoother and an orthogonalizer, where the nodal contributions are
computed by separate classes, as opposed to the original Fortran
implementation. Essentially, the algorithm executes two loops:

-   An initialization step: The original mesh boundaries are saved. In
    case the mesh needs to be snapped to the land boundaries, the
    closest mesh edges to each land boundary are found
    (FindNearestMeshBoundary).

-   An outer loop, which itself is composed of the following steps:

    1.  Computation of the orthogonalizer contributions.

    2.  Computation of the smoother contributions.

    3.  Allocation of the linear system to be solved.

    4.  Summation of the two contributions (matrix assembly). The two
        contributions are weighted based on the desired smoothing to
        orthogonality ratio. OpenMP thread parallelization is used when
        summing the terms (loop iterations are independent).

-   An inner iteration: the resulting linear system is solved
    explicitly. The nodal coordinates are updated and the nodes moving
    on the mesh boundary are projected to the original mesh boundary
    (SnapMeshToOriginalMeshBoundary). In case a projection to land
    boundary is requested, the mesh nodes are projected to the land
    boundaries. An OpenMP parallelization is used in
    OrthogonalizationAndSmoothing::InnerIteration() because the update
    of the nodal coordinates was made iteration-independent.

MeshKernel API has five functions to enable the client to display the
mesh during the computations (interactivity). These functions are:

-   ggeo_orthogonalize_initialize

-   ggeo_orthogonalize_prepare_outer_iteration

-   ggeo_orthogonalize_inner_iteration.

-   ggeo_orthogonalize_finalize_outer_iteration

-   ggeo_orthogonalize_delete

The execution flow of these functions is shown in Figure A1 of the
Appendix. Additional details about these functions can be retrieved from
the API documentation.

# The MeshRefinement class{#chap:refinement}

Mesh refinement is based on iteratively splitting the edges until the
desired level of refinement or the maximum number of iterations is
reached. Refinement can be based on samples or based on a polygon. The
refinement based on samples uses the averaging interpolation algorithm
to compute the level of refinement from the samples to the centers of
the edges. At a high level, the mesh refinement is performed as follow:

-   Flag the nodes inside the refinement polygon.

-   Flag all face nodes of the faces not fully included in the polygon.

-   Execute the refinement iterations

    1.  For each edge store the index of its neighboring edge sharing a
        hanging node (the so-called brother edge). This is required for
        the following steps because edges with hanging nodes will not be
        divided further.

    2.  Compute edge and face refinement masks from the samples.

    3.  Compute if a face should be divided based on the computed
        refinement value.

    4.  Split the face by dividing the edges.

-   Connect the hanging nodes if required, thus forming triangular faces
    in the transition area.

As with OrthogonalizationAndSmoothing, MeshRefinement modifies an
existing mesh instance.

# The spline class{#chap:splines}

The spline class stores the corner points of each spline. Besides the
corner points, the derivatives at the corner points are also stored. The
coordinates of the points between the corner points are computed in the
static method Splines::Interpolate.

# The CurvilinearGridFromSplines class{#chap:curvifromsplines}

In this class, the algorithm to gradually develop a mesh from the
centreline of the channel towards the boundaries is implemented. It is
the most complex algorithm in the library. The curvilinear mesh is
developed from the center spline by the following steps:

-   Initialization step

    -   The splines are labelled (central or transversal spline) based
        on the number of corner points and the intersecting angles.

    -   The canal heights at a different position along the central
        spline are computed from the crossing splines.

    -   The normal vectors of each m part are computed, as these
        determine the growing front directions.

    -   The edge velocities to apply to each normal direction are
        computed.

-   Iteration step, where the mesh is grown of one layer at the time
    from the left and right sides of the central spline:

    -   Compute the node velocities from the edge velocities.

    -   Find the nodes at the front (the front might miss some faces and
        be irregular).

    -   Compute the maximum growth time to avoid faces with intersecting
        edges.

    -   Grow the grid by translating the nodes at the front by an amount
        equal to the product of the nodal velocity by the maximum grow
        time.

-   Post-processing

    -   Remove the skewed faces whose aspect rati0 exceeds a prescribed
        value.

    -   Compute the resulting CurvilinearGrid from the internal table of
        computed nodes (m_gridPoints).

to support interactivity with the client, the original Fortran algorithm
was divided into separate API calls:

-   ggeo_curvilinear_mesh_from_splines_ortho_initialize, corresponding
    to the initialization step above.

-   ggeo_curvilinear_mesh_from_splines_iteration, corresponding to the
    iteration step above.

-   ggeo_curvilinear_mesh_from_splines_ortho_refresh_mesh, corresponding
    to the post-processing above, plus the conversion of the
    CurvilinearGrid to an unstructured mesh.

-   ggeo_curvilinear_mesh_from_splines_ortho_delete, necessary to delete
    the CurvilinearGridFromSplines instance used in the previous API
    calls.

# The AveragingInterpolation class{#chap:averagingInterp}

The averaging interpolation operates on three specific locations: Faces
(m face mass centers), Nodes, and Edges(m edge centers). The idea is to
collect all samples close to the locations and perform a mathematical
operation on their values. The available operations are:

-   Simple averaging: computes a simple mean.

-   Closest: takes the value of the closest sample to the interpolation
    location.

-   Max: takes the maximum sample value.

-   Min: takes the minimum sample value.

-   InverseWeightedDistance: computes the inverse weighted sample mean.

-   MinAbsValue: computes the minimum absolute value.

The algorithm operates as follow:

-   The samples are ordered in an RTree for a fast search

-   The search area around the location is constructed as follow:

    1.  For face locations, the sample areas correspond to the faces,
        increased/decreased by the relativeSearchRadius parameter
        (relativeSearchRadius \> 1 increased, relativeSearchRadius \< 1
        decreased).

    2.  For Nodes and Edge locations, the dual face around the node is
        constructed by connecting the mid-points of all edges connected
        to the node. As above, the resulting polygon can be
        increased/decreased by the relativeSearchRadius parameter.

-   The search radius is computed from the constructed polygon nodes and
    the locations, the maximum value is taken.

-   The sample RTree is inquired to retrieve the indices of the samples
    within the search radius.

-   The operations described above are executed on the found samples.

-   For the Edges location, the interpolated values at the node are
    averaged.

# The TriangulationInterpolation class{#chap:triangInterp}

As for averaging, the triangle interpolation operates at three specific
locations: Faces, Nodes, and Edges. The idea is to triangulate the
samples and identify for each location the triangle that fully contains
it. Only the values at the nodes of the identified triangle are used in
the computation of each location. The algorithm operates as follow:

-   The triangulation of the samples is computed

-   For each triangle, the circumcentre is computed

-   The triangle circumcentres are ordered in an RTree for a fast search

-   For each location, the closest circumcentre is found

-   If the corresponding triangle contains the location then the linear
    interpolation is performed, otherwise, the next neighbouring
    triangle is searched. The next neighbouring triangle is the first
    triangle that satisfies these two conditions:

    1.  shares one of the current triangle edges.

    2.  the crossing of the edge and the line connecting the location
        the current triangle circumcentre exists.

-   If the next triangle does not contain the location, repeat the step
    above for a maximum number of times equal to two times the number of
    triangles.

When a triangle enclosing a specific location is not found, the
interpolated value at that location is invalid. The handling of
spherical accurate projection occurs at low-level geometrical functions
( IsPointInPolygonNodes, AreLineCrossing). Therefore, the algorithm is
independent of the implementation details that occur at the level of the
geometrical functions.

# Appendix{#chap:appendix}

![Sequence diagram for orthogonalization and
smoothing.](figures/sequence_diagram_orthogonalization.jpg){width="100%"}
