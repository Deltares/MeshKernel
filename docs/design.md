
# Design {#Design}

The choices made for the representation of the geometric entities have
a strong influence on library design and performance. In MeshKernel an
unstructured mesh is uniquely defined by two entities:

-   The nodes vector: represented using `std::vector<Point>`, where
    `Point` is a structure containing the 2D coordinates of a point (cartesian or
    spherical).

-   The edges vector: represented using an `std::vector<std::pair<size_t,size_t>>`,
    containing the start and the end indices of the
    edges in the node vector described above.

All other mesh properties are computed from these two entities, such as
the face nodes, the face edges, the faces mass centers, and the faces
circumcenters. See \ref meshkernel::Mesh for some more details.

MeshKernel has an API namespace (\ref meshkernelapi)
and a back-end namespace (\ref meshkernel), where the classes implementing the algorithms 
are included. \ref meshkernelapi contains the library
API methods and several structures are used for communicating with the clients. 
These structures must be replicated in the client applications 
and filled with appropriate values.

When the client application creates a new mesh two API calls are required: in the first call (\ref meshkernelapi::mkernel_allocate_state) a
new entry is created in the `meshInstances` container, and in the second call
(e.g. \ref meshkernelapi::mkernel_mesh2d_set) the entry is set using the nodes and edges information
obtained from the client. After this call, the mesh with all computed mappings 
is stored in memory and is ready to be used by the algorithms (the library keeps the mesh state).

The client now calls the \ref meshkernelapi::mkernel_mesh2d_refine_based_on_samples function. In
the local scope of the \ref meshkernelapi::mkernel_mesh2d_refine_based_on_samples function 
an instance of the \ref meshkernel::MeshRefinement class is created, the Refine method is executed and the 
resulting mesh is saved in the `meshInstances` vector. 
The client retrieves the dimensions of the last state of the mesh using the \ref meshkernelapi::mkernel_get_dimensions_mesh2d function.
After allocating the necessary memory, the meshkernelapi::mkernel_get_mesh2d_data function should be called, where all information
required for rendering the new mesh (the nodes, the edges, and the faces) is copied 
from the library state to flat arrays. This extra copy is required for clients that
cannot communicate an array of structures at the API level (C# and C++ clients could potentially do this, FORTRAN clients not ).

By using this design only the mesh instances are saved throughout the
API calls and all other algorithm classes act as mesh modifiers, 
and are destroyed automatically after the API call is completed. 
Exceptions to this rule are the algorithms supporting interactivity.
In these cases, the algorithms are divided into several API methods, 
and their instances survive until an explicit "delete" method is invoked (e.g. \ref meshkernelapi::mkernel_mesh2d_delete_orthogonalization).
