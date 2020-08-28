# Gridgeom mesh editing library

The Deltares library to generate and edit 2D unstructured and curvilinear meshes, suited for DFlowFM simulator.

The library is separated in an API namespace (GridGeomApi) used for communication with the client and a backend namespace (GridGeom), where the classes implementing the algorithms are included. 
The API namespace contains the library API methods (Gridgeom.cpp) and several structures used for communicating with the clients. These structures are mirrored in the client application and filled with appropriate values. 