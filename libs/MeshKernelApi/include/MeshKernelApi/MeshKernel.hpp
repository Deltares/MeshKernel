//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2021.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// contact: delft3d.support@deltares.nl
// Stichting Deltares
// P.O. Box 177
// 2600 MH Delft, The Netherlands
//
// All indications and logos of, and references to, "Delft3D" and "Deltares"
// are registered trademarks of Stichting Deltares, and remain the property of
// Stichting Deltares. All rights reserved.
//
//------------------------------------------------------------------------------

#pragma once

#include <MeshKernel/Parameters.hpp>

#include <MeshKernelApi/Contacts.hpp>
#include <MeshKernelApi/CurvilinearGrid.hpp>
#include <MeshKernelApi/GeometryList.hpp>
#include <MeshKernelApi/GriddedSamples.hpp>
#include <MeshKernelApi/Mesh1D.hpp>
#include <MeshKernelApi/Mesh2D.hpp>

#if defined(_WIN32)
#if !defined(MKERNEL_API)
#define MKERNEL_API __declspec(dllexport)
#endif
#else
#define MKERNEL_API __attribute__((visibility("default")))
#endif

/// \namespace meshkernelapi
/// @brief Contains all structs and functions exposed at the API level
namespace meshkernelapi
{
    /// @brief Enumeration for api error types
    enum MeshKernelApiErrors
    {
        Success = 0,
        NotImplemented = 1,
        MeshGeometryError = 2,
        AlgorithmError = 3,
        MeshKernelError = 4,
        StadardLibraryException = 5,
        UnknownException = 6
    };

#ifdef __cplusplus
    extern "C"
    {
#endif
        /// @brief Creates a new mesh state and returns the generated \p meshKernelId
        /// @param[in] projectionType  Cartesian (0), spherical (1) or spherical accurate(2) state
        /// @param[out] meshKernelId The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_allocate_state(int projectionType, int& meshKernelId);

        /// @brief Deallocate mesh state
        /// @param[in] meshKernelId The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_deallocate_state(int meshKernelId);

        /// @brief Deletes a mesh in a polygon using several options
        /// @param[in] meshKernelId   The id of the mesh state
        /// @param[in] polygon        The polygon where to perform the operation
        /// @param[in] deletionOption The deletion option \ref meshkernel::Mesh2D::DeleteMeshOptions
        /// @param[in] invertDeletion Whether to invert the deletion of selected features (0 no, 1 yes)
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_delete(int meshKernelId,
                                              const GeometryList& polygon,
                                              int deletionOption,
                                              int invertDeletion);

        /// @brief Sets the meshkernel::Mesh2D state
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] mesh2d       The Mesh2D data
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_set(int meshKernelId, const Mesh2D& mesh2d);

        /// @brief Sets the meshkernel::Mesh1D state
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] mesh1d       The Mesh1D data
        /// @returns Error code
        MKERNEL_API int mkernel_mesh1d_set(int meshKernelId, const Mesh1D& mesh1d);

        /// @brief Sets the meshkernel::Network1D state
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] polylines    The polylines describing the network
        /// @returns Error code
        MKERNEL_API int mkernel_network1d_set(int meshKernelId, const GeometryList& polylines);

        /// @brief Compute the network chainages from fixed point locations
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] fixedChainages The fixed chainages for each polyline. Chunks are separated by the separator, each chunk corresponds to a polyline
        /// @param[in] sizeFixedChainages The size of fixed chainages vector
        /// @param[in] minFaceSize The minimum face size. The distance between two chainages must be no less than this length
        /// @param[in] fixedChainagesOffset The offset to use for fixed chainages
        /// @return Error code
        MKERNEL_API int mkernel_network1d_compute_fixed_chainages(int meshKernelId,
                                                                  double* fixedChainages,
                                                                  int sizeFixedChainages,
                                                                  double minFaceSize,
                                                                  double fixedChainagesOffset);

        /// @brief Compute the network chainages at a regular offset
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] offset The offset between points
        /// @return Error code
        MKERNEL_API int mkernel_network1d_compute_offsetted_chainages(int meshKernelId, double offset);

        /// @brief Convert network chainages to mesh1d nodes and edges
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] minFaceSize The minimum face size below which two nodes will be merged
        /// @return Error code
        MKERNEL_API int mkernel_network1d_to_mesh1d(int meshKernelId, double minFaceSize);

        /// @brief Gets the Mesh2D dimensions
        ///
        /// The integer parameters of the Mesh2D struct are set to the corresponding dimensions
        /// The pointers are set to null, and must be set to correctly sized memory
        /// before passing the struct to `mkernel_mesh2d_get_data`.
        /// @param[in]  meshKernelId The id of the mesh state
        /// @param[out] mesh2d       The Mesh2D dimensions
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_get_dimensions(int meshKernelId,
                                                      Mesh2D& mesh2d);

        /// @brief Gets the Mesh2D dimensions data
        ///
        /// This function ought to be called after `mkernel_mesh2d_get_dimensions` has been called
        /// and the pointers have been set to correctly sized memory.
        /// @param[in]     meshKernelId The id of the mesh state
        /// @param[in,out] mesh2d       The Mesh2D data
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_get_data(int meshKernelId, Mesh2D& mesh2d);

        /// @brief Gets the curvilinear grid dimensions as a CurvilinearGrid struct (converted as set of edges and nodes).
        ///
        /// The integer parameters of the CurvilinearGrid struct are set to the corresponding dimensions
        /// The pointers are set to null, and must be set to correctly sized memory
        /// before passing the struct to `mkernel_curvilinear_get_data`.
        /// @param[in]  meshKernelId    The id of the mesh state.
        /// @param[out] curvilinearGrid The structure containing the dimensions of the curvilinear grid.
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_get_dimensions(int meshKernelId, CurvilinearGrid& curvilinearGrid);

        /// @brief Gets the curvilinear grid data as a CurvilinearGrid struct (converted as set of edges and nodes)
        ///
        /// This function ought to be called after `mkernel_get_curvilinear_dimension` has been called
        /// and the pointers have been set to correctly sized memory.
        /// @param[in]  meshKernelId    The id of the mesh state
        /// @param[out] curvilinearGrid The structure containing the curvilinear grid arrays.
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_get_data(int meshKernelId, CurvilinearGrid& curvilinearGrid);

        /// @brief Gets the Mesh1D dimensions
        ///
        /// The integer parameters of the Mesh1D struct are set to the corresponding dimensions
        /// The pointers are set to null, and must be set to correctly sized memory
        /// before passing the struct to `mkernel_mesh1d_get_dimensions`
        /// @param[in]  meshKernelId The id of the mesh state
        /// @param[out] mesh1d       The structure containing the dimensions of the Mesh1D.
        /// @returns Error code
        MKERNEL_API int mkernel_mesh1d_get_dimensions(int meshKernelId, Mesh1D& mesh1d);

        /// @brief Gets the Mesh1D dimensions data
        ///
        /// This function ought to be called after `mkernel_mesh1d_get_dimensions` has been called
        /// and the pointers have been set to correctly sized memory
        /// @param[in]     meshKernelId The id of the mesh state
        /// @param[in,out] mesh1d       The structure containing the Mesh1D arrays.
        /// @returns Error code
        MKERNEL_API int mkernel_mesh1d_get_data(int meshKernelId, Mesh1D& mesh1d);

        /// @brief Gets the number of 1d-2d contacts
        /// @param[in]  meshKernelId           The id of the mesh state
        /// @param[out] contacts               Contacts data
        /// @returns                           Error code
        MKERNEL_API int mkernel_contacts_get_dimensions(int meshKernelId, Contacts& contacts);

        /// @brief Gets the 1d-2d contacts indices (from index / to indices)
        /// @param[in]  meshKernelId           The id of the mesh state
        /// @param[out] contacts               Contacts data
        /// @returns                           Error code
        MKERNEL_API int mkernel_contacts_get_data(int meshKernelId, Contacts& contacts);

        /// @brief Count the number of hanging edges in a mesh2d.
        /// An hanging edge is an edge where one of the two nodes is not connected.
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[out] numEdges The number of hanging edges
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_count_hanging_edges(int meshKernelId, int& numEdges);

        /// @brief Gets the indices of hanging edges. An hanging edge is an edge where one of the two nodes is not connected.
        /// @param[in]     meshKernelId        The id of the mesh state
        /// @param[in,out] edges Pointer to memory where the indices of the hanging edges will be stored
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_get_hanging_edges(int meshKernelId, int* edges);

        /// @brief Deletes all hanging edges. An hanging edge is an edge where one of the two nodes is not connected.
        /// @param[in] meshKernelId The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_delete_hanging_edges(int meshKernelId);

        /// @brief Mesh2d orthogonalization.
        ///
        /// The function modifies the mesh for achieving orthogonality between the edges and the segments connecting the face circumcenters.
        /// The amount of orthogonality is traded against the mesh smoothing (in this case the equality of face areas).
        /// The parameter to regulate the amount of orthogonalization is contained in  \ref meshkernel::OrthogonalizationParameters::orthogonalization_to_smoothing_factor
        /// @param[in] meshKernelId                The id of the mesh state
        /// @param[in] projectToLandBoundaryOption The option to determine how to snap to land boundaries
        /// @param[in] orthogonalizationParameters The structure containing the orthogonalization parameters \ref meshkernel::OrthogonalizationParameters
        /// @param[in] selectingPolygon            The polygon where to perform the orthogonalization  (num_coordinates = 0 for an empty polygon)
        /// @param[in] landBoundaries              The land boundaries to account for in the orthogonalization process (num_coordinates = 0 for no land boundaries)
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_compute_orthogonalization(int meshKernelId,
                                                                 int projectToLandBoundaryOption,
                                                                 const meshkernel::OrthogonalizationParameters& orthogonalizationParameters,
                                                                 const GeometryList& selectingPolygon,
                                                                 const GeometryList& landBoundaries);

        /// @brief Initialization of the \ref meshkernel::OrthogonalizationAndSmoothing algorithm.
        ///
        /// This is the first function to call when using orthogonalization in interactive mode (visualizing the grid while it is orthogonalizing),
        /// in order to set the internal state of the algorithm reused during the iterations.
        /// @param[in] meshKernelId                The id of the mesh state
        /// @param[in] projectToLandBoundaryOption The option to determine how to snap to land boundaries
        /// @param[in] orthogonalizationParameters The structure containing the user defined orthogonalization parameters
        /// @param[in] selectingPolygon            The polygon where to perform the orthogonalization (num_coordinates = 0 for an empty polygon)
        /// @param[in] landBoundaries              The land boundaries to account for in the orthogonalization process  (num_coordinates = 0 for no land boundaries)
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_initialize_orthogonalization(int meshKernelId,
                                                                    int projectToLandBoundaryOption,
                                                                    meshkernel::OrthogonalizationParameters& orthogonalizationParameters,
                                                                    const GeometryList& selectingPolygon,
                                                                    const GeometryList& landBoundaries);

        /// @brief Prepares an outer orthogonalization iteration, computing the new orthogonalization and smoothing weights from the modified mesh geometry (in interactive mode).
        ///
        /// `mkernel_mesh2d_initialize_orthogonalization` function must be called before.
        /// @param[in] meshKernelId The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_prepare_outer_iteration_orthogonalization(int meshKernelId);

        /// @brief Performs inner orthogonalization iteration, by slowly moving the mesh nodes to new optimal positions (interactive mode).
        ///
        /// `mkernel_mesh2d_prepare_outer_iteration_orthogonalization` function must be called before.
        /// @param[in] meshKernelId The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_compute_inner_ortogonalization_iteration(int meshKernelId);

        /// @brief Finalizes the orthogonalization outer iteration, computing the new coefficients for grid adaption and the new face circumcenters (interactive mode).
        /// @param[in] meshKernelId  The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_finalize_inner_ortogonalization_iteration(int meshKernelId);

        /// @brief Cleans the orthogonalization algorithm state, allocated in `mkernel_mesh2d_initialize_orthogonalization` (interactive mode)
        /// @param[in] meshKernelId The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_delete_orthogonalization(int meshKernelId);

        /// @brief Gets the mesh orthogonality, expressed as the ratio between the edges and the segments connecting the face circumcenters.
        /// @param[in]  meshKernelId The id of the mesh state
        /// @param[out] geometryList The orthogonality values of each edge
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_get_orthogonality(int meshKernelId, GeometryList& geometryList);

        /// @brief Gets the smoothness, expressed as the ratio between the values of two neighboring faces areas.
        /// @param[in]  meshKernelId The id of the mesh state
        /// @param[out] geometryList The smoothness values at each edge
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_get_smoothness(int meshKernelId, GeometryList& geometryList);

        /// @brief Get the computed spline points between two corner nodes
        /// @param[in]  geometryListIn  The input corner nodes of the splines
        /// @param[out] geometryListOut The output spline
        /// @param[out] numberOfPointsBetweenNodes The number of spline points to generate between two corner nodes.
        /// @returns Error code
        MKERNEL_API int mkernel_get_splines(const GeometryList& geometryListIn,
                                            GeometryList& geometryListOut,
                                            int numberOfPointsBetweenNodes);

        /// @brief Generates a triangular mesh2d grid within a polygon. The size of the triangles is determined from the length of the polygon edges.
        /// @param[in] meshKernelId  The id of the mesh state
        /// @param[in] polygonPoints The polygon where to triangulate
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_make_mesh_from_polygon(int meshKernelId, const GeometryList& polygonPoints);

        /// @brief Makes a triangular mesh from  a set of samples, triangulating the sample points.
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] samples The samples where to triangulate
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_make_mesh_from_samples(int meshKernelId, const GeometryList& samples);

        /// @brief Makes uniform meshes from a series of polygons
        /// @param[in] meshKernelId       The id of the mesh state
        /// @param[in] makeGridParameters The structure containing the make grid parameters
        /// @param[in] geometryList       The polygons to account for
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_make_uniform(int meshKernelId,
                                                    const meshkernel::MakeGridParameters& makeGridParameters,
                                                    const GeometryList& geometryList);

        /// @brief Makes uniform mesh on a defined extension
        /// @param[in] meshKernelId       The id of the mesh state
        /// @param[in] makeGridParameters The structure containing the make grid parameters
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_make_uniform_on_extension(int meshKernelId,
                                                                 const meshkernel::MakeGridParameters& makeGridParameters);

        /// @brief Retrieves the boundaries of a mesh as a series of separated polygons.
        ///
        /// For example, if a mesh has an single inner hole, two polygons will be generated, one for the inner boundary and one for the outer boundary.
        /// @param[in]  meshKernelId The id of the mesh state
        /// @param[out] boundaryPolygons The output network boundary polygon
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_get_mesh_boundaries_as_polygons(int meshKernelId, GeometryList& boundaryPolygons);

        /// @brief Counts the number of polygon nodes contained in the mesh boundary polygons computed in function `mkernel_mesh2d_get_mesh_boundaries_as_polygons`
        /// @param[in]  meshKernelId         The id of the mesh state
        /// @param[out] numberOfPolygonNodes The number of polygon nodes
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_count_mesh_boundaries_as_polygons(int meshKernelId, int& numberOfPolygonNodes);

        /// @brief Gets the edges intersected by a polyline, with additional information on the intersections
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] boundaryPolyLine An input polyline, defined as a series of points
        /// @param[out] edgeNodes The indices of the intersected edge nodes. The first node of the edge is on the left (the virtual node), the second node of the edge is on the right (the inner node)
        /// @param[out] edgeIndex For each intersected edge, the edge index
        /// @param[out] edgeDistances For each intersection, the location of the intersection expressed as adimensional distance from the edge starting node
        /// @param[out] segmentDistances For each intersection, the location of the intersection expressed as adimensional distance from the polyline segment start
        /// @param[out] segmentIndexes For each intersection, the segment index
        /// @param[out] faceIndexes For each intersection, the face index
        /// @param[out] faceNumEdges For each intersection, the number of intersections
        /// @param[out] faceEdgeIndex For each intersection, the index of the intersected edge
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_intersections_from_polyline(int meshKernelId,
                                                                   const GeometryList& boundaryPolyLine,
                                                                   int* edgeNodes,
                                                                   int* edgeIndex,
                                                                   double* edgeDistances,
                                                                   double* segmentDistances,
                                                                   int* segmentIndexes,
                                                                   int* faceIndexes,
                                                                   int* faceNumEdges,
                                                                   int* faceEdgeIndex);

        /// @brief Refines the polygon perimeter between two nodes. This interval is refined to achieve a target edge length.
        ///
        /// The function is often used before `mkernel_mesh2d_make_mesh_from_polygon`, for generating a triangular mesh where edges have a desired length.
        /// @param[in]  meshKernelId       The id of the mesh state
        /// @param[in]  polygonToRefine    The input polygon to refine
        /// @param[in]  firstNodeIndex     The first index of the refinement interval
        /// @param[in]  secondNodeIndex    The second index of the refinement interval
        /// @param[in]  targetEdgeLength   The target interval edge length
        /// @param[out] refinedPolygon     The refined polygon
        /// @returns Error code
        MKERNEL_API int mkernel_polygon_refine(int meshKernelId,
                                               const GeometryList& polygonToRefine,
                                               int firstNodeIndex,
                                               int secondNodeIndex,
                                               double targetEdgeLength,
                                               GeometryList& refinedPolygon);

        /// @brief Counts the number of polygon nodes resulting from polygon refinement with `mkernel_polygon_refine`.
        ///
        /// This function should be used by clients before `mkernel_polygon_refine` for allocating \ref GeometryList containing the refinement result.
        /// @param[in] meshKernelId          The id of the mesh state
        /// @param[in] polygonToRefine       The input polygon to refine
        /// @param[in] firstIndex            The first index of the refinement interval
        /// @param[in] secondIndex           The second index of the refinement interval
        /// @param[in] distance              The target interval edge length
        /// @param[out] numberOfPolygonNodes The number of nodes after refinement
        /// @returns Error code
        MKERNEL_API int mkernel_polygon_count_refine(int meshKernelId,
                                                     const GeometryList& polygonToRefine,
                                                     int firstIndex,
                                                     int secondIndex,
                                                     double distance,
                                                     int& numberOfPolygonNodes);

        /// @brief Merges the mesh2d nodes within a distance of 0.001 m, effectively removing all small edges
        /// @param[in] meshKernelId   The id of the mesh state
        /// @param[in] geometryListIn The polygon defining the area where the operation will be performed
        /// @param[in] mergingDistance The distance below which two nodes will be merged
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_merge_nodes(int meshKernelId, const GeometryList& geometryListIn, double mergingDistance);

        /// @brief Merges two mesh2d nodes into one.
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] firstNode    The index of the first node to merge
        /// @param[in] secondNode   The index of the second node to merge
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_merge_two_nodes(int meshKernelId, int firstNode, int secondNode);

        /// @brief Gets the indices of the mesh2d nodes selected with a polygon.
        /// @param[in]  meshKernelId   The id of the mesh state
        /// @param[in]  geometryListIn The input polygon
        /// @param[in]  inside         Selection of the nodes inside the polygon (1) or outside (0)
        /// @param[out] selectedNodes  The selected nodes indices
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_get_nodes_in_polygons(int meshKernelId,
                                                             const GeometryList& geometryListIn,
                                                             int inside,
                                                             int* selectedNodes);

        /// @brief Counts the number of selected mesh node indices.
        ///
        /// This function should be used by clients before `mkernel_mesh2d_get_nodes_in_polygons` for allocating an integer array storing the selection results.
        /// @param[in]  meshKernelId      The id of the mesh state
        /// @param[in]  geometryListIn    The input polygon
        /// @param[in]  inside            Selection of the nodes inside the polygon (1) or outside (0)
        /// @param[out] numberOfMeshNodes The number of selected nodes
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_count_nodes_in_polygons(int meshKernelId,
                                                               const GeometryList& geometryListIn,
                                                               int inside,
                                                               int& numberOfMeshNodes);

        /// @brief Insert a new mesh2d edge connecting two nodes
        /// @param[in]  meshKernelId The id of the mesh state
        /// @param[in]  startNode    The index of the first node to connect
        /// @param[in]  endNode      The index of the second node to connect
        /// @param[out] newEdgeIndex The index of the new generated edge
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_insert_edge(int meshKernelId, int startNode, int endNode, int& newEdgeIndex);

        /// @brief Insert a new mesh2d node at a specific coordinate.
        /// @param[in]  meshKernelId The id of the mesh state
        /// @param[in]  xCoordinate  The x coordinate of the node to insert
        /// @param[in]  yCoordinate  The y coordinate of the node to insert
        /// @param[out] nodeIndex    The index of the new mesh node
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_insert_node(int meshKernelId, double xCoordinate, double yCoordinate, int& nodeIndex);

        /// @brief Deletes a mesh2d node
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] nodeIndex    The index of the node to delete
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_delete_node(int meshKernelId, int nodeIndex);

        /// @brief Moves a mesh2d node to a new position
        /// @param[in] meshKernelId    The id of the mesh state
        /// @param[in] xCoordinate     The new x coordinate of the node to move
        /// @param[in] yCoordinate     The new y coordinate of the node to move
        /// @param[in] nodeIndex       The index of the mesh2d node to be moved
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_move_node(int meshKernelId, double xCoordinate, double yCoordinate, int nodeIndex);

        /// @brief Deletes the closest mesh2d edge to a point.
        /// The coordinates of the edge middle points are used for calculating the distances to the point.
        /// @param[in] meshKernelId   The id of the mesh state
        /// @param[in] xCoordinate    The x coordinate of the point
        /// @param[in] yCoordinate    The y coordinate of the point
        /// @param[in]  xLowerLeftBoundingBox  The x coordinate of the lower left corner of the bounding box
        /// @param[in]  yLowerLeftBoundingBox  The y coordinate of the lower left corner of the bounding box
        /// @param[in]  xUpperRightBoundingBox The x coordinate of the upper right corner of the bounding box
        /// @param[in]  yUpperRightBoundingBox The y coordinate of the upper right corner of the bounding box
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_delete_edge(int meshKernelId,
                                                   double xCoordinate,
                                                   double yCoordinate,
                                                   double xLowerLeftBoundingBox,
                                                   double yLowerLeftBoundingBox,
                                                   double xUpperRightBoundingBox,
                                                   double yUpperRightBoundingBox);

        /// @brief Gets the closest mesh2d edge to a point in a bounding box
        /// @param[in]  meshKernelId  The id of the mesh state
        /// @param[in]  xCoordinate   The x coordinate of the point
        /// @param[in]  yCoordinate   The y coordinate of the point
        /// @param[in]  xLowerLeftBoundingBox  The x coordinate of the lower left corner of the bounding box
        /// @param[in]  yLowerLeftBoundingBox  The y coordinate of the lower left corner of the bounding box
        /// @param[in]  xUpperRightBoundingBox The x coordinate of the upper right corner of the bounding box
        /// @param[in]  yUpperRightBoundingBox The y coordinate of the upper right corner of the bounding box
        /// @param[out] edgeIndex   The found edge index
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_get_edge(int meshKernelId,
                                                double xCoordinate,
                                                double yCoordinate,
                                                double xLowerLeftBoundingBox,
                                                double yLowerLeftBoundingBox,
                                                double xUpperRightBoundingBox,
                                                double yUpperRightBoundingBox,
                                                int& edgeIndex);

        /// @brief Generate a new polygon from an existing one by offsetting the perimeter by a given distance.
        ///
        /// Offsetting can be done inward or outward the existing polygon.
        /// @param[in] meshKernelId     The id of the mesh state
        /// @param[in] geometryListIn   The polygon to offset
        /// @param[in] inWard           Compute the inner offset (1) or outer offset offset (0)
        /// @param[in] distance         The offset distance
        /// @param[out] geometryListOut The resulting offset polygon
        /// @returns Error code
        MKERNEL_API int mkernel_polygon_get_offset(int meshKernelId,
                                                   const GeometryList& geometryListIn,
                                                   int inWard,
                                                   double distance,
                                                   GeometryList& geometryListOut);

        /// @brief Counts the number of polygon nodes resulting from polygon offset.
        ///
        /// This function should be used by clients before `mkernel_polygon_get_offset` for allocating the \ref GeometryList containing the offset result.
        /// @param[in] meshKernelId          The id of the mesh state
        /// @param[in] geometryListIn        The polygon to offset
        /// @param[in] innerPolygon          Whether to compute inner (1) or outer offset (0)
        /// @param[in] distance              The offset distance
        /// @param[out] numberOfPolygonNodes The number of nodes in the offset polygon
        /// @returns Error code
        MKERNEL_API int mkernel_polygon_count_offset(int meshKernelId,
                                                     const GeometryList& geometryListIn,
                                                     int innerPolygon,
                                                     double distance,
                                                     int& numberOfPolygonNodes);

        /// @brief Refines a mesh2d based on samples. Refinement is achieved by successive splits of the face edges.
        ///
        /// The number of successive splits is indicated on the sample value.
        /// For example a value of 0 means no split and hence no refinement, a value of 1 a single split (a quadrilateral face generates 4 faces),
        /// a value of 2 two splits (a quadrilateral face generates 16 faces).
        /// @param[in] meshKernelId             The id of the mesh state
        /// @param[in] samples                  The sample set
        /// @param[in] relativeSearchRadius     The relative search radius relative to the face size, used for some interpolation algorithms
        /// @param[in] minimumNumSamples        The minimum number of samples used for some averaging algorithms
        /// @param[in] meshRefinementParameters The mesh refinement parameters
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_refine_based_on_samples(int meshKernelId,
                                                               const GeometryList& samples,
                                                               double relativeSearchRadius,
                                                               int minimumNumSamples,
                                                               const meshkernel::MeshRefinementParameters& meshRefinementParameters);

        /// @brief Refine based on gridded samples
        ///
        /// The number of successive splits is indicated on the sample value.
        /// For example a value of 0 means no split and hence no refinement, a value of 1 a single split (a quadrilateral face generates 4 faces),
        /// a value of 2 two splits (a quadrilateral face generates 16 faces).
        /// @param[in] meshKernelId             The id of the mesh state
        /// @param[in] griddedSamples                  The gridded samples
        /// @param[in] meshRefinementParameters The mesh refinement parameters
        /// @param[in] useNodalRefinement       Use nodal refinement
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_refine_based_on_gridded_samples(int meshKernelId,
                                                                       const GriddedSamples& griddedSamples,
                                                                       const meshkernel::MeshRefinementParameters& meshRefinementParameters,
                                                                       bool useNodalRefinement);

        /// @brief Refines a mesh2d within a polygon. Refinement is achieved by splitting the edges contained in the polygon by two.
        /// @param[in] meshKernelId             The id of the mesh state
        /// @param[in] geometryList             The closed polygon where to perform the refinement
        /// @param[in] meshRefinementParameters The mesh refinement parameters
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_refine_based_on_polygon(int meshKernelId,
                                                               const GeometryList& geometryList,
                                                               const meshkernel::MeshRefinementParameters& meshRefinementParameters);

        /// @brief Finds the mesh2d node closest to a point, within a search radius.
        /// @param[in]  meshKernelId  The id of the mesh state
        /// @param[in]  xCoordinate   The x coordinate of the point
        /// @param[in]  yCoordinate   The y coordinate of the point
        /// @param[in]  searchRadius  The search radius
        /// @param[in]  xLowerLeftBoundingBox  The x coordinate of the lower left corner of the bounding box
        /// @param[in]  yLowerLeftBoundingBox  The y coordinate of the lower left corner of the bounding box
        /// @param[in]  xUpperRightBoundingBox The x coordinate of the upper right corner of the bounding box
        /// @param[in]  yUpperRightBoundingBox The y coordinate of the upper right corner of the bounding box
        /// @param[out] nodeIndex     The index of the found node
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_get_node_index(int meshKernelId,
                                                      double xCoordinate,
                                                      double yCoordinate,
                                                      double searchRadius,
                                                      double xLowerLeftBoundingBox,
                                                      double yLowerLeftBoundingBox,
                                                      double xUpperRightBoundingBox,
                                                      double yUpperRightBoundingBox,
                                                      int& nodeIndex);

        /// @brief Gets the closest mesh2d node coordinates to a point, searching within a radius.
        /// @param[in]  meshKernelId    Id of the grid state
        /// @param[in]  xCoordinateIn   The x coordinate of the node to insert
        /// @param[in]  yCoordinateIn   The y coordinate of the node to insert
        /// @param[in]  searchRadius    The radii where to search for mesh nodes
        /// @param[in]  xLowerLeftBoundingBox  The x coordinate of the lower left corner of the bounding box
        /// @param[in]  yLowerLeftBoundingBox  The y coordinate of the lower left corner of the bounding box
        /// @param[in]  xUpperRightBoundingBox The x coordinate of the upper right corner of the bounding box
        /// @param[in]  yUpperRightBoundingBox The y coordinate of the upper right corner of the bounding box
        /// @param[out] xCoordinateOut  The x coordinate of the found Mesh2D node
        /// @param[out] yCoordinateOut  The y coordinate of the found Mesh2D node
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_get_closest_node(int meshKernelId,
                                                        double xCoordinateIn,
                                                        double yCoordinateIn,
                                                        double searchRadius,
                                                        double xLowerLeftBoundingBox,
                                                        double yLowerLeftBoundingBox,
                                                        double xUpperRightBoundingBox,
                                                        double yUpperRightBoundingBox,
                                                        double& xCoordinateOut,
                                                        double& yCoordinateOut);

        /// @brief Selects the polygon nodes within another polygon.
        /// @param[in]  meshKernelId   The id of the mesh state
        /// @param[in]  selectingPolygon   The selecting polygon (num_coordinates = 0 for an empty polygon)
        /// @param[in]  polygonToSelect    The polygon to select
        /// @param[out] selectionResults   The selection result, contained in the in the values field of \ref GeometryList (0.0 not selected, 1.0 selected).
        /// Note that the selection selectionResults variable must be allocated by the client.
        /// @returns Error code
        MKERNEL_API int mkernel_polygon_get_included_points(int meshKernelId,
                                                            const GeometryList& selectingPolygon,
                                                            const GeometryList& polygonToSelect,
                                                            GeometryList& selectionResults);

        /// @brief Flips mesh2d edges, to optimize the mesh smoothness. This operation is usually performed after `mkernel_mesh2d_refine_based_on_samples` or `mkernel_mesh2d_refine_based_on_polygon`.
        ///
        /// Nodes that are connected to more than six other nodes are typically enclosed by faces of highly non-uniform shape and wildly varying areas.
        /// @param[in] meshKernelId                  The id of the mesh state
        /// @param[in] isTriangulationRequired       The option to triangulate also non triangular cells (if activated squares becomes triangles)
        /// @param[in] projectToLandBoundaryRequired The option to determine how to snap to land boundaries
        /// @param[in] selectingPolygon              The polygon where to perform the edge flipping (num_coordinates = 0 for an empty polygon)
        /// @param[in] landBoundaries                The land boundaries to account for when flipping the edges (num_coordinates = 0 for no land boundaries)
        /// @returns Error code
        MKERNEL_API int mkernel_mesh2d_flip_edges(int meshKernelId,
                                                  int isTriangulationRequired,
                                                  int projectToLandBoundaryRequired,
                                                  const GeometryList& selectingPolygon,
                                                  const GeometryList& landBoundaries);

        /// @brief Gets the number of obtuse mesh2d triangles. Obtuse triangles are those having one edge longer than the sum of the other two.
        /// @param[in]  meshKernelId       The id of the mesh state
        /// @param[out] numObtuseTriangles The number of obtuse triangles
        /// @return Error code
        MKERNEL_API int mkernel_mesh2d_count_obtuse_triangles(int meshKernelId, int& numObtuseTriangles);

        /// @brief Gets the mass centers of obtuse mesh2d triangles. Obtuse triangles are those having one edge longer than the sum of the other two.
        /// @param[in]  meshKernelId  The id of the mesh state
        /// @param[out] result        The coordinates of the obtuse triangles mass centers stored in coordinates_x and coordinates_y of a \ref GeometryList
        /// @return Error code
        MKERNEL_API int mkernel_mesh2d_get_obtuse_triangles_mass_centers(int meshKernelId, GeometryList& result);

        /// @brief Counts the number of small mesh2d flow edges. The flow edges are the edges connecting faces circumcenters.
        /// @param[in] meshKernelId                  The id of the mesh state
        /// @param[in] smallFlowEdgesLengthThreshold The configurable length for detecting a small flow edge
        /// @param[out] numSmallFlowEdges            The number of the small flow edges
        /// @return Error code
        MKERNEL_API int mkernel_mesh2d_count_small_flow_edge_centers(int meshKernelId,
                                                                     double smallFlowEdgesLengthThreshold,
                                                                     int& numSmallFlowEdges);

        /// @brief Gets the small mesh2d flow edges. The flow edges are the edges connecting faces circumcenters.
        /// @param[in]  meshKernelId            The id of the mesh state
        /// @param[in]  smallFlowEdgesThreshold The configurable threshold for detecting a small flow edge
        /// @param[out] result                  The middle points of the small flow edges, stored in coordinates_x and coordinates_y of a \ref GeometryList
        /// @return Error code
        MKERNEL_API int mkernel_mesh2d_get_small_flow_edge_centers(int meshKernelId,
                                                                   double smallFlowEdgesThreshold,
                                                                   GeometryList& result);

        /// @brief Deletes all small mesh2d flow edges and small triangles. The flow edges are the edges connecting faces circumcenters.
        /// @param[in] meshKernelId               The id of the mesh state
        /// @param[in] smallFlowEdgesThreshold    The configurable threshold for detecting the small flow edges
        /// @param[in] minFractionalAreaTriangles The ratio of the face area to the average area of neighboring non triangular faces.
        ///                                       This parameter is used for determining if a triangular face is small.
        /// @return Error code
        MKERNEL_API int mkernel_mesh2d_delete_small_flow_edges_and_small_triangles(int meshKernelId,
                                                                                   double smallFlowEdgesThreshold,
                                                                                   double minFractionalAreaTriangles);

        /// @brief Computes 1d-2d contacts, where each single 1d node is connected to one mesh2d face circumcenter (ggeo_make1D2Dinternalnetlinks_dll)
        ///
        /// \see meshkernel::Contacts::ComputeSingleContacts
        /// @param[in]  meshKernelId  The id of the mesh state
        /// @param[in]  oneDNodeMask  The mask to apply to 1d nodes (1 = connect node, 0 = do not connect)
        /// @param[in]  polygons      The polygons selecting the area where the 1d-2d contacts will be generated.
        /// @param[in] projectionFactor     The projection factor used for generating links when 1d nodes are not inside mesh2d
        /// @return Error code
        MKERNEL_API int mkernel_contacts_compute_single(int meshKernelId,
                                                        const int* oneDNodeMask,
                                                        const GeometryList& polygons,
                                                        double projectionFactor);

        /// @brief Computes 1d-2d contacts, where a single 1d node is connected to multiple 2d face circumcenters (ggeo_make1D2Dembeddedlinks_dll)
        ///
        /// \see meshkernel::Contacts::ComputeMultipleContacts
        /// @param[in] meshKernelId  The id of the mesh state
        /// @param[in] oneDNodeMask  The mask to apply to 1d nodes (1 = generate a connection, 0 = do not generate a connection)
        /// @return Error code
        MKERNEL_API int mkernel_contacts_compute_multiple(int meshKernelId,
                                                          const int* oneDNodeMask);

        /// @brief Computes 1d-2d contacts, where a 2d face per polygon is connected to the closest 1d node (ggeo_make1D2Droofgutterpipes_dll)
        ///
        /// \see meshkernel::Contacts::ComputeContactsWithPolygons
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] oneDNodeMask The mask to apply to 1d nodes (1 = generate a connection, 0 = do not generate a connection)
        /// @param[in] polygons     The polygons to connect
        /// @return Error code
        MKERNEL_API int mkernel_contacts_compute_with_polygons(int meshKernelId,
                                                               const int* oneDNodeMask,
                                                               const GeometryList& polygons);

        /// @brief Computes 1d-2d contacts, where 1d nodes are connected to the 2d faces mass centers containing the input point (ggeo_make1D2Dstreetinletpipes_dll)
        ///
        /// \see meshkernel::Contacts::ComputeContactsWithPoints
        /// @param[in] meshKernelId  The id of the mesh state
        /// @param[in] oneDNodeMask  The mask to apply to 1d nodes (1 = generate a connection, 0 = do not generate a connection)
        /// @param[in] points        The points selecting the faces to connect
        /// @return Error code
        MKERNEL_API int mkernel_contacts_compute_with_points(int meshKernelId,
                                                             const int* oneDNodeMask,
                                                             const GeometryList& points);

        /// @brief Computes 1d-2d contacts, where 1d nodes are connected to the closest 2d faces at the boundary (ggeo_make1D2DRiverLinks_dll)
        ///
        /// \see meshkernel::Contacts::ComputeBoundaryContacts
        /// @param[in] meshKernelId The id of the mesh state.
        /// @param[in] oneDNodeMask The mask to apply to 1d nodes (1 = generate a connection, 0 = do not generate a connection)
        /// @param[in] polygons     The points selecting the faces to connect.
        /// @param[in] searchRadius The radius used for searching neighboring faces, if equal to constants::missing::doubleValue, the search radius will be calculated internally.
        /// @return Error code
        MKERNEL_API int mkernel_contacts_compute_boundary(int meshKernelId,
                                                          const int* oneDNodeMask,
                                                          const GeometryList& polygons,
                                                          double searchRadius);

        /// @brief Directional curvilinear grid refinement. Additional gridlines are added perpendicularly to the segment defined by lowerLeftCorner and xUpperRightCorner.
        ///
        /// \p firstPoint and \p secondPoint must lie on the same grid line.
        /// @param[in] meshKernelId      The id of the mesh state.
        /// @param[in] xLowerLeftCorner  The x coordinate of the lower left corner of the block to refine
        /// @param[in] yLowerLeftCorner  The y coordinate of the lower left corner of the block to refine
        /// @param[in] xUpperRightCorner The x coordinate of the upper right corner of the block to refine
        /// @param[in] yUpperRightCorner The y coordinate of the upper right corner of the block to refine
        /// @param[in] refinement        The number of grid lines to add between \p firstPoint and \p secondPoint
        /// @return                            Error code
        MKERNEL_API int mkernel_curvilinear_refine(int meshKernelId,
                                                   double xLowerLeftCorner,
                                                   double yLowerLeftCorner,
                                                   double xUpperRightCorner,
                                                   double yUpperRightCorner,
                                                   int refinement);

        /// @brief Directional curvilinear grid de-refinement. Grid lines are removed perpendicularly to the segment defined by lowerLeftCorner and xUpperRightCorner.
        ///
        /// \p firstPoint and \p secondPoint must lie on the same grid line.
        /// @param meshKernelId          The id of the mesh state.
        /// @param[in] xLowerLeftCorner  The x coordinate of the lower left corner of the block to de-refine
        /// @param[in] yLowerLeftCorner  The y coordinate of the lower left corner of the block to de-refine
        /// @param[in] xUpperRightCorner The x coordinate of the upper right corner of the block to de-refine
        /// @param[in] yUpperRightCorner The y coordinate of the upper right corner of the block to de-refine
        /// @return Error code
        MKERNEL_API int mkernel_curvilinear_derefine(int meshKernelId,
                                                     double xLowerLeftCorner,
                                                     double yLowerLeftCorner,
                                                     double xUpperRightCorner,
                                                     double yUpperRightCorner);

        /// @brief Generates curvilinear grid from splines with transfinite interpolation
        /// @param[in] meshKernelId          The id of the mesh state
        /// @param[in] splines               The splines to use for curvilinear grid generation
        /// @param[in] curvilinearParameters The curvilinear parameters
        /// @returns                         Error code
        MKERNEL_API int mkernel_curvilinear_compute_transfinite_from_splines(int meshKernelId,
                                                                             const GeometryList& splines,
                                                                             const meshkernel::CurvilinearParameters& curvilinearParameters);

        /// @brief Computes a curvilinear mesh in a polygon. 3 separate polygon nodes need to be selected.
        /// @param[in] meshKernelId  The id of the mesh state
        /// @param[in] polygons      The input polygons
        /// @param[in] firstNode     The first selected node
        /// @param[in] secondNode    The second selected node
        /// @param[in] thirdNode     The third selected node
        /// @param[in] useFourthSide Use the fourth polygon side to compute the curvilinear grid (0 no, 1 yes)
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_compute_transfinite_from_polygon(int meshKernelId,
                                                                             const GeometryList& polygons,
                                                                             int firstNode,
                                                                             int secondNode,
                                                                             int thirdNode,
                                                                             int useFourthSide);

        /// @brief Computes a curvilinear mesh in a triangle. 3 separate polygon nodes need to be selected.
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] polygon      The input polygons
        /// @param[in] firstNode    The first selected node
        /// @param[in] secondNode   The second selected node
        /// @param[in] thirdNode    The third node
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_compute_transfinite_from_triangle(int meshKernelId,
                                                                              const GeometryList& polygon,
                                                                              int firstNode,
                                                                              int secondNode,
                                                                              int thirdNode);

        /// @brief Generates curvilinear grid from splines with the advancing front method.
        /// @param[in] meshKernelId                   The id of the mesh state
        /// @param[in] geometryList                   The input splines corners
        /// @param[in] curvilinearParameters          The input parameters to generate the curvilinear grid
        /// @param[in] splinesToCurvilinearParameters The parameters of the advancing front algorithm
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_compute_orthogonal_grid_from_splines(int meshKernelId,
                                                                                 const GeometryList& geometryList,
                                                                                 const meshkernel::CurvilinearParameters& curvilinearParameters,
                                                                                 const meshkernel::SplinesToCurvilinearParameters& splinesToCurvilinearParameters);

        /// @brief Generates a curvilinear grid from splines with the advancing front method. Initialization step (interactive)
        /// @param[in] meshKernelId                    The id of the mesh state
        /// @param[in] geometryList                    The input splines corners
        /// @param[in] curvilinearParameters           The input parameters to generate the curvilinear grid
        /// @param[in] splinesToCurvilinearParameters  The parameters of the advancing front algorithm
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_initialize_orthogonal_grid_from_splines(int meshKernelId,
                                                                                    const GeometryList& geometryList,
                                                                                    const meshkernel::CurvilinearParameters& curvilinearParameters,
                                                                                    const meshkernel::SplinesToCurvilinearParameters& splinesToCurvilinearParameters);

        /// @brief One advancement of the front in curvilinear grid from splines (interactive)
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] layer        The layer index
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_iterate_orthogonal_grid_from_splines(int meshKernelId, int layer);

        /// @brief Converts curvilinear grid to mesh and refreshes the state (interactive)
        /// @param[in] meshKernelId The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_refresh_orthogonal_grid_from_splines(int meshKernelId);

        /// @brief Finalizes curvilinear grid from splines algorithm
        /// @param[in] meshKernelId The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_delete_orthogonal_grid_from_splines(int meshKernelId);

        /// @brief Makes a new curvilinear grid. If polygons is not empty, the first polygon will be used
        /// @param[in] meshKernelId       The id of the mesh state
        /// @param[in] makeGridParameters The structure containing the make grid parameters
        /// @param[in] geometryList       The polygons to account for
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_make_uniform(int meshKernelId,
                                                         const meshkernel::MakeGridParameters& makeGridParameters,
                                                         const GeometryList& geometryList);

        /// @brief Makes a new curvilinear grid. If polygons is not empty, the first polygon will be used
        /// @param[in] meshKernelId       The id of the mesh state
        /// @param[in] makeGridParameters The structure containing the make grid parameters
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_make_uniform_on_extension(int meshKernelId,
                                                                      const meshkernel::MakeGridParameters& makeGridParameters);

        /// @brief Initializes the orthogonal curvilinear algorithm
        /// @param[in] meshKernelId                The id of the mesh state
        /// @param[in] orthogonalizationParameters The orthogonalization parameters to use in the algorithm
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_initialize_orthogonalize(int meshKernelId,
                                                                     const meshkernel::OrthogonalizationParameters& orthogonalizationParameters);

        /// @brief Freezes a line in the curvilinear orthogonalization process
        /// @param[in] meshKernelId        The id of the mesh state
        /// @param[in] xFirstGridLineNode  The x coordinate of the first point of the line to freeze
        /// @param[in] yFirstGridLineNode  The y coordinate of the first point of the line to freeze
        /// @param[in] xSecondGridLineNode The x coordinate of the second point of the line to freeze
        /// @param[in] ySecondGridLineNode The y coordinate of the second point of the line to freeze
        /// @returns  Error code
        MKERNEL_API int mkernel_curvilinear_set_frozen_lines_orthogonalize(int meshKernelId,
                                                                           double xFirstGridLineNode,
                                                                           double yFirstGridLineNode,
                                                                           double xSecondGridLineNode,
                                                                           double ySecondGridLineNode);

        /// @brief Define a block on the curvilinear grid where to perform orthogonalization
        /// @param[in] meshKernelId      The id of the mesh state
        /// @param[in] xLowerLeftCorner  The x coordinate of the lower left corner of the block to orthogonalize
        /// @param[in] yLowerLeftCorner  The y coordinate of the lower left corner of the block to orthogonalize
        /// @param[in] xUpperRightCorner The x coordinate of the upper right corner of the block to orthogonalize
        /// @param[in] yUpperRightCorner The y coordinate of the upper right corner of the block to orthogonalize
        /// @return  Error code
        MKERNEL_API int mkernel_curvilinear_set_block_orthogonalize(int meshKernelId,
                                                                    double xLowerLeftCorner,
                                                                    double yLowerLeftCorner,
                                                                    double xUpperRightCorner,
                                                                    double yUpperRightCorner);

        /// @brief Orthogonalize a curvilinear grid
        /// @param[in] meshKernelId       The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_orthogonalize(int meshKernelId);

        /// @brief Resets the CurvilinearGridOrthogonalization instance in MeshKernelState
        /// @param[in] meshKernelId The id of the mesh state
        /// @return  Error code
        MKERNEL_API int mkernel_curvilinear_finalize_orthogonalize(int meshKernelId);

        /// @brief Smooths a curvilinear grid
        /// @param[in] meshKernelId        The id of the mesh state
        /// @param[in] smoothingIterations The number of smoothing iterations to perform
        /// @param[in] xLowerLeftCorner    The x coordinate of the lower left corner of the block to smooth
        /// @param[in] yLowerLeftCorner    The y coordinate of the lower left corner of the block to smooth
        /// @param[in] xUpperRightCorner   The x coordinate of the right corner of the block to smooth
        /// @param[in] yUpperRightCorner   The y coordinate of the upper right corner of the block to smooth
        /// @return Error code
        MKERNEL_API int mkernel_curvilinear_smoothing(int meshKernelId,
                                                      int smoothingIterations,
                                                      double xLowerLeftCorner,
                                                      double yLowerLeftCorner,
                                                      double xUpperRightCorner,
                                                      double yUpperRightCorner);

        /// @brief Smooths a curvilinear grid along the direction specified by a segment
        /// @param[in] meshKernelId                  The id of the mesh state
        /// @param[in] smoothingIterations           The number of smoothing iterations to perform
        /// @param[in] xFirstGridlineNode            The x coordinate of the first curvilinear grid node
        /// @param[in] yFirstGridlineNode            The y coordinate of the first curvilinear grid node
        /// @param[in] xSecondGridLineNode           The x coordinate of the second curvilinear grid node
        /// @param[in] ySecondGridLineNode           The y coordinate of the second curvilinear grid node
        /// @param[in] xLowerLeftCornerSmoothingArea The x coordinate of the lower left corner of the smoothing area
        /// @param[in] yLowerLeftCornerSmoothingArea The x coordinate of the lower left corner of the smoothing area
        /// @param[in] xUpperRightCornerSmootingArea The x coordinate of the upper right corner of the smoothing area
        /// @param[in] yUpperRightCornerSmootingArea The x coordinate of the upper right corner of the smoothing area
        /// @return Error code
        MKERNEL_API int mkernel_curvilinear_smoothing_directional(int meshKernelId,
                                                                  int smoothingIterations,
                                                                  double xFirstGridlineNode,
                                                                  double yFirstGridlineNode,
                                                                  double xSecondGridLineNode,
                                                                  double ySecondGridLineNode,
                                                                  double xLowerLeftCornerSmoothingArea,
                                                                  double yLowerLeftCornerSmoothingArea,
                                                                  double xUpperRightCornerSmootingArea,
                                                                  double yUpperRightCornerSmootingArea);

        /// @brief Initializes the curvilinear line shift algorithm
        /// @param[in] meshKernelId The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_initialize_line_shift(int meshKernelId);

        /// @brief Sets the start and end nodes of the line to shift
        /// @param[in] meshKernelId        The id of the mesh state
        /// @param[in] xFirstGridLineNode  The x coordinate of the first curvilinear grid node to shift
        /// @param[in] yFirstGridLineNode  The y coordinate of the first curvilinear grid node to shift
        /// @param[in] xSecondGridLineNode The x coordinate of the second curvilinear grid node to shift
        /// @param[in] ySecondGridLineNode The y coordinate of the second curvilinear grid node to shift
        /// @returns  Error code
        MKERNEL_API int mkernel_curvilinear_set_line_line_shift(int meshKernelId,
                                                                double xFirstGridLineNode,
                                                                double yFirstGridLineNode,
                                                                double xSecondGridLineNode,
                                                                double ySecondGridLineNode);

        /// @brief Defines a block on the curvilinear where the shifting is distributed
        /// @param[in] meshKernelId      The id of the mesh state
        /// @param[in] xLowerLeftCorner  The x coordinate of the lower left corner of the block
        /// @param[in] yLowerLeftCorner  The y coordinate of the lower left corner of the block
        /// @param[in] xUpperRightCorner The x coordinate of the upper right corner of the block
        /// @param[in] yUpperRightCorner The y coordinate of the upper right corner of the block
        /// @return  Error code
        MKERNEL_API int mkernel_curvilinear_set_block_line_shift(int meshKernelId,
                                                                 double xLowerLeftCorner,
                                                                 double yLowerLeftCorner,
                                                                 double xUpperRightCorner,
                                                                 double yUpperRightCorner);

        /// @brief Moves a node of the line to shift, the operation can be performed multiple times.
        /// @param[in] meshKernelId    The id of the mesh state
        /// @param[in] xFromCoordinate The x coordinate of the node to move (the closest curvilinear grid node will be found)
        /// @param[in] yFromCoordinate The y coordinate of the node to move (the closest curvilinear grid node will be found)
        /// @param[in] xToCoordinate   The x coordinate of the new node position
        /// @param[in] yToCoordinate   The y coordinate of the new node position
        /// @return  Error code
        MKERNEL_API int mkernel_curvilinear_move_node_line_shift(int meshKernelId,
                                                                 double xFromCoordinate,
                                                                 double yFromCoordinate,
                                                                 double xToCoordinate,
                                                                 double yToCoordinate);

        /// @brief Computes the new grid, shifting the line towards the moved nodes and distributing the shifting in block specified before
        /// @param[in] meshKernelId The id of the mesh state
        /// @return  Error code
        MKERNEL_API int mkernel_curvilinear_line_shift(int meshKernelId);

        /// @brief Resets the instance of the line shift algorithm in MeshKernelState
        /// @param[in] meshKernelId The id of the mesh state
        /// @return  Error code
        MKERNEL_API int mkernel_curvilinear_finalize_line_shift(int meshKernelId);

        /// @brief Inserts a new face on a curvilinear grid. The new face will be inserted on top of the closest edge by linear extrapolation.
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] xCoordinate  The x coordinate of the point used for finding the closest face.
        /// @param[in] yCoordinate  The y coordinate of the point used for finding the closest face.
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_insert_face(int meshKernelId, double xCoordinate, double yCoordinate);

        /// @brief Converts a curvilinear grid to an unstructured mesh
        MKERNEL_API int mkernel_curvilinear_convert_to_mesh2d(int meshKernelId);

        /// @brief Attracts/repulses grid lines in a block towards another set grid line
        /// @param[in] meshKernelId         The id of the mesh state
        /// @param[in] repulsionParameter   The attraction/repulsion parameter. If positive the grid lines will be attracted towards the set line, if negative the lines will be repulsed
        /// @param[in] xFirstNodeOnTheLine  The x coordinate of the first node of the set line
        /// @param[in] yFirstNodeOnTheLine  The y coordinate of the first node of the set line
        /// @param[in] xSecondNodeOnTheLine The x coordinate of the second node of the set line
        /// @param[in] ySecondNodeOnTheLine The y coordinate of the second node of the set line
        /// @param[in] xLowerLeftCorner     The x coordinate of the lower left corner of the block where the operation is performed
        /// @param[in] yLowerLeftCorner     The y coordinate of the lower left corner of the block where the operation is performed
        /// @param[in] xUpperRightCorner    The x coordinate of the upper right corner of the block where the operation is performed
        /// @param[in] yUpperRightCorner    The y coordinate of the upper right corner of the block where the operation is performed
        /// @return  Error code
        MKERNEL_API int mkernel_curvilinear_line_attraction_repulsion(int meshKernelId,
                                                                      double repulsionParameter,
                                                                      double xFirstNodeOnTheLine,
                                                                      double yFirstNodeOnTheLine,
                                                                      double xSecondNodeOnTheLine,
                                                                      double ySecondNodeOnTheLine,
                                                                      double xLowerLeftCorner,
                                                                      double yLowerLeftCorner,
                                                                      double xUpperRightCorner,
                                                                      double yUpperRightCorner);

        /// @brief Mirrors a boundary gridline outwards. The boundary grid line is defined by its starting and ending points
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] mirroringFactor       The mirroring factor
        /// @param[in] xFirstGridLineNode    The x coordinate of the first grid line point
        /// @param[in] yFirstGridLineNode    The y coordinate of the first grid line point
        /// @param[in] xSecondGridLineNode   The x coordinate of the second grid line point
        /// @param[in] ySecondGridLineNode   The y coordinate of the second grid line point
        /// @return Error code
        MKERNEL_API int mkernel_curvilinear_line_mirror(int meshKernelId,
                                                        double mirroringFactor,
                                                        double xFirstGridLineNode,
                                                        double yFirstGridLineNode,
                                                        double xSecondGridLineNode,
                                                        double ySecondGridLineNode);

        /// @brief Delete the node closest to a point
        /// @param meshKernelId The id of the mesh state
        /// @param[in] xPointCoordinate The x coordinate of the point
        /// @param[in] yPointCoordinate The y coordinate of the point
        /// @return  Error code
        MKERNEL_API int mkernel_curvilinear_delete_node(int meshKernelId,
                                                        double xPointCoordinate,
                                                        double yPointCoordinate);

        /// @brief Moves a point of a curvilinear grid from one location to another
        /// @param meshKernelId The id of the mesh state
        /// @param[in] xFromPoint The x coordinate of point to move
        /// @param[in] yFromPoint The y coordinate of point to move
        /// @param[in] xToPoint The new x coordinate of the point
        /// @param[in] yToPoint The new y coordinate of the point
        /// @return Error code
        MKERNEL_API int mkernel_curvilinear_move_node(int meshKernelId,
                                                      double xFromPoint,
                                                      double yFromPoint,
                                                      double xToPoint,
                                                      double yToPoint);

        /// @brief Gets the double value used in the back-end library as separator and missing value
        /// @return The double missing value used in mesh kernel
        MKERNEL_API double mkernel_get_separator();

        /// @brief Gets the double value used to separate the inner part of a polygon from its outer part
        /// @return The double missing value used in mesh kernel
        MKERNEL_API double mkernel_get_inner_outer_separator();

        /// @brief Triangle interpolation (ec_module)
        ///
        /// \see meshkernel::TriangulationInterpolation
        /// @param[in]  meshKernelId       The id of the mesh state
        /// @param[in]  samples            The samples coordinates and values
        /// @param[in]  locationType       The location type (see \ref meshkernel::Mesh::Location enum)
        /// @param[in]  results            The interpolation results with x and y coordinates
        /// @return Error code
        MKERNEL_API int mkernel_mesh2d_triangulation_interpolation(int meshKernelId,
                                                                   const GeometryList& samples,
                                                                   int locationType,
                                                                   GeometryList& results);

        /// @brief AveragingInterpolation interpolation (ec_module)
        ///
        /// \see meshkernel::AveragingInterpolation
        /// @param[in] meshKernelId           The id of the mesh state
        /// @param[in] samples                The samples coordinates and values
        /// @param[in] locationType           The location type (see \ref meshkernel::Mesh::Location enum)
        /// @param[in] averagingMethodType    The averaging method (see Method enum)
        /// @param[in] relativeSearchSize     The relative search size around the location (larger increases the number of samples considered)
        /// @param[in] minNumSamples          The minimum number of samples used for some interpolation algorithms to perform a valid interpolation
        /// @param[in] results                The interpolation results with x and y coordinates
        /// @return Error code
        MKERNEL_API int mkernel_mesh2d_averaging_interpolation(int meshKernelId,
                                                               const GeometryList& samples,
                                                               int locationType,
                                                               int averagingMethodType,
                                                               double relativeSearchSize,
                                                               size_t minNumSamples,
                                                               GeometryList& results);

        /// @brief Gets pointer to error message.
        /// @param[out] error_message The pointer to the latest error message
        /// @returns Error code
        MKERNEL_API int mkernel_get_error(const char*& error_message);

        /// @brief Gets pointer to version string.
        /// @param[out] version Version string
        /// @returns Error code
        MKERNEL_API int mkernel_get_version(const char*& version);

        /// @brief Gets the index of the erroneous entity.
        /// @param[out] invalidIndex The index of the erroneous entity
        /// @param[out] type         The entity type (node, edge or face, see MeshLocations)
        /// @returns Error code
        MKERNEL_API int mkernel_get_geometry_error(int& invalidIndex, int& type);

        /// @brief Gets an int indicating the edge location type
        /// @param[out] type The int indicating the edge location type
        /// @returns Error code
        MKERNEL_API int mkernel_get_edges_location_type(int& type);

        /// @brief Gets an int indicating the node location type
        /// @param[out] type The int indicating the node location type
        /// @returns Error code
        MKERNEL_API int mkernel_get_nodes_location_type(int& type);

        /// @brief Gets an int indicating the faces location type
        /// @param[out] type The int indicating the face location type
        /// @returns Error code
        MKERNEL_API int mkernel_get_faces_location_type(int& type);

        /// @brief Gets an int indicating the simple averaging method type
        /// @param[out] method The int indicating the averaging method type
        /// @returns Error code
        MKERNEL_API int mkernel_get_averaging_method_simple_averaging(int& method);

        /// @brief Gets an int indicating the closest point averaging method type
        /// @param[out] method The int indicating the closest point averaging method type
        /// @returns Error code
        MKERNEL_API int mkernel_get_averaging_method_closest_point(int& method);

        /// @brief Gets an int indicating the max value averaging method type
        /// @param[out] method The int indicating the max value averaging method type
        /// @returns Error code
        MKERNEL_API int mkernel_get_averaging_method_max(int& method);

        /// @brief Gets an int indicating the minimum averaging method type
        /// @param[out] method The int indicating the minimum averaging method type
        /// @returns Error code
        MKERNEL_API int mkernel_get_averaging_method_min(int& method);

        /// @brief Gets an int indicating the inverse distance weights averaging method type
        /// @param[out] method The int indicating the inverse weight distance averaging method type
        /// @returns Error code
        MKERNEL_API int mkernel_get_averaging_method_inverse_distance_weighting(int& method);

        /// @brief Gets an int indicating the minimum absolute value averaging method type
        /// @param[out] method The int indicating the minimum absolute value averaging method type
        /// @returns Error code
        MKERNEL_API int mkernel_get_averaging_method_min_absolute_value(int& method);

        /// @brief Gets an int indicating the cartesian projection
        /// @param[out] projection The int indicating the cartesian projection
        /// @return Error code
        MKERNEL_API int mkernel_get_projection_cartesian(int& projection);

        /// @brief Gets an int indicating the spherical projection
        /// @param[out] projection The int indicating the spherical projection
        /// @return Error code
        MKERNEL_API int mkernel_get_projection_spherical(int& projection);

        /// @brief Gets an int indicating the spherical accurate projection
        /// @param[out] projection The int indicating the spherical accurate projection
        /// @return Error code
        MKERNEL_API int mkernel_get_projection_spherical_accurate(int& projection);

        /// @brief Gets the coordinate projection of the meshkernel state
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[out] projection The int indicating the projection type
        /// @return Error code
        MKERNEL_API int mkernel_get_projection(int meshKernelId, int& projection);

#ifdef __cplusplus
    }
#endif
} // namespace meshkernelapi
