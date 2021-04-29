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

#include <MeshKernelApi/Contacts.hpp>
#include <MeshKernelApi/CurvilinearGrid.hpp>
#include <MeshKernelApi/CurvilinearParameters.hpp>
#include <MeshKernelApi/GeometryList.hpp>
#include <MeshKernelApi/InterpolationParameters.hpp>
#include <MeshKernelApi/MakeMeshParameters.hpp>
#include <MeshKernelApi/Mesh1D.hpp>
#include <MeshKernelApi/Mesh2D.hpp>
#include <MeshKernelApi/OrthogonalizationParameters.hpp>
#include <MeshKernelApi/SampleRefineParameters.hpp>
#include <MeshKernelApi/SplinesToCurvilinearParameters.hpp>

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
        Exception = 1,
        InvalidGeometry = 2
    };

#ifdef __cplusplus
    extern "C"
    {
#endif
        /// @brief Creates a new mesh state and returns the generated \p meshKernelId
        /// @param[in] isGeographic  Cartesian (0) or spherical (1) mesh
        /// @param[out] meshKernelId The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_allocate_state(int isGeographic, int& meshKernelId);

        /// @brief Deallocate mesh state
        /// @param[in] meshKernelId The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_deallocate_state(int meshKernelId);

        /// @brief Deletes a mesh in a polygon using several options
        /// @param[in] meshKernelId   The id of the mesh state
        /// @param[in] polygon        The polygon where to perform the operation
        /// @param[in] deletionOption The deletion option \ref meshkernel::Mesh2D::DeleteMeshOptions
        /// @param[in] invertDeletion Inverts the deletion of selected features
        /// @returns Error code
        MKERNEL_API int mkernel_delete_mesh2d(int meshKernelId,
                                              const GeometryList& polygon,
                                              int deletionOption,
                                              bool invertDeletion);

        /// @brief Sets the meshkernel::Mesh2D state
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] mesh2d       The Mesh2D data
        /// @returns Error code
        MKERNEL_API int mkernel_set_mesh2d(int meshKernelId,
                                           const Mesh2D& mesh2d);

        /// @brief Sets the meshkernel::Mesh1D state
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] mesh1d       The Mesh1D data
        /// @returns Error code
        MKERNEL_API int mkernel_set_mesh1d(int meshKernelId,
                                           const Mesh1D& mesh1d);

        /// @brief Gets the Mesh2D dimensions
        ///
        /// The integer parameters of the Mesh2D struct are set to the corresponding dimensions
        /// The pointers are set to null, and must be set to correctly sized memory
        /// before passing the struct to `mkernel_get_data_mesh2d`.
        /// @param[in]  meshKernelId The id of the mesh state
        /// @param[out] mesh2d       The Mesh2D dimensions
        /// @returns Error code
        MKERNEL_API int mkernel_get_dimensions_mesh2d(int meshKernelId,
                                                      Mesh2D& mesh2d);

        /// @brief Gets the Mesh2D dimensions data
        ///
        /// This function ought to be called after `mkernel_get_dimensions_mesh2d` has been called
        /// and the pointers have been set to correctly sized memory.
        /// @param[in]     meshKernelId The id of the mesh state
        /// @param[in,out] mesh2d       The Mesh2D data
        /// @returns Error code
        MKERNEL_API int mkernel_get_data_mesh2d(int meshKernelId, Mesh2D& mesh2d);

        /// @brief Gets the curvilinear grid dimensions as a CurvilinearGrid struct (converted as set of edges and nodes)
        ///
        /// The integer parameters of the CurvilinearGrid struct are set to the corresponding dimensions
        /// The pointers are set to null, and must be set to correctly sized memory
        /// before passing the struct to `mkernel_get_data_curvilinear`.
        /// @param[in]  meshKernelId    The id of the mesh state.
        /// @param[out] curvilinearGrid The structure containing the dimensions of the curvilinear grid.
        /// @returns Error code
        MKERNEL_API int mkernel_get_dimensions_curvilinear(int meshKernelId, CurvilinearGrid& curvilinearGrid);

        /// @brief Gets the curvilinear grid data as a CurvilinearGrid struct (converted as set of edges and nodes)
        ///
        /// This function ought to be called after `mkernel_get_curvilinear_dimension` has been called
        /// and the pointers have been set to correctly sized memory.
        /// @param[in]  meshKernelId    The id of the mesh state
        /// @param[out] curvilinearGrid The structure containing the curvilinear grid arrays.
        /// @returns Error code
        MKERNEL_API int mkernel_get_data_curvilinear(int meshKernelId, CurvilinearGrid& curvilinearGrid);

        /// @brief Gets the Mesh1D dimensions
        ///
        /// The integer parameters of the Mesh1D struct are set to the corresponding dimensions
        /// The pointers are set to null, and must be set to correctly sized memory
        /// before passing the struct to `mkernel_get_dimensions_mesh1d`
        /// @param[in]  meshKernelId The id of the mesh state
        /// @param[out] mesh1d       The structure containing the dimensions of the Mesh1D.
        /// @returns Error code
        MKERNEL_API int mkernel_get_dimensions_mesh1d(int meshKernelId, Mesh1D& mesh1d);

        /// @brief Gets the Mesh1D dimensions data
        ///
        /// This function ought to be called after `mkernel_get_dimensions_mesh1d` has been called
        /// and the pointers have been set to correctly sized memory
        /// @param[in]     meshKernelId The id of the mesh state
        /// @param[in,out] mesh1d       The structure containing the Mesh1D arrays.
        /// @returns Error code
        MKERNEL_API int mkernel_get_data_mesh1d(int meshKernelId, Mesh1D& mesh1d);

        /// @brief Gets the number of contacts
        /// @param[in]  meshKernelId           The id of the mesh state
        /// @param[out] contacts               Contacts data
        /// @returns                           Error code
        MKERNEL_API int mkernel_get_dimensions_contacts(int meshKernelId, Contacts& contacts);

        /// @brief Gets the contacts indices (from index / to indices)
        /// @param[in]  meshKernelId           The id of the mesh state
        /// @param[out] contacts               Contacts data
        /// @returns                           Error code
        MKERNEL_API int mkernel_get_data_contacts(int meshKernelId, Contacts& contacts);

        /// @brief Count the number of hanging edges in a mesh2d.
        /// An hanging edge is an edge where one of the two nodes is not connected.
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[out] numEdges The number of hanging edges
        /// @returns Error code
        MKERNEL_API int mkernel_count_hanging_edges_mesh2d(int meshKernelId, int& numEdges);

        /// @brief Gets the indices of hanging edges. An hanging edge is an edge where one of the two nodes is not connected.
        /// @param[in]     meshKernelId        The id of the mesh state
        /// @param[in,out] edges Pointer to memory where the indices of the hanging edges will be stored
        /// @returns Error code
        MKERNEL_API int mkernel_get_hanging_edges_mesh2d(int meshKernelId, int* edges);

        /// @brief Deletes all hanging edges. An hanging edge is an edge where one of the two nodes is not connected.
        /// @param[in] meshKernelId The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_delete_hanging_edges_mesh2d(int meshKernelId);

        /// @brief Mesh2d orthogonalization.
        ///
        /// The function modifies the mesh for achieving orthogonality between the edges and the segments connecting the face circumcenters.
        /// The amount of orthogonality is traded against the mesh smoothing (in this case the equality of face areas).
        /// The parameter to regulate the amount of orthogonalization is contained in  \ref meshkernelapi::OrthogonalizationParameters::orthogonalization_to_smoothing_factor
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] projectToLandBoundaryOption The option to determine how to snap to land boundaries
        /// @param[in] orthogonalizationParameters The structure containing the orthogonalization parameters \ref meshkernelapi::OrthogonalizationParameters
        /// @param[in] polygons                    The polygon where to perform the orthogonalization
        /// @param[in] landBoundaries              The land boundaries to account for in the orthogonalization process
        /// @returns Error code
        MKERNEL_API int mkernel_compute_orthogonalization_mesh2d(int meshKernelId,
                                                                 int projectToLandBoundaryOption,
                                                                 const OrthogonalizationParameters& orthogonalizationParameters,
                                                                 const GeometryList& polygons,
                                                                 const GeometryList& landBoundaries);

        /// @brief Initialization of the \ref meshkernel::OrthogonalizationAndSmoothing algorithm.
        ///
        /// This is the first function to call when using orthogonalization in interactive mode (visualizing the grid while it is orthogonalizing),
        /// in order to set the internal state of the algorithm reused during the iterations.
        /// @param[in] meshKernelId                The id of the mesh state
        /// @param[in] projectToLandBoundaryOption The option to determine how to snap to land boundaries
        /// @param[in] orthogonalizationParameters The structure containing the user defined orthogonalization parameters
        /// @param[in] geometryListPolygon         The polygon where to perform the orthogonalization
        /// @param[in] geometryListLandBoundaries  The land boundaries to account for in the orthogonalization process
        /// @returns Error code
        MKERNEL_API int mkernel_initialize_orthogonalization_mesh2d(int meshKernelId,
                                                                    int projectToLandBoundaryOption,
                                                                    OrthogonalizationParameters& orthogonalizationParameters,
                                                                    const GeometryList& geometryListPolygon,
                                                                    const GeometryList& geometryListLandBoundaries);

        /// @brief Prepares an outer orthogonalization iteration, computing the new orthogonalization and smoothing weights from the modified mesh geometry (in interactive mode).
        ///
        /// `mkernel_initialize_orthogonalization_mesh2d` function must be called before.
        /// @param[in] meshKernelId The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_prepare_outer_iteration_orthogonalization_mesh2d(int meshKernelId);

        /// @brief Performs inner orthogonalization iteration, by slowly moving the mesh nodes to new optimal positions (interactive mode).
        ///
        /// `mkernel_prepare_outer_iteration_orthogonalization_mesh2d` function must be called before.
        /// @param[in] meshKernelId The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_compute_inner_ortogonalization_iteration_mesh2d(int meshKernelId);

        /// @brief Finalizes the orthogonalization outer iteration, computing the new coefficients for grid adaption and the new face circumcenters (interactive mode).
        /// @param[in] meshKernelId  The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_finalize_inner_ortogonalization_iteration_mesh2d(int meshKernelId);

        /// @brief Cleans the orthogonalization algorithm state, allocated in `mkernel_initialize_orthogonalization_mesh2d` (interactive mode)
        /// @param[in] meshKernelId The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_delete_orthogonalization_mesh2d(int meshKernelId);

        /// @brief Gets the mesh orthogonality, expressed as the ratio between the edges and the segments connecting the face circumcenters.
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[out] geometryList The orthogonality values of each edge
        /// @returns Error code
        MKERNEL_API int mkernel_get_orthogonality_mesh2d(int meshKernelId, GeometryList& geometryList);

        /// @brief Gets the smoothness, expressed as the ratio between the values of two neighboring faces areas.
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[out] geometryList The smoothness values at each edge
        /// @returns Error code
        MKERNEL_API int mkernel_get_smoothness_mesh2d(int meshKernelId, GeometryList& geometryList);

        /// @brief Get the computed spline points between two corner nodes
        /// @param[in] geometryListIn The input corner nodes of the splines
        /// @param[out] geometryListOut The output spline
        /// @param[out] numberOfPointsBetweenNodes The number of spline points to generate between two corner nodes.
        /// @returns Error code
        MKERNEL_API int mkernel_get_splines(const GeometryList& geometryListIn,
                                            GeometryList& geometryListOut,
                                            int numberOfPointsBetweenNodes);

        /// @brief Gets the closest mesh2d node coordinates to a point, searching within a radius.
        /// @param[in]  meshKernelId    Id of the grid state
        /// @param[in]  point           The point coordinate
        /// @param[in]  searchRadius    The radii where to search for mesh nodes
        /// @param[out] node            The found Mesh2D node coordinates
        /// @returns Error code
        MKERNEL_API int mkernel_get_closest_node_mesh2d(int meshKernelId,
                                                        const GeometryList& point,
                                                        double searchRadius,
                                                        GeometryList& node);

        /// @brief Generates a triangular mesh2d grid within a polygon. The size of the triangles is determined from the length of the polygon edges.
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] polygonPoints The polygon where to triangulate
        /// @returns Error code
        MKERNEL_API int mkernel_make_mesh_from_polygon_mesh2d(int meshKernelId, const GeometryList& polygonPoints);

        /// @brief Makes a triangular mesh from  a set of samples, triangulating the sample points.
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] samples The samples where to triangulate
        /// @returns Error code
        MKERNEL_API int mkernel_make_mesh_from_samples_mesh2d(int meshKernelId, const GeometryList& samples);

        /// @brief Retrieves the boundaries of a mesh as a series of separated polygons.
        ///
        /// For example, if a mesh has an single inner hole, two polygons will be generated, one for the inner boundary and one for the outer boundary.
        /// @param[in]  meshKernelId The id of the mesh state
        /// @param[out] boundaryPolygons The output network boundary polygon
        /// @returns Error code
        MKERNEL_API int mkernel_get_mesh_boundaries_as_polygons_mesh2d(int meshKernelId, GeometryList& boundaryPolygons);

        /// @brief Counts the number of polygon nodes contained in the mesh boundary polygons computed in function `mkernel_get_mesh_boundaries_as_polygons_mesh2d`
        /// @param[in]  meshKernelId         The id of the mesh state
        /// @param[out] numberOfPolygonNodes The number of polygon nodes
        /// @returns Error code
        MKERNEL_API int mkernel_count_mesh_boundaries_as_polygon_mesh2d(int meshKernelId, int& numberOfPolygonNodes);

        /// @brief Refines the polygon perimeter between two nodes. This interval is refined to achieve a target edge length.
        ///
        /// The function is often used before `mkernel_make_mesh_from_polygon_mesh2d`, for generating a triangular mesh where edges have a desired length.
        /// @param[in]  meshKernelId       The id of the mesh state
        /// @param[in]  polygonToRefine    The input polygon to refine
        /// @param[in]  firstNodeIndex     The first index of the refinement interval
        /// @param[in]  secondNodeIndex    The second index of the refinement interval
        /// @param[in]  targetEdgeLength   The target interval edge length
        /// @param[out] refinedPolygon     The refined polygon
        /// @returns Error code
        MKERNEL_API int mkernel_refine_polygon(int meshKernelId,
                                               const GeometryList& polygonToRefine,
                                               int firstNodeIndex,
                                               int secondNodeIndex,
                                               double targetEdgeLength,
                                               GeometryList& refinedPolygon);

        /// @brief Counts the number of polygon nodes resulting from polygon refinement with `mkernel_refine_polygon`.
        ///
        /// This function should be used by clients before `mkernel_refine_polygon` for allocating \ref GeometryList containing the refinement result.
        /// @param[in] meshKernelId   The id of the mesh state
        /// @param[in] geometryListIn The input polygon to refine
        /// @param[in] firstIndex     The first index of the refinement interval
        /// @param[in] secondIndex    The second index of the refinement interval
        /// @param[in] distance       The target interval edge length
        /// @param[out] numberOfPolygonNodes The number of nodes after refinement
        /// @returns Error code
        MKERNEL_API int mkernel_count_refine_polygon(int meshKernelId,
                                                     const GeometryList& geometryListIn,
                                                     int firstIndex,
                                                     int secondIndex,
                                                     double distance,
                                                     int& numberOfPolygonNodes);

        /// @brief Merges the mesh2d nodes within a distance of 0.001 m, effectively removing all small edges
        /// @param[in] meshKernelId   The id of the mesh state
        /// @param[in] geometryListIn The polygon defining the area where to perform will be performed
        /// @returns Error code
        MKERNEL_API int mkernel_merge_nodes_mesh2d(int meshKernelId, const GeometryList& geometryListIn);

        /// @brief Merges two mesh2d nodes in one.
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] firstNode    The index of the first node to merge
        /// @param[in] secondNode   The index of the second node to merge
        /// @returns Error code
        MKERNEL_API int mkernel_merge_two_nodes_mesh2d(int meshKernelId, int firstNode, int secondNode);

        /// @brief Gets the indices of the mesh2d nodes selected with a polygon.
        /// @param[in]  meshKernelId   The id of the mesh state
        /// @param[in]  geometryListIn The input polygon
        /// @param[in]  inside         Selection of the nodes inside the polygon (1) or outside (0)
        /// @param[out] selectedNodes  The selected nodes indices
        /// @returns Error code
        MKERNEL_API int mkernel_nodes_in_polygons_mesh2d(int meshKernelId,
                                                         const GeometryList& geometryListIn,
                                                         int inside,
                                                         int* selectedNodes);

        /// @brief Counts the number of selected mesh node indices.
        ///
        /// This function should be used by clients before `mkernel_nodes_in_polygons_mesh2d` for allocating an integer array storing the selection results.
        /// @param[in]  meshKernelId      The id of the mesh state
        /// @param[in]  geometryListIn    The input polygon
        /// @param[in]  inside            Selection of the nodes inside the polygon (1) or outside (0)
        /// @param[out] numberOfMeshNodes The number of selected nodes
        /// @returns Error code
        MKERNEL_API int mkernel_count_nodes_in_polygons_mesh2d(int meshKernelId,
                                                               const GeometryList& geometryListIn,
                                                               int inside,
                                                               int& numberOfMeshNodes);

        /// @brief Insert a new mesh2d edge connecting two nodes
        /// @param[in]  meshKernelId The id of the mesh state
        /// @param[in]  startNode    The index of the first node to connect
        /// @param[in]  endNode      The index of the second node to connect
        /// @param[out] newEdgeIndex The index of the new generated edge
        /// @returns Error code
        MKERNEL_API int mkernel_insert_edge_mesh2d(int meshKernelId, int startNode, int endNode, int& newEdgeIndex);

        /// @brief Insert a new mesh2d node at a specific coordinate.
        /// @param[in]  meshKernelId The id of the mesh state
        /// @param[in]  nodeCoordinate  The coordinate of the new node to insert
        /// @param[out] nodeIndex    The index of the new mesh node
        /// @returns Error code
        MKERNEL_API int mkernel_insert_node_mesh2d(int meshKernelId, GeometryList const& nodeCoordinate, int& nodeIndex);

        /// @brief Deletes a mesh2d node
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] nodeIndex    The index of the node to delete
        /// @returns Error code
        MKERNEL_API int mkernel_delete_node_mesh2d(int meshKernelId, int nodeIndex);

        /// @brief Moves a mesh2d node to a new position
        /// @param[in] meshKernelId    The id of the mesh state
        /// @param[in] newNodePosition The new coordinate of the moved node
        /// @param[in] nodeIndex       The index of the mesh2d node to be moved
        /// @returns Error code
        MKERNEL_API int mkernel_move_node_mesh2d(int meshKernelId, const GeometryList& newNodePosition, int nodeIndex);

        /// @brief Deletes the closest mesh2d edge to a point.
        /// The coordinates of the edge middle points are used for calculating the distances to the point.
        /// @param[in] meshKernelId   The id of the mesh state
        /// @param[in] point          The coordinate of the point
        /// @returns Error code
        MKERNEL_API int mkernel_delete_edge_mesh2d(int meshKernelId, const GeometryList& point);

        /// @brief Gets the closest mesh2d edge to a point.
        /// @param[in] meshKernelId   The id of the mesh state
        /// @param[in] point          The coordinate of the point
        /// @param[out] edgeIndex     The found edge index
        /// @returns Error code
        MKERNEL_API int mkernel_get_edge_mesh2d(int meshKernelId, const GeometryList& point, int& edgeIndex);

        /// @brief Generate a new polygon from an existing one by offsetting the perimeter by a given distance.
        ///
        /// Offsetting can be done inward or outward the existing polygon.
        /// @param[in] meshKernelId     The id of the mesh state
        /// @param[in] geometryListIn   The polygon to offset
        /// @param[in] inWard           Compute the inner offset (true) or outer offset offset (false)
        /// @param[in] distance         The offset distance
        /// @param[out] geometryListOut The resulting offset polygon
        /// @returns Error code
        MKERNEL_API int mkernel_get_offset_polygon(int meshKernelId,
                                                   const GeometryList& geometryListIn,
                                                   bool inWard,
                                                   double distance,
                                                   GeometryList& geometryListOut);

        /// @brief Counts the number of polygon nodes resulting from polygon offset.
        ///
        /// This function should be used by clients before `mkernel_get_offset_polygon` for allocating the \ref GeometryList containing the offset result.
        /// @param[in] meshKernelId          The id of the mesh state
        /// @param[in] geometryListIn        The polygon to offset
        /// @param[in] innerPolygon          Compute inner (true) or outer offset (false)
        /// @param[in] distance              The offset distance
        /// @param[out] numberOfPolygonNodes The number of nodes in the offset polygon
        /// @returns Error code
        MKERNEL_API int mkernel_count_offset_polygon(int meshKernelId,
                                                     const GeometryList& geometryListIn,
                                                     bool innerPolygon,
                                                     double distance,
                                                     int& numberOfPolygonNodes);

        /// @brief Refines a mesh2d based on samples. Refinement is achieved by successive splits of the face edges.
        ///
        /// The number of successive splits is indicated on the sample value.
        /// For example a value of 0 means no split and hence no refinement, a value of 1 a single split (a quadrilateral face generates 4 faces),
        /// a value of 2 two splits (a quadrilateral face generates 16 faces).
        /// @param[in] meshKernelId            The id of the mesh state
        /// @param[in] samples                 The sample set
        /// @param[in] interpolationParameters The interpolation parameters
        /// @param[in] sampleRefineParameters  The interpolation settings related to the samples
        /// @returns Error code
        MKERNEL_API int mkernel_refine_based_on_samples_mesh2d(int meshKernelId,
                                                               const GeometryList& samples,
                                                               const InterpolationParameters& interpolationParameters,
                                                               const SampleRefineParameters& sampleRefineParameters);

        /// @brief Refines a mesh2d within a polygon. Refinement is achieved by splitting the edges contained in the polygon by two.
        /// @param[in] meshKernelId            The id of the mesh state
        /// @param[in] geometryList            The closed polygon where to perform the refinement
        /// @param[in] interpolationParameters The interpolation parameters
        /// @returns Error code
        MKERNEL_API int mkernel_refine_based_on_polygon_mesh2d(int meshKernelId,
                                                               const GeometryList& geometryList,
                                                               const InterpolationParameters& interpolationParameters);

        /// @brief Finds the mesh2d node closest to a point, within a search radius.
        /// @param[in] meshKernelId   The id of the mesh state
        /// @param[in] point          The coordinate of the point
        /// @param[in] searchRadius   The search radius
        /// @param[out] nodeIndex     The index of the found node
        /// @returns Error code
        MKERNEL_API int mkernel_get_node_index_mesh2d(int meshKernelId,
                                                      const GeometryList& point,
                                                      double searchRadius,
                                                      int& nodeIndex);

        /// @brief Selects the polygon nodes within another polygon.
        /// @param[in]  meshKernelId   The id of the mesh state
        /// @param[in]  selectingPolygon   The selecting polygon
        /// @param[in]  polygonToSelect    The polygon to select
        /// @param[out] selectionResults   The selection result, contained in the in the values field of \ref GeometryList (0.0 not selected, 1.0 selected).
        /// Note that the selection selectionResults variable must be allocated by the client.
        /// @returns Error code
        MKERNEL_API int mkernel_get_points_in_polygon(int meshKernelId,
                                                      const GeometryList& selectingPolygon,
                                                      const GeometryList& polygonToSelect,
                                                      GeometryList& selectionResults);

        /// @brief Flips mesh2d edges, to optimize the mesh smoothness. This operation is usually performed after `mkernel_refine_based_on_samples_mesh2d` or `mkernel_refine_based_on_polygon_mesh2d`.
        ///
        /// Nodes that are connected to more than six other nodes are typically enclosed by faces of highly non-uniform shape and wildly varying areas.
        /// @param[in] meshKernelId                The id of the mesh state
        /// @param[in] isTriangulationRequired     The option to triangulate also non triangular cells (if activated squares becomes triangles)
        /// @param[in] projectToLandBoundaryOption The option to determine how to snap to land boundaries
        /// @returns Error code
        MKERNEL_API int mkernel_flip_edges_mesh2d(int meshKernelId,
                                                  int isTriangulationRequired,
                                                  int projectToLandBoundaryOption);

        /// @brief Gets the number of obtuse mesh2d triangles. Obtuse triangles are those having one edge longer than the sum of the other two.
        /// @param[in]  meshKernelId       The id of the mesh state
        /// @param[out] numObtuseTriangles The number of obtuse triangles
        /// @return Error code
        MKERNEL_API int mkernel_count_obtuse_triangles_mesh2d(int meshKernelId, int& numObtuseTriangles);

        /// @brief Gets the mass centers of obtuse mesh2d triangles. Obtuse triangles are those having one edge longer than the sum of the other two.
        /// @param[in]  meshKernelId  The id of the mesh state
        /// @param[out] result        The coordinates of the obtuse triangles mass centers stored in coordinates_x and coordinates_y of a \ref GeometryList
        /// @return Error code
        MKERNEL_API int mkernel_get_obtuse_triangles_mass_centers_mesh2d(int meshKernelId, GeometryList& result);

        /// @brief Counts the number of small mesh2d flow edges. The flow edges are the edges connecting faces circumcenters.
        /// @param[in] meshKernelId                  The id of the mesh state
        /// @param[in] smallFlowEdgesLengthThreshold The configurable length for detecting a small flow edge
        /// @param[out] numSmallFlowEdges            The number of the small flow edges
        /// @return Error code
        MKERNEL_API int mkernel_count_small_flow_edge_centers_mesh2d(int meshKernelId,
                                                                     double smallFlowEdgesLengthThreshold,
                                                                     int& numSmallFlowEdges);

        /// @brief Gets the small mesh2d flow edges. The flow edges are the edges connecting faces circumcenters.
        /// @param[in] meshKernelId            The id of the mesh state
        /// @param[in] smallFlowEdgesThreshold The configurable threshold for detecting a small flow edge
        /// @param[out] result                 The middle points of the small flow edges, stored in coordinates_x and coordinates_y of a \ref GeometryList
        /// @return Error code
        MKERNEL_API int mkernel_get_small_flow_edge_centers_mesh2d(int meshKernelId,
                                                                   double smallFlowEdgesThreshold,
                                                                   GeometryList& result);

        /// @brief Deletes all small mesh2d flow edges and small triangles. The flow edges are the edges connecting faces circumcenters.
        /// @param[in] meshKernelId               The id of the mesh state
        /// @param[in] smallFlowEdgesThreshold    The configurable threshold for detecting the small flow edges
        /// @param[in] minFractionalAreaTriangles The ratio of the face area to the average area of neighboring non triangular faces.
        /// This parameter is used for determining if a triangular face is small.
        /// @return Error code
        MKERNEL_API int mkernel_delete_small_flow_edges_and_small_triangles_mesh2d(int meshKernelId,
                                                                                   double smallFlowEdgesThreshold,
                                                                                   double minFractionalAreaTriangles);

        /// @brief Computes 1d-2d contacts, where each single 1d node is connected to one mesh2d face circumcenter (ggeo_make1D2Dinternalnetlinks_dll)
        ///
        /// \see meshkernel::Contacts::ComputeSingleContacts
        /// @param[in]  meshKernelId  The id of the mesh state
        /// @param[in]  oneDNodeMask  The mask to apply to 1d nodes (1 = connect node, 0 = do not connect)
        /// @param[in]  polygons      The polygons selecting the area where the 1d-2d contacts will be generated.
        /// @return                   Error code
        MKERNEL_API int mkernel_compute_single_contacts(int meshKernelId,
                                                        const int* oneDNodeMask,
                                                        const GeometryList& polygons);

        /// @brief Computes 1d-2d contacts, where a single 1d node is connected to multiple 2d face circumcenters (ggeo_make1D2Dembeddedlinks_dll)
        ///
        /// \see meshkernel::Contacts::ComputeMultipleContacts
        /// @param[in]  meshKernelId  The id of the mesh state
        /// @param[in]  oneDNodeMask  The mask to apply to 1d nodes (1 = generate a connection, 0 = do not generate a connection)
        /// @return                   Error code
        MKERNEL_API int mkernel_compute_multiple_contacts(int meshKernelId,
                                                          const int* oneDNodeMask);

        /// @brief Computes 1d-2d contacts, where a 2d face per polygon is connected to the closest 1d node (ggeo_make1D2Droofgutterpipes_dll)
        ///
        /// \see meshkernel::Contacts::ComputeContactsWithPolygons
        /// @param[in]  meshKernelId  The id of the mesh state
        /// @param[in]  oneDNodeMask  The mask to apply to 1d nodes (1 = generate a connection, 0 = do not generate a connection)
        /// @param[in]  polygons      The polygons to connect
        /// @return                   Error code
        MKERNEL_API int mkernel_compute_with_polygons_contacts(int meshKernelId,
                                                               const int* oneDNodeMask,
                                                               const GeometryList& polygons);

        /// @brief Computes 1d-2d contacts, where 1d nodes are connected to the 2d faces mass centers containing the input point (ggeo_make1D2Dstreetinletpipes_dll)
        ///
        /// \see meshkernel::Contacts::ComputeContactsWithPoints
        /// @param[in]  meshKernelId  The id of the mesh state
        /// @param[in]  oneDNodeMask  The mask to apply to 1d nodes (1 = generate a connection, 0 = do not generate a connection)
        /// @param[in]  points        The points selecting the faces to connect
        /// @return                   Error code
        MKERNEL_API int mkernel_compute_with_points_contacts(int meshKernelId,
                                                             const int* oneDNodeMask,
                                                             const GeometryList& points);

        /// @brief Computes 1d-2d contacts, where 1d nodes are connected to the closest 2d faces at the boundary (ggeo_make1D2DRiverLinks_dll)
        ///
        /// \see meshkernel::Contacts::ComputeBoundaryContacts
        /// @param[in]  meshKernelId The id of the mesh state.
        /// @param[in]  oneDNodeMask The mask to apply to 1d nodes (1 = generate a connection, 0 = do not generate a connection)
        /// @param[in]  polygons     The points selecting the faces to connect.
        /// @param[in]  searchRadius The radius used for searching neighboring faces, if equal to doubleMissingValue, the search radius will be calculated internally.
        /// @return                  Error code
        MKERNEL_API int mkernel_compute_boundary_contacts(int meshKernelId,
                                                          const int* oneDNodeMask,
                                                          const GeometryList& polygons,
                                                          double searchRadius);

        /// @brief Directional curvilinear grid refinement. Additional gridlines are added perpendicularly to the segment defined by \p firstPoint and \p secondPoint.
        ///
        /// \p firstPoint and \p secondPoint must lie on the same grid line.
        /// @param[in] meshKernelId            The id of the mesh state.
        /// @param[in] firstPoint              The first point defining the refinement zone.
        /// @param[in] secondPoint             The second point defining the refinement zone.
        /// @param[in] refinement              The number of grid lines to add between \p firstPoint and \p secondPoint
        /// @return                            Error code
        MKERNEL_API int mkernel_refine_curvilinear(int meshKernelId,
                                                   const GeometryList& firstPoint,
                                                   const GeometryList& secondPoint,
                                                   int refinement);

        /// @brief Directional curvilinear grid derefinement. Grid lines are removed perpendicularly to the segment defined by \p firstPoint and \p secondPoint.
        ///
        /// \p firstPoint and \p secondPoint must lie on the same grid line.
        /// @param meshKernelId            The id of the mesh state.
        /// @param firstPoint              The first point defining the de-refinement zone.
        /// @param secondPoint             The second point defining the de-refinement zone.
        /// @return Error code
        MKERNEL_API int mkernel_derefine_curvilinear(int meshKernelId,
                                                     const GeometryList& firstPoint,
                                                     const GeometryList& secondPoint);

        /// @brief Generates curvilinear grid from splines with transfinite interpolation
        /// @param[in] meshKernelId          The id of the mesh state
        /// @param[in] splines               The splines
        /// @param[in] curvilinearParameters The curvilinear parameters
        /// @returns                         Error code
        MKERNEL_API int mkernel_compute_transfinite_from_splines_curvilinear(int meshKernelId,
                                                                             const GeometryList& splines,
                                                                             const CurvilinearParameters& curvilinearParameters);

        /// @brief Computes a curvilinear mesh in a polygon. 3 separate polygon nodes need to be selected.
        /// @param[in] meshKernelId  The id of the mesh state
        /// @param[in] polygons      The input polygons
        /// @param[in] firstNode     The first selected node
        /// @param[in] secondNode    The second selected node
        /// @param[in] thirdNode     The third node
        /// @param[in] useFourthSide Use (true/false) the fourth polygon side to compute the curvilinear grid
        /// @returns Error code
        MKERNEL_API int mkernel_compute_transfinite_from_polygon_curvilinear(int meshKernelId,
                                                                             const GeometryList& polygons,
                                                                             int firstNode,
                                                                             int secondNode,
                                                                             int thirdNode,
                                                                             bool useFourthSide);

        /// @brief Computes a curvilinear mesh in a triangle. 3 separate polygon nodes need to be selected.
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] polygon      The input polygons
        /// @param[in] firstNode    The first selected node
        /// @param[in] secondNode   The second selected node
        /// @param[in] thirdNode    The third node
        /// @returns Error code
        MKERNEL_API int mkernel_compute_transfinite_from_triangle_curvilinear(int meshKernelId,
                                                                              const GeometryList& polygon,
                                                                              int firstNode,
                                                                              int secondNode,
                                                                              int thirdNode);

        /// @brief Makes curvilinear grid from splines with an advancing front.
        /// @param[in] meshKernelId                   The id of the mesh state
        /// @param[in] geometryList                   The input splines corners
        /// @param[in] curvilinearParameters          The input parameters to generate the curvilinear grid
        /// @param[in] splinesToCurvilinearParameters The parameters of the advancing front algorithm
        /// @returns Error code
        MKERNEL_API int mkernel_compute_orthogonal_grid_from_splines_curvilinear(int meshKernelId,
                                                                                 const GeometryList& geometryList,
                                                                                 const CurvilinearParameters& curvilinearParameters,
                                                                                 const SplinesToCurvilinearParameters& splinesToCurvilinearParameters);

        /// @brief Generates a curvilinear grid from splines with the advancing front method. Initialization step (interactive)
        /// @param[in] meshKernelId                    The id of the mesh state
        /// @param[in] geometryList                    The input splines corners
        /// @param[in] curvilinearParameters           The input parameters to generate the curvilinear grid
        /// @param[in] splinesToCurvilinearParameters  The parameters of the advancing front algorithm
        /// @returns Error code
        MKERNEL_API int mkernel_initialize_orthogonal_grid_from_splines_curvilinear(int meshKernelId,
                                                                                    const GeometryList& geometryList,
                                                                                    const CurvilinearParameters& curvilinearParameters,
                                                                                    const SplinesToCurvilinearParameters& splinesToCurvilinearParameters);

        /// @brief One advancement of the front in curvilinear grid from splines (interactive)
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] layer        The layer index
        /// @returns Error code
        MKERNEL_API int mkernel_iterate_orthogonal_grid_from_splines_curvilinear(int meshKernelId, int layer);

        /// @brief Converts curvilinear grid to mesh and refreshes the state (interactive)
        /// @param[in] meshKernelId The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_refresh_orthogonal_grid_from_splines_curvilinear(int meshKernelId);

        /// @brief Finalizes curvilinear grid from splines algorithm
        /// @param[in] meshKernelId The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_delete_orthogonal_grid_from_splines_curvilinear(int meshKernelId);

        /// @brief Makes a new mesh
        /// @param[in] meshKernelId       The id of the mesh state
        /// @param[in] makeGridParameters The structure containing the make grid parameters
        /// @param[in] geometryList       The polygon to account for
        /// @returns Error code
        MKERNEL_API int mkernel_make_uniform_curvilinear(int meshKernelId,
                                                         const MakeMeshParameters& makeGridParameters,
                                                         const GeometryList& geometryList);

        /// @brief Initializes the orthogonal curvilinear algorithm
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] orthogonalizationParameters The orthogonalization parameters to use in the algorithm
        /// @returns Error code
        MKERNEL_API int mkernel_initialize_orthogonalize_curvilinear(int meshKernelId,
                                                                     const OrthogonalizationParameters& orthogonalizationParameters);

        /// @brief Freezes a line in the curvilinear orthogonalization process
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] firstGridLineNode The geometry list containing the first point of the line to freeze
        /// @param[in] secondGridLineNode The geometry list containing the second point of the line to freeze
        /// @returns  Error code
        MKERNEL_API int mkernel_set_frozen_lines_orthogonalize_curvilinear(int meshKernelId,
                                                                           const GeometryList& firstGridLineNode,
                                                                           const GeometryList& secondGridLineNode);

        /// @brief Define a block on the curvilinear grid where to perform orthogonalization
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] lowerLeftCorner  The geometry list containing the lower left corner of the block to orthogonalize
        /// @param[in] upperRightCorner The geometry list containing the upper left corner of the block to orthogonalize
        /// @return  Error code
        MKERNEL_API int mkernel_set_block_orthogonalize_curvilinear(int meshKernelId,
                                                                    const GeometryList& lowerLeftCorner,
                                                                    const GeometryList& upperRightCorner);

        /// @brief Orthogonalize a curvilinear grid
        /// @param[in] meshKernelId       The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_orthogonalize_curvilinear(int meshKernelId);

        /// @brief Resets the CurvilinearGridOrthogonalization instance in MeshKernelState
        /// @param[in] meshKernelId The id of the mesh state
        /// @return  Error code
        MKERNEL_API int mkernel_finalize_orthogonalize_curvilinear(int meshKernelId);

        /// @brief Smooths a curvilinear grid
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] smoothingIterations The number of smoothing iterations to perform
        /// @param[in] lowerLeftCorner The geometry list containing the lower left corner of the block to smooth
        /// @param[in] upperRightCorner The geometry list containing the upper right corner of the block to smooth
        /// @return Error code
        MKERNEL_API int mkernel_smoothing_curvilinear(int meshKernelId,
                                                      int smoothingIterations,
                                                      const GeometryList& lowerLeftCorner,
                                                      const GeometryList& upperRightCorner);

        /// @brief Smooths a curvilinear grid along the direction specified by a segment
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] smoothingIterations The number of smoothing iterations to perform
        /// @param[in] firstGridlineNode The first point of the segment
        /// @param[in] secondGridLineNode The second point of the segment
        /// @param[in] lowerLeftCornerSmoothingArea The geometry list containing the lower left corner of the smoothing area
        /// @param[in] upperRightCornerSmootingArea The geometry list containing the upper right corner of the smoothing area
        /// @return Error code
        MKERNEL_API int mkernel_smoothing_directional_curvilinear(int meshKernelId,
                                                                  int smoothingIterations,
                                                                  GeometryList const& firstGridlineNode,
                                                                  GeometryList const& secondGridLineNode,
                                                                  GeometryList const& lowerLeftCornerSmoothingArea,
                                                                  GeometryList const& upperRightCornerSmootingArea);

        /// @brief Instantiates the curvilinear line shift algorithm
        /// @param[in] meshKernelId The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_initialize_line_shift_curvilinear(int meshKernelId);

        /// @brief Sets the start and end nodes of the line to shift
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] firstGridLineNode The geometry list containing the first point of the line to shift
        /// @param[in] secondGridLineNode The geometry list containing the second point of the line to shift
        /// @returns  Error code
        MKERNEL_API int mkernel_set_line_line_shift_curvilinear(int meshKernelId,
                                                                GeometryList const& firstGridLineNode,
                                                                GeometryList const& secondGridLineNode);

        /// @brief Defines a block on the curvilinear where the shifting is distributed
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] lowerLeftCorner  The geometry list containing the lower left corner of the block
        /// @param[in] upperRightCorner The geometry list containing the upper right corner of the block
        /// @return  Error code
        MKERNEL_API int mkernel_set_block_line_shift_curvilinear(int meshKernelId,
                                                                 GeometryList const& lowerLeftCorner,
                                                                 GeometryList const& upperRightCorner);

        /// @brief Moves a node of the line to shift
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] fromCoordinate The geometry list containing the Cartesian coordinates of the node to move (the closest curvilinear grid node will be found)
        /// @param[in] toCoordinate  The geometry list containing the Cartesian coordinates of the new node position
        /// @return  Error code
        MKERNEL_API int mkernel_move_node_line_shift_curvilinear(int meshKernelId,
                                                                 GeometryList const& fromCoordinate,
                                                                 GeometryList const& toCoordinate);

        /// @brief Computes the new grid, shifting the line towards the moved nodes and distributing the shifting in block specified before
        /// @param[in] meshKernelId The id of the mesh state
        /// @return  Error code
        MKERNEL_API int mkernel_line_shift_curvilinear(int meshKernelId);

        /// @brief Resets the instance of the line shift algorithm in MeshKernelState
        /// @param[in] meshKernelId The id of the mesh state
        /// @return  Error code
        MKERNEL_API int mkernel_finalize_line_shift_curvilinear(int meshKernelId);

        /// @brief Inserts a new face on a curvilinear grid. The new face will be inserted on top of the closest edge by linear extrapolation.
        /// @param[in] meshKernelId       The id of the mesh state
        /// @param[in] point  The geometry list containing the point used for finding the closest edge.
        /// @returns Error code
        MKERNEL_API int mkernel_insert_face_curvilinear(int meshKernelId, const GeometryList& point);

        /// @brief Converts a curvilinear grid to an unstructured mesh
        MKERNEL_API int mkernel_convert_curvilinear_to_mesh2d(int meshKernelId);

        /// @brief Gets the double value used in the back-end library as separator and missing value
        /// @return The double missing value used in mesh kernel
        MKERNEL_API double mkernel_get_separator();

        /// @brief Gets the double value used to separate the inner part of a polygon from its outer part
        /// @return The double missing value used in mesh kernel
        MKERNEL_API double mkernel_get_inner_outer_separator();

        /// @brief Triangle interpolation (ec_module)
        /// @param[in]  mesh2d             The Mesh2D data
        /// @param[in]  samplesXCoordinate The sample x-coordinates
        /// @param[in]  samplesYCoordinate The sample y-coordinates
        /// @param[in]  samplesValue       The sample values
        /// @param[in]  numSamples         The number of samples
        /// @param[in]  locationType       The location type
        /// @param[in]  spherical          Current projection (0 cartesian, 1 spherical)
        /// @param[in]  sphericalAccurate  Accurate spherical projection (0 default spherical, 1 spherical accurate)
        /// @param[out] results            The interpolation results
        /// @return Error code
        MKERNEL_API int triangulation(const Mesh2D& mesh2d,
                                      const double** samplesXCoordinate,
                                      const double** samplesYCoordinate,
                                      const double** samplesValue,
                                      const int& numSamples,
                                      const int& locationType,
                                      const int& spherical,
                                      const int& sphericalAccurate,
                                      double** results);

        /// @brief AveragingInterpolation interpolation (ec_module)
        /// @param[in] meshKernelId           The id of the mesh state
        /// @param[in] samples                The samples coordinates and values
        /// @param[in] locationType           The location type (see MeshLocations enum)
        /// @param[in] averagingMethodType    The averaging method (see Method enum)
        /// @param[in] relativeSearchSize     The relative search size around the location (larger increases the number of samples considered)
        /// @param[in] results                The interpolation results with x and y coordinates
        /// @return Error code
        MKERNEL_API int mkernel_averaging_interpolation_mesh2d(int meshKernelId,
                                                               const GeometryList& samples,
                                                               const int& locationType,
                                                               const int& averagingMethodType,
                                                               const double& relativeSearchSize,
                                                               GeometryList& results);

        /// @brief Gets pointer to error message.
        /// @param[out] error_message
        /// @returns Error code
        MKERNEL_API int mkernel_get_error(const char*& error_message);

        /// @brief Gets the index of the erroneous entity.
        /// @param[out] invalidIndex The index of the erroneous entity
        /// @param[out] type         The entity type (node, edge or face, see MeshLocations)
        /// @returns Error code
        MKERNEL_API int mkernel_get_geometry_error(int& invalidIndex, int& type);

#ifdef __cplusplus
    }
#endif
} // namespace meshkernelapi
