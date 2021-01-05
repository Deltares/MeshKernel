//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2020.
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

#include <MeshKernelApi/CurvilinearParameters.hpp>
#include <MeshKernelApi/GeometryList.hpp>
#include <MeshKernelApi/InterpolationParameters.hpp>
#include <MeshKernelApi/MakeMeshParameters.hpp>
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

        /// @param[out] meshKernelId Identifier for the created grid state
        /// @returns Error code
        MKERNEL_API int mkernel_new_mesh(int& meshKernelId);

        /// @brief Deallocates mesh state
        /// @param[in] meshKernelId Id of the grid state
        /// @returns Error code
        MKERNEL_API int mkernel_deallocate_state(int meshKernelId);

        /// @brief Deletes a mesh in a polygon using several options
        /// @param[in] meshKernelId Id of the grid state
        /// @param[in] disposableGeometryList The polygon where to perform the operation
        /// @param[in] deletionOption The deletion option (to be detailed)
        /// @param[in] invertDeletion Inverts the deletion of selected features
        /// @returns Error code
        MKERNEL_API int mkernel_delete_mesh(int meshKernelId, const GeometryList& disposableGeometryList, int deletionOption, bool invertDeletion);

        /// @brief Sets the grid state
        /// @param[in] meshKernelId Id of the grid state
        /// @param[in] meshGeometryDimensions Mesh dimensions
        /// @param[in] meshGeometry Mesh data
        /// @param[in] isGeographic Cartesian (false) or spherical (true) mesh
        /// @returns Error code
        MKERNEL_API int mkernel_set_state(int meshKernelId, const Mesh2D& meshGeometry, bool isGeographic);

        /// @brief Gets the mesh state as a <see cref="Mesh2D"/> structure
        /// @param[in] meshKernelId Id of the grid state
        /// @param[out] meshGeometry Grid data
        /// @returns Error code
        MKERNEL_API int mkernel_get_mesh(int meshKernelId, Mesh2D& meshGeometry);

        /// @brief Gets the mesh faces
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[out] meshGeometry Mesh data (including face information)
        /// @returns Error code
        MKERNEL_API int mkernel_find_faces(int meshKernelId, Mesh2D& meshGeometry);

        /// @brief Count the number of hanging edges
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[out] numHangingEdgesIndices
        /// @returns Error code
        MKERNEL_API int mkernel_count_hanging_edges(int meshKernelId, int& numHangingEdgesIndices);

        /// @brief Gets the indices of hanging edges
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in,out] hangingEdgesIndices Pointer to memory where the hanging edge indices will be stored
        /// @returns Error code
        MKERNEL_API int mkernel_get_hanging_edges(int meshKernelId, int** hangingEdgesIndices);

        /// @brief Deletes the hanging edges
        /// @param[in] meshKernelId Id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_delete_hanging_edges(int meshKernelId);

        /// @brief Orthogonalization
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] projectToLandBoundaryOption The option to determine how to snap to land boundaries
        /// @param[in] orthogonalizationParameters The structure containing the orthogonalization parameters
        /// @param[in] geometryListPolygon The polygon where to perform the orthogonalization
        /// @param[in] geometryListLandBoundaries The land boundaries to account for in the orthogonalization process
        /// @returns Error code
        MKERNEL_API int mkernel_orthogonalize(int meshKernelId,
                                              int projectToLandBoundaryOption,
                                              const OrthogonalizationParameters& orthogonalizationParameters,
                                              const GeometryList& geometryListPolygon,
                                              const GeometryList& geometryListLandBoundaries);

        /// @brief Orthogonalization initialization (first function to use in interactive mode)
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] projectToLandBoundaryOption The option to determine how to snap to land boundaries
        /// @param[in] orthogonalizationParameters The structure containing the user defined orthogonalization parameters
        /// @param[in] geometryListPolygon The polygon where to perform the orthogonalization
        /// @param[in] geometryListLandBoundaries The land boundaries to account for in the orthogonalization process
        /// @returns Error code
        MKERNEL_API int mkernel_orthogonalize_initialize(int meshKernelId,
                                                         int projectToLandBoundaryOption,
                                                         OrthogonalizationParameters& orthogonalizationParameters,
                                                         const GeometryList& geometryListPolygon,
                                                         const GeometryList& geometryListLandBoundaries);

        /// @brief Prepares outer orthogonalization iteration (interactive mode)
        /// @param[in] meshKernelId Id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_orthogonalize_prepare_outer_iteration(int meshKernelId);

        /// @brief Performs inner orthogonalization iteration (interactive mode)
        /// @param[in] meshKernelId Id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_orthogonalize_inner_iteration(int meshKernelId);

        /// @brief Finalizes orthogonalization outer iteration (interactive mode)
        /// @param[in] meshKernelId
        /// @returns Error code
        MKERNEL_API int mkernel_orthogonalize_finalize_outer_iteration(int meshKernelId);

        /// @brief Cleans up back-end orthogonalization algorithm (interactive mode)
        /// @param[in] meshKernelId Id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_orthogonalize_delete(int meshKernelId);

        /// @brief Gets the orthogonality
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[out] geometryList The orthogonality values of each edge
        /// @returns Error code
        MKERNEL_API int mkernel_get_orthogonality(int meshKernelId, GeometryList& geometryList);

        /// @brief Gets the smoothness
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[out] geometryList The smoothness values of each edge
        /// @returns Error code
        MKERNEL_API int mkernel_get_smoothness(int meshKernelId, GeometryList& geometryList);

        /// @brief Get spline intermediate points
        /// @param[in] geometryListIn The input corner nodes of the splines
        /// @param[out] geometryListOut The output spline
        /// @param[out] numberOfPointsBetweenNodes The number of spline nodes between the corners points
        /// @returns Error code
        MKERNEL_API int mkernel_get_splines(const GeometryList& geometryListIn, GeometryList& geometryListOut, int numberOfPointsBetweenNodes);

        /// @brief Gets the coordinates of the closest existing node
        /// @param[in] meshKernelId Id of the grid state
        /// @param[in] geometryListIn Node coordinates
        /// @param[in] searchRadius the radius where to search for the node
        /// @param[out] geometryListOut Mesh node coordinates
        /// @returns Error code
        MKERNEL_API int mkernel_get_node_coordinate(int meshKernelId, const GeometryList& geometryListIn, double searchRadius, GeometryList& geometryListOut);

        /// @brief Makes curvilinear grid from splines with an advancing front.
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryList The input splines corners
        /// @param[in] curvilinearParameters The input parameters to generate the curvilinear grid
        /// @param[in] splinesToCurvilinearParameters The parameters of the advancing front algorithm
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_mesh_from_splines_ortho(int meshKernelId, const GeometryList& geometryList, const CurvilinearParameters& curvilinearParameters, const SplinesToCurvilinearParameters& splinesToCurvilinearParameters);

        /// @brief Generates a curvilinear grid from splines with the advancing front method. Initialization step (interactive)
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryList The input splines corners
        /// @param[in] curvilinearParameters The input parameters to generate the curvilinear grid
        /// @param[in] splinesToCurvilinearParameters The parameters of the advancing front algorithm
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_mesh_from_splines_ortho_initialize(int meshKernelId, const GeometryList& geometryList, const CurvilinearParameters& curvilinearParameters, const SplinesToCurvilinearParameters& splinesToCurvilinearParameters);

        /// @brief One advancment of the front in curvilinear grid from splines (interactive)
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] layer The layer index
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_mesh_from_splines_ortho_iteration(int meshKernelId, int layer);

        /// @brief Converts curvilinear grid to mesh and refreshes the state (interactive)
        /// @param[in] meshKernelId
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_mesh_from_splines_ortho_refresh_mesh(int meshKernelId);

        /// @brief Finalizes curvilinear grid from splines algorithm
        /// @param[in] meshKernelId
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_mesh_from_splines_ortho_delete(int meshKernelId);

        /// @brief Makes a new mesh
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] makeGridParameters The structure containing the make grid parameters
        /// @param[in] geometryList The polygon to account for
        /// @returns Error code
        MKERNEL_API int mkernel_make_mesh(int meshKernelId, const MakeMeshParameters& makeGridParameters, const GeometryList& geometryList);

        /// @brief Makes a triangular grid in a polygon
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryList The polygon where to triangulate
        /// @returns Error code
        MKERNEL_API int mkernel_make_mesh_from_polygon(int meshKernelId, const GeometryList& geometryList);

        /// @brief Makes a triangular mesh from samples
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryList The samples where to triangulate
        /// @returns Error code
        MKERNEL_API int mkernel_make_mesh_from_samples(int meshKernelId, GeometryList& geometryList);

        /// @brief Retrieves the mesh boundary polygon
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[out] geometryList The output network boundary polygon
        /// @returns Error code
        MKERNEL_API int mkernel_copy_mesh_boundaries_to_polygon(int meshKernelId, GeometryList& geometryList);

        /// @brief Counts the number of polygon nodes contained in the mesh boundary polygon
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[out] numberOfPolygonNodes The number of polygon points
        /// @returns Error code
        MKERNEL_API int mkernel_copy_mesh_boundaries_to_polygon_count_nodes(int meshKernelId, int& numberOfPolygonNodes);

        /// @brief Gets the refined polygon
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListIn The input polygons
        /// @param[in] firstIndex The index of the first node
        /// @param[in] secondIndex The index of the second node
        /// @param[in] distance The refinement distance
        /// @param[out] geometryListOut
        /// @returns Error code
        MKERNEL_API int mkernel_refine_polygon(int meshKernelId, const GeometryList& geometryListIn, int firstIndex, int secondIndex, double distance, GeometryList& geometryListOut);

        /// @brief Counts the number of nodes after polygon refinement
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListIn The input polygon
        /// @param[in] firstIndex The index of the first node
        /// @param[in] secondIndex The index of the second node
        /// @param[in] distance The refinement distance
        /// @param[out] numberOfPolygonNodes The number of nodes after refinement
        /// @returns Error code
        MKERNEL_API int mkernel_refine_polygon_count(int meshKernelId, GeometryList& geometryListIn, int firstIndex, int secondIndex, double distance, int& numberOfPolygonNodes);

        /// @brief Merges nodes within a distance of 0.001 m, effectively removing small edges
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListIn The polygon where to perform the operation
        /// @returns Error code
        MKERNEL_API int mkernel_merge_nodes(int meshKernelId, const GeometryList& geometryListIn);

        /// @brief Merges node \p startNode to \p endNode
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] startNode The index of the first node to merge
        /// @param[in] endNode The index of the second node to merge
        /// @returns Error code
        MKERNEL_API int mkernel_merge_two_nodes(int meshKernelId, int startNode, int endNode);

        /// @brief Gets the selected mesh node indices
        ///        \see how to pass arrays in https://www.mono-project.com/docs/advanced/pinvoke/#memory-management
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListIn The input polygons
        /// @param[in] inside Count nodes indices inside (1) or outside (0) the polygon
        /// @param[out] selectedNodes The selected nodes indices
        /// @returns Error code
        MKERNEL_API int mkernel_nodes_in_polygons(int meshKernelId, const GeometryList& geometryListIn, int inside, int** selectedNodes);

        /// @brief Counts the number of selected mesh node indices
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListIn The input polygons
        /// @param[in] inside Count nodes inside (1) or outside (0) the polygon
        /// @param[out] numberOfMeshNodes The number of selected nodes
        /// @returns Error code
        MKERNEL_API int mkernel_count_nodes_in_polygons(int meshKernelId, const GeometryList& geometryListIn, int inside, int& numberOfMeshNodes);

        /// @brief Insert a new edge connecting \p startNode and \p endNode
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] startNode The index of the first node to connect
        /// @param[in] endNode The index of the second node to connect
        /// @param[out] newEdgeIndex The index of the new edge
        /// @returns Error code
        MKERNEL_API int mkernel_insert_edge(int meshKernelId, int startNode, int endNode, int& newEdgeIndex);

        /// @brief Inserts a new node
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] xCoordinate X-coordinate of the new node
        /// @param[in] yCoordinate y-coordinate of the new node
        /// @param[out] nodeIndex The index of the new mesh node
        /// @returns Error code
        MKERNEL_API int mkernel_insert_node(int meshKernelId, double xCoordinate, double yCoordinate, int& nodeIndex);

        /// @brief Deletes a node with specified \p nodeIndex
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] nodeIndex The nodeIndex to delete
        /// @returns Error code
        MKERNEL_API int mkernel_delete_node(int meshKernelId, int nodeIndex);

        /// @brief Moves a selected node to a new position
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListIn The new coordinate
        /// @param[in] nodeIndex The node index (to be detailed)
        /// @returns Error code
        MKERNEL_API int mkernel_move_node(int meshKernelId, const GeometryList& geometryListIn, int nodeIndex);

        /// @brief Deletes the closest mesh edge within the search radius from the input point
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListIn The input point coordinates
        /// @returns Error code
        MKERNEL_API int mkernel_delete_edge(int meshKernelId, const GeometryList& geometryListIn);

        /// @brief Deletes the closest mesh edge within the search radius from the input point
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListIn The input point coordinates
        /// @param[out] edgeIndex The edge index
        /// @returns Error code
        MKERNEL_API int mkernel_find_edge(int meshKernelId, const GeometryList& geometryListIn, int& edgeIndex);

        /// @brief Offsets a polygon
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListIn The polygon to be offsetted
        /// @param[in] innerPolygon Compute inner (true) or outer (false) polygon
        /// @param[in] distance The offset distance
        /// @param[out] geometryListOut The offsetted polygon
        /// @returns Error code
        MKERNEL_API int mkernel_offsetted_polygon(int meshKernelId, const GeometryList& geometryListIn, bool innerPolygon, double distance, GeometryList& geometryListOut);

        /// @brief Gets the number of nodes of the offsetted polygon  Count the number of nodes after polygon refinement
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListIn The polygon to be offsetted
        /// @param[in] innerPolygon Compute inner (true) or outer (false) polygon
        /// @param[in] distance The offset distance
        /// @param[out] numberOfPolygonNodes The number of nodes of the generated polygon
        /// @returns Error code
        MKERNEL_API int mkernel_offsetted_polygon_count(int meshKernelId, const GeometryList& geometryListIn, bool innerPolygon, double distance, int& numberOfPolygonNodes);

        /// @brief Refines a grid based on the samples contained in the geometry list
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryList The sample set
        /// @param[in] interpolationParameters The interpolation parameters
        /// @param[in] sampleRefineParameters The interpolation settings related to the samples
        /// @returns Error code
        MKERNEL_API int mkernel_refine_mesh_based_on_samples(int meshKernelId,
                                                             const GeometryList& geometryList,
                                                             const InterpolationParameters& interpolationParameters,
                                                             const SampleRefineParameters& sampleRefineParameters);

        /// @brief Refines a grid based on polygon
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryList The closed polygon where to perform the refinement
        /// @param[in] interpolationParameters The interpolation parameters
        /// @returns Error code
        MKERNEL_API int mkernel_refine_mesh_based_on_polygon(int meshKernelId, const GeometryList& geometryList, const InterpolationParameters& interpolationParameters);

        /// @brief Finds the node index closest to the input point
        /// @param[in] meshKernelId
        /// @param[in] geometryListIn
        /// @param[in] searchRadius
        /// @param[out] nodeIndex
        /// @returns Error code
        MKERNEL_API int mkernel_get_node_index(int meshKernelId, const GeometryList& geometryListIn, double searchRadius, int& nodeIndex);

        /// @brief Selects points in polygons
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] inputPolygon The polygon(s) used for selection
        /// @param[in] inputPoints The points to select
        /// @param[out] selectedPoints The selected points in the zCoordinates field (0.0 not selected, 1.0 selected)
        /// @returns Error code
        MKERNEL_API int mkernel_points_in_polygon(int meshKernelId, const GeometryList& inputPolygon, const GeometryList& inputPoints, GeometryList& selectedPoints);

        /// @brief Flips the edges
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] isTriangulationRequired The option to triangulate also non triangular cells (if activated squares becomes triangles)
        /// @param[in] projectToLandBoundaryOption The option to determine how to snap to land boundaries
        /// @returns Error code
        MKERNEL_API int mkernel_flip_edges(int meshKernelId, int isTriangulationRequired, int projectToLandBoundaryOption);

        /// @brief Generates curvilinear grid from splines with transfinite interpolation
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListIn
        /// @param[in] curvilinearParameters
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_mesh_from_splines(int meshKernelId,
                                                              const GeometryList& geometryListIn,
                                                              const CurvilinearParameters& curvilinearParameters);

        /// @brief Computes a curvilinear mesh in a polygon. 3 separate polygon nodes need to be selected.
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] polygon The input polygon
        /// @param[in] firstNode The first selected node
        /// @param[in] secondNode The second selected node
        /// @param[in] thirdNode The third node
        /// @param[in] useFourthSide Use (true/false) the fourth polygon side to compute the curvilinear grid
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_from_polygon(int meshKernelId,
                                                         const GeometryList& polygon,
                                                         int firstNode,
                                                         int secondNode,
                                                         int thirdNode,
                                                         bool useFourthSide);

        /// @brief Computes a curvilinear mesh in a triangle. 3 separate polygon nodes need to be selected.
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] polygon The input polygons
        /// @param[in] firstNode The first selected node
        /// @param[in] secondNode The second selected node
        /// @param[in] thirdNode The third node
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_from_triangle(int meshKernelId,
                                                          const GeometryList& polygon,
                                                          int firstNode,
                                                          int secondNode,
                                                          int thirdNode);

        /// @brief Gets the number of obtuse triangles (those having one edge longer than the sum of the other two)
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[out] numObtuseTriangles The number of obtuse triangles
        /// @return
        MKERNEL_API int mkernel_get_obtuse_triangles_count(int meshKernelId, int& numObtuseTriangles);

        /// @brief Gets the obtuse triangle mass centers (those having one edge longer than the sum of the other two)
        /// @param[in] meshKernelId  Id of the mesh state
        /// @param[out] result The obtuse triangles mass centers
        /// @return Error code (0 Successful)
        MKERNEL_API int mkernel_get_obtuse_triangles(int meshKernelId, GeometryList& result);

        /// @brief Counts the small flow edges (flow edges are the edges connecting the face circumcenters)
        /// @param[in] meshKernelId  Id of the mesh state
        /// @param[in] smallFlowEdgesThreshold The configurable threshold for detecting the small flow edges
        /// @param[out] numSmallFlowEdges The number of the small flow edges
        /// @return
        MKERNEL_API int mkernel_get_small_flow_edge_centers_count(int meshKernelId, double smallFlowEdgesThreshold, int& numSmallFlowEdges);

        /// @brief Gets the small flow edges (flow edges are the edges connecting the face circumcenters)
        /// @param[in] meshKernelId  Id of the mesh state
        /// @param[in] smallFlowEdgesThreshold The configurable threshold for detecting the small flow edges
        /// @param[out] result The center points of the small flow edges
        /// @return Error code (0 Successful)
        MKERNEL_API int mkernel_get_small_flow_edge_centers(int meshKernelId, double smallFlowEdgesThreshold, GeometryList& result);

        /// @brief Deletes the small flow edges (flow edges are the edges connecting the face circumcenters)
        /// @param[in] meshKernelId  Id of the mesh state
        /// @param[in] smallFlowEdgesThreshold The configurable threshold for detecting the small flow edges
        /// @param[in] minFractionalAreaTriangles The ratio of the face area to the average area of neighboring non triangular faces
        /// @return Error code (0 Successful)
        MKERNEL_API int mkernel_delete_small_flow_edges(int meshKernelId, double smallFlowEdgesThreshold, double minFractionalAreaTriangles);

        /// @brief Gets the double value used in the back-end library as separator and missing value
        /// @return The double missing value used in mesh kernel
        MKERNEL_API double mkernel_get_separator();

        /// @brief Gets the double value used to separate the inner part of a polygon from its outer part
        /// @return The double missing value used in mesh kernel
        MKERNEL_API double mkernel_get_inner_outer_separator();

        /// @brief Triangle interpolation (ec_module)
        /// @param[in] meshGeometry Mesh data
        /// @param[in] startIndex start index (not used)
        /// @param[in] samplesXCoordinate The sample x coordinates
        /// @param[in] samplesYCoordinate The sample y coordinates
        /// @param[in] samplesValue The sample values
        /// @param[in] numSamples The number of samples
        /// @param[out] results The interpolation results
        /// @param[in] locationType The location type
        /// @param[in] spherical Current projection (0 cartesian, 1 spherical)
        /// @param[in] sphericalAccurate Accurate spherical projection (0 default spherical, 1 spherical accurate)
        /// @return Error code (0 Successful)
        MKERNEL_API int triangulation(const Mesh2D& meshGeometry,
                                      int& startIndex,
                                      const double** samplesXCoordinate,
                                      const double** samplesYCoordinate,
                                      const double** samplesValue,
                                      int& numSamples,
                                      double** results,
                                      int& locationType,
                                      int& spherical,
                                      int& sphericalAccurate);

        /// @brief AveragingInterpolation interpolation (ec_module)
        /// @param[in] meshGeometry Mesh data
        /// @param[in] startIndex Mesh data start index (not used)
        /// @param[in] samplesXCoordinate The sample x coordinates
        /// @param[in] samplesYCoordinate The sample y coordinates
        /// @param[in] samplesValue The sample values
        /// @param[in] numSamples The number of samples
        /// @param[out] results The interpolation results
        /// @param[in] locationType The location type (see MeshLocations enum)
        /// @param[in] Wu1Duni A setting for 1d meshes (not used)
        /// @param[in] averagingMethod The averaging method (see Method enum)
        /// @param[in] minNumberOfSamples The minimum amount of samples (not used)
        /// @param[in] relativeSearchSize The relative search size around the location (larger increases the number of samples considered)
        /// @param[in] spherical Current projection (0 cartesian, 1 spherical)
        /// @param[in] sphericalAccurate Accurate spherical computations (0 default spherical, 1 spherical accurate)
        /// @return Error code (0 Successful)
        MKERNEL_API int averaging(const Mesh2D& meshGeometry,
                                  const int& startIndex,
                                  const double** samplesXCoordinate,
                                  const double** samplesYCoordinate,
                                  const double** samplesValue,
                                  const int& numSamples,
                                  double** results,
                                  const int& locationType,
                                  const double& Wu1Duni,
                                  const int& averagingMethod,
                                  const int& minNumberOfSamples,
                                  const double& relativeSearchSize,
                                  const int& spherical,
                                  const int& sphericalAccurate);

        /// @brief Gets pointer to error message.
        /// @param[out] error_message
        /// @returns Error code
        MKERNEL_API int mkernel_get_error(const char*& error_message);

        /// @brief Gets the index of the erroneous entity.
        /// @param[out] invalidIndex The index of the erroneous entity
        /// @param[out] type The entity type (node, edge or face, see MeshLocations)
        /// @returns Error code
        MKERNEL_API int mkernel_get_geometry_error(int& invalidIndex, int& type);

#ifdef __cplusplus
    }
#endif
} // namespace meshkernelapi
