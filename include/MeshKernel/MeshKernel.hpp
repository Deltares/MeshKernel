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
#include <MeshKernel/MeshGeometryDimensions.hpp>
#include <MeshKernel/MeshGeometry.hpp>
#include <MeshKernel/GeometryListNative.hpp>
#include <MeshKernel/OrthogonalizationParametersNative.hpp>
#include <MeshKernel/CurvilinearParametersNative.hpp>
#include <MeshKernel/SplinesToCurvilinearParametersNative.hpp>
#include <MeshKernel/MakeGridParametersNative.hpp>
#include <MeshKernel/InterpolationParametersNative.hpp>
#include <MeshKernel/SampleRefineParametersNative.hpp>

#if defined(_WIN32)
#if !defined(MKERNEL_API)
#define MKERNEL_API __declspec(dllexport)
#endif
#else
#define MKERNEL_API __attribute__((visibility("default")))
#endif

// contains all mesh instances
namespace meshkernelapi
{
#ifdef __cplusplus
    extern "C"
    {
#endif
        /// @brief Create a new mesh state and return the generated <param name="meshKernelId"/>

        /// @param[out] meshKernelId Identifier for the created grid state
        /// @returns Error code
        MKERNEL_API int mkernel_new_mesh(int& meshKernelId);

        /// @brief Deallocate mesh state
        /// @param[in] meshKernelId Id of the grid state
        /// @returns Error code
        MKERNEL_API int mkernel_deallocate_state(int meshKernelId);

        /// @brief Deletes a mesh in a polygon using several options
        /// @param[in] meshKernelId Id of the grid state
        /// @param[in] disposableGeometryList The polygon where to perform the operation
        /// @param[in] deletionOption The deletion option (to be detailed)
        /// @param[in] invertDeletion Inverts the deletion of selected features
        /// @returns Error code
        MKERNEL_API int mkernel_delete_mesh(int meshKernelId, GeometryListNative& disposableGeometryList, int deletionOption, bool invertDeletion);

        /// @brief Set the grid state
        /// @param[in] meshKernelId Id of the grid state
        /// @param[in] meshGeometryDimensions Mesh dimensions
        /// @param[in] meshGeometry Mesh data
        /// @param[in] isGeographic Cartesian or spherical mesh
        /// @returns Error code
        MKERNEL_API int mkernel_set_state(int meshKernelId, const MeshGeometryDimensions& meshGeometryDimensions, const MeshGeometry& meshGeometry, bool isGeographic);

        /// @brief Gets the mesh state as a <see cref="MeshGeometry"/> structure
        /// @param[in] meshKernelId Id of the grid state
        /// @param[in,out] meshGeometryDimensions Mesh dimensions
        /// @param[in,out] meshGeometry Grid data
        /// @returns Error code
        MKERNEL_API int mkernel_get_mesh(int meshKernelId, MeshGeometryDimensions& meshGeometryDimensions, MeshGeometry& meshGeometry);

        /// @brief Gets the mesh faces
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in,out] meshGeometryDimensions Grid dimensions
        /// @param[in,out] meshGeometry Mesh data (including face information)
        /// @returns Error code
        MKERNEL_API int mkernel_find_faces(int meshKernelId, MeshGeometryDimensions& meshGeometryDimensions, MeshGeometry& meshGeometry);

        /// @brief Orthogonalization
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] isTriangulationRequired The option to triangulate also non triangular cells (if activated squares becomes triangles)
        /// @param[in] isAccountingForLandBoundariesRequired The option to account for land boundaries
        /// @param[in] projectToLandBoundaryOption The option to determine how to snap to land boundaries
        /// @param[in] orthogonalizationParametersNative The structure containing the orthogonalization parameters
        /// @param[in] geometryListNativePolygon The polygon where to perform the orthogonalization
        /// @param[in] geometryListNativeLandBoundaries The land boundaries to account for in the orthogonalization process
        /// @returns Error code
        MKERNEL_API int mkernel_orthogonalize(int meshKernelId,
                                              int isTriangulationRequired,
                                              int isAccountingForLandBoundariesRequired,
                                              int projectToLandBoundaryOption,
                                              const OrthogonalizationParametersNative& orthogonalizationParametersNative,
                                              const GeometryListNative& geometryListNativePolygon,
                                              const GeometryListNative& geometryListNativeLandBoundaries);

        /// @brief Orthogonalization initialization (first function to use in interactive mode)
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] isTriangulationRequired The option to triangulate also non triangular cells (if activated squares becomes triangles)
        /// @param[in] isAccountingForLandBoundariesRequired The option to account for land boundaries
        /// @param[in] projectToLandBoundaryOption The option to determine how to snap to land boundaries
        /// @param[in] orthogonalizationParameters The structure containing the user defined orthogonalization parameters
        /// @param[in] geometryListNativePolygon The polygon where to perform the orthogonalization
        /// @param[in] geometryListNativeLandBoundaries The land boundaries to account for in the orthogonalization process
        /// @returns Error code
        MKERNEL_API int mkernel_orthogonalize_initialize(int meshKernelId,
                                                         int isTriangulationRequired,
                                                         int isAccountingForLandBoundariesRequired,
                                                         int projectToLandBoundaryOption,
                                                         OrthogonalizationParametersNative& orthogonalizationParametersNative,
                                                         GeometryListNative& geometryListNativePolygon,
                                                         GeometryListNative& geometryListNativeLandBoundaries);

        /// @brief Prepare outer orthogonalization iteration (interactive mode)
        /// @param[in] meshKernelId Id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_orthogonalize_prepare_outer_iteration(int meshKernelId);

        /// @brief Perform inner orthogonalization iteration (interactive mode)
        /// @param[in] meshKernelId Id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_orthogonalize_inner_iteration(int meshKernelId);

        /// @brief Finalize orthogonalization outer iteration (interactive mode)
        /// @param[in] meshKernelId
        /// @returns Error code
        MKERNEL_API int mkernel_orthogonalize_finalize_outer_iteration(int meshKernelId);

        /// @brief Clean up back-end orthogonalization algorithm (interactive mode)
        /// @param[in] meshKernelId Id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_orthogonalize_delete(int meshKernelId);

        /// @brief Gets the orthogonality
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in,out] geometryListIn The orthogonality values of each edge
        /// @returns Error code
        MKERNEL_API int mkernel_get_orthogonality(int meshKernelId, GeometryListNative& geometryList);

        /// @brief Gets the smoothness
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in,out] geometryListIn The smoothness values of each edge
        /// @returns Error code
        MKERNEL_API int mkernel_get_smoothness(int meshKernelId, GeometryListNative& geometryList);

        /// @brief Get spline intermediate points
        /// @param[in] disposableGeometryListIn The input corner vertices of the splines
        /// @param[in,out] disposableGeometryListOut The output spline
        /// @param[in,out] numberOfPointsBetweenVertices The number of spline vertices between the corners points
        /// @returns Error code
        MKERNEL_API int mkernel_get_splines(const GeometryListNative& geometryListIn, GeometryListNative& geometry_list_out, int number_of_points_between_vertices);

        /// @brief Get the coordinates of the closest existing vertex
        /// @param[in] meshKernelId Id of the grid state
        /// @param[in] geometryListIn Vertex coordinates
        /// @param[in] searchRadius the radius where to search for the vertex
        /// @param[in,out] geometryListOut Mesh vertex coordinates
        /// @returns Error code
        MKERNEL_API int mkernel_get_node_coordinate(int meshKernelId, GeometryListNative& geometryListIn, double searchRadius, GeometryListNative& geometryListOut);

        /// @brief Make curvilinear grid from splines with an advancing front.
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListNative The input splines corners
        /// @param[in] curvilinearParametersNative The input parameters to generate the curvilinear grid
        /// @param[in] splinesToCurvilinearParametersNative The parameters of the advancing front algorithm
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_mesh_from_splines_ortho(int meshKernelId, const GeometryListNative& geometryListNative, const CurvilinearParametersNative& curvilinearParameters, const SplinesToCurvilinearParametersNative& splineToCurvilinearParameters);

        /// @brief Generate a curvilinear grid from splines with the advancing front method. Initialization step (interactive)
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListNative The input splines corners
        /// @param[in] curvilinearParametersNative The input parameters to generate the curvilinear grid
        /// @param[in] splinesToCurvilinearParametersNative The parameters of the advancing front algorithm
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_mesh_from_splines_ortho_initialize(int meshKernelId, const GeometryListNative& geometryListNative, const CurvilinearParametersNative& curvilinearParametersNative, const SplinesToCurvilinearParametersNative& splinesToCurvilinearParametersNative);

        /// @brief One advancment of the front in curvilinear grid from splines (interactive)
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] layer The layer index
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_mesh_from_splines_iteration(int meshKernelId, int layer);

        /// @brief Converts curvilinear grid to mesh and refreshes the state (interactive)
        /// @param[in] meshKernelId
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_mesh_from_splines_ortho_refresh_mesh(int meshKernelId);

        /// @brief Finalize curvilinear grid from splines algorithm
        /// @param[in] meshKernelId
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_mesh_from_splines_ortho_delete(int meshKernelId);

        /// @brief Make a new mesh
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] makeGridParameters The structure containing the make grid parameters
        /// @param[in] geometryListNative The polygon to account for
        /// @returns Error code
        MKERNEL_API int mkernel_make_mesh(int meshKernelId, const MakeGridParametersNative& makeGridParameters, const GeometryListNative& geometryListNative);

        /// @brief Make a triangular grid in a polygon
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListNative The polygon where to triangulate
        /// @returns Error code
        MKERNEL_API int mkernel_make_mesh_from_polygon(int meshKernelId, const GeometryListNative& geometryListNative);

        /// @brief Make a triangular mesh from samples
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListNative The samples where to triangulate
        /// @returns Error code
        MKERNEL_API int mkernel_make_mesh_from_samples(int meshKernelId, GeometryListNative& geometryListNative);

        /// @brief Retrieves the mesh boundary polygon
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in,out] geometryListNative The output network boundary polygon
        /// @returns Error code
        MKERNEL_API int mkernel_copy_mesh_boundaries_to_polygon(int meshKernelId, GeometryListNative& geometryListNative);

        /// @brief Counts the number of polygon vertices contained in the mesh boundary polygon
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in,out] numberOfPolygonVertices The number of polygon points
        /// @returns Error code
        MKERNEL_API int mkernel_copy_mesh_boundaries_to_polygon_count_vertices(int meshKernelId, int& numberOfPolygonVertices);

        /// @brief Gets the refined polygon
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListIn The input polygons
        /// @param[in] firstIndex The index of the first vertex
        /// @param[in] secondIndex The index of the second vertex
        /// @param[in] distance The refinement distance
        /// @param[in,out] geometryListOut
        /// @returns Error code
        MKERNEL_API int mkernel_refine_polygon(int meshKernelId, const GeometryListNative& geometryListIn, int firstIndex, int secondIndex, double distance, GeometryListNative& geometryListOut);

        /// @brief Count the number of vertices after polygon refinement
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListIn The input polygon
        /// @param[in] firstIndex The index of the first vertex
        /// @param[in] secondIndex The index of the second vertex
        /// @param[in] distance The refinement distance
        /// @param[in,out] numberOfPolygonVertices The number of vertices after refinement
        /// @returns Error code
        MKERNEL_API int mkernel_refine_polygon_count(int meshKernelId, GeometryListNative& geometryListIn, int firstIndex, int secondIndex, double distance, int& numberOfPolygonVertices);

        /// @brief Merges vertices within a distance of 0.001 m, effectively removing small edges
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListIn The polygon where to perform the operation
        /// @returns Error code
        MKERNEL_API int mkernel_merge_nodes(int meshKernelId, GeometryListNative& geometryListIn);

        /// @brief Merges vertex <param name="startVertexIndex"/> to <param name="endVertexIndex"/>
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] startNode The index of the first vertex to merge
        /// @param[in] endNode The index of the second vertex to merge
        /// @returns Error code
        MKERNEL_API int mkernel_merge_two_nodes(int meshKernelId, int startNode, int endNode);

        /// @brief Gets the selected mesh node indexes  (see how to pass arrays in https://www.mono-project.com/docs/advanced/pinvoke/#memory-management)
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListIn The input polygons
        /// @param[in] numberOfMeshVertices The number of selected nodes
        /// @param[in,out] selectedVerticesPtr The selected vertices indexes
        /// @returns Error code
        MKERNEL_API int mkernel_nodes_in_polygons(int meshKernelId, GeometryListNative& geometryListIn, int inside, int numberOfMeshVertices, int** selectedVertices);

        /// @brief Counts the number of selected mesh node indexes
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListIn The input polygons
        /// @param[in,out] numberOfMeshVertices The number of selected nodes
        /// @returns Error code
        MKERNEL_API int mkernel_count_nodes_in_polygons(int meshKernelId, GeometryListNative& geometryListIn, int inside, int& numberOfMeshVertices);

        /// @brief Insert a new edge connecting <param name="startVertexIndex"/> and <param name="endVertexIndex"/>
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] startNode The index of the first node to connect
        /// @param[in] endNode The index of the second node to connect
        /// @param[in,out] newEdgeIndex The index of the new edge
        /// @returns Error code
        MKERNEL_API int mkernel_insert_edge(int meshKernelId, int startNode, int endNode, int& newEdgeIndex);

        /// @brief Inserts a new node
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] disposableGeometryList The polygon where to perform the operation
        /// @param[in,out] vertexIndex The index of the new mesh node
        /// @returns Error code
        MKERNEL_API int mkernel_insert_node(int meshKernelId, double xCoordinate, double yCoordinate, double zCoordinate, int& vertexIndex);

        /// @brief Deletes a node with specified <param name="nodeIndex"/>
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] nodeIndex The nodeIndex to delete
        /// @returns Error code
        MKERNEL_API int mkernel_delete_node(int meshKernelId, int nodeIndex);

        /// @brief Function to move a selected node to a new position
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListIn The new coordinate
        /// @param[in] nodeIndex The node index (to be detailed)
        /// @returns Error code
        MKERNEL_API int mkernel_move_node(int meshKernelId, GeometryListNative& geometryListIn, int nodeIndex);

        /// @brief Deletes the closest mesh edge within the search radius from the input point
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListIn The input point coordinates
        /// @param[in] searchRadius The search radius
        /// @returns Error code
        MKERNEL_API int mkernel_delete_edge(int meshKernelId, GeometryListNative& geometryListIn);

        /// @brief Deletes the closest mesh edge within the search radius from the input point
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListIn The input point coordinates
        /// @param[in] searchRadius The search radius
        /// @param[in,out] edgeIndex The edge index
        /// @returns Error code
        MKERNEL_API int mkernel_find_edge(int meshKernelId, GeometryListNative& geometryListIn, int& edgeIndex);

        /// @brief Offset a polygon
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListIn The coordinate of the offset point
        /// @param[in] innerPolygon Compute inner/outer polygon
        /// @param[in] distance The offset distance
        /// @param[in,out] geometryListOut The offsetted polygon
        /// @returns Error code
        MKERNEL_API int mkernel_offsetted_polygon(int meshKernelId, GeometryListNative& geometryListIn, bool innerPolygon, double distance, GeometryListNative& geometryListOut);

        /// @brief Get the number of vertices of the offsetted polygon  Count the number of vertices after polygon refinement
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListIn The coordinate of the offset point
        /// @param[in] innerPolygon Compute inner/outer polygon
        /// @param[in] distance The offset distance
        /// @param[in,out] numberOfPolygonVertices The number of vertices of the generated polygon
        /// @returns Error code
        MKERNEL_API int mkernel_offsetted_polygon_count(int meshKernelId, GeometryListNative& geometryListIn, bool innerPolygon, double distance, int& numberOfPolygonVertices);

        /// @brief Refine a grid based on the samples contained in the geometry list
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListNative The sample set
        /// @param[in] interpolationParametersNative The interpolation parameters
        /// @param[in] sampleRefineParametersNative The interpolation settings related to the samples
        /// @returns Error code
        MKERNEL_API int mkernel_refine_mesh_based_on_samples(int meshKernelId, GeometryListNative& geometryListNative, InterpolationParametersNative& interpolationParametersNative, SampleRefineParametersNative& sampleRefineParametersNative);

        /// @brief Refine a grid based on polygon
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListNative The closed polygon where to perform the refinement
        /// @param[in] interpolationParametersNative The interpolation parameters
        /// @returns Error code
        MKERNEL_API int mkernel_refine_mesh_based_on_polygon(int meshKernelId, GeometryListNative& geometryListNative, InterpolationParametersNative& interpolationParametersNative);

        /// @brief Finds the vertex index closest to the input point
        /// @param[in] meshKernelId
        /// @param[in] geometryListIn
        /// @param[in] searchRadius
        /// @param[in,out] vertexIndex
        /// @returns Error code
        MKERNEL_API int mkernel_get_node_index(int meshKernelId, GeometryListNative& geometryListIn, double searchRadius, int& vertexIndex);

        /// @brief Selects points in polygons
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] inputPolygon The polygon(s) used for selection
        /// @param[in] inputPoints The points to select
        /// @param[in,out] selectedPoints The selected points in the zCoordinates field (0.0 not selected, 1.0 selected)
        /// @returns Error code
        MKERNEL_API int mkernel_points_in_polygon(int meshKernelId, GeometryListNative& inputPolygon, GeometryListNative& inputPoints, GeometryListNative& selectedPoints);

        /// @brief Flip the edges
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] isTriangulationRequired The option to triangulate also non triangular cells (if activated squares becomes triangles)
        /// @param[in] isAccountingForLandBoundariesRequired The option to account for land boundaries
        /// @param[in] projectToLandBoundaryOption The option to determine how to snap to land boundaries
        /// @returns Error code
        MKERNEL_API int mkernel_flip_edges(int meshKernelId, int isTriangulationRequired, int isAccountingForLandBoundariesRequired, int projectToLandBoundaryOption);

        /// @brief Generates curvilinear grid from splines with transfinite interpolation
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] geometryListNativeIn
        /// @param[in] curvilinearParametersNative
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_mesh_from_splines(int meshKernelId, GeometryListNative& geometryListNativeIn, CurvilinearParametersNative& curvilinearParametersNative);

        /// @brief Computes a curvilinear mesh in a polygon. 3 separate polygon nodes need to be selected.
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] polygonNative The input polygon
        /// @param[in] firstNode The first selected node
        /// @param[in] secondNode The second selected node
        /// @param[in] thirdNode The third node
        /// @param[in] useFourthSide Use (true/false) the fourth polygon side to compute the curvilinear grid
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_from_polygon(int meshKernelId, GeometryListNative& polygonNative, int firstNode, int secondNode, int thirdNode, bool useFourthSide);

        /// @brief Computes a curvilinear mesh in a triangle. 3 separate polygon nodes need to be selected.
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in] polygonNative The input polygons
        /// @param[in] firstNode The first selected node
        /// @param[in] secondNode The second selected node
        /// @param[in] thirdNode The third node
        /// @returns Error code
        MKERNEL_API int mkernel_curvilinear_from_triangle(int meshKernelId, GeometryListNative& polygonNative, int firstNode, int secondNode, int thirdNode);

        /// @brief Gets the number of obtuse triangles (those having one edge longer than the sum of the other two)
        /// @param[in] meshKernelId Id of the mesh state
        /// @param[in,out] numObtuseTriangles The number of obtuse triangles
        /// @return
        MKERNEL_API int mkernel_get_obtuse_triangles_count(int meshKernelId, int& numObtuseTriangles);

        /// @brief Gets the obtuse triangle mass centers (those having one edge longer than the sum of the other two)
        /// @param[in] meshKernelId  Id of the mesh state
        /// @param[in,out] result The obtuse triangles mass centers
        /// @return Error code (0 Successful)
        MKERNEL_API int mkernel_get_obtuse_triangles(int meshKernelId, GeometryListNative& result);

        /// @brief Count the small flow edges (flow edges are the edges connecting the face circumcenters)
        /// @param[in] meshKernelId  Id of the mesh state
        /// @param[in] smallFlowEdgesThreshold The configurable threshold for detecting the small flow edges
        /// @param[out] numSmallFlowEdges The number of the small flow edges
        /// @return
        MKERNEL_API int mkernel_get_small_flow_edge_centers_count(int meshKernelId, double smallFlowEdgesThreshold, int& numSmallFlowEdges);

        /// @brief Gets the small flow edges (flow edges are the edges connecting the face circumcenters)
        /// @param[in] meshKernelId  Id of the mesh state
        /// @param[in] smallFlowEdgesThreshold The configurable threshold for detecting the small flow edges
        /// @param[in,out] result The center points of the small flow edges
        /// @return Error code (0 Successful)
        MKERNEL_API int mkernel_get_small_flow_edge_centers(int meshKernelId, double smallFlowEdgesThreshold, GeometryListNative& result);

        /// @brief Removes the small flow edges (flow edges are the edges connecting the face circumcenters)
        /// @param[in] meshKernelId  Id of the mesh state
        /// @param[in] smallFlowEdgesThreshold The configurable threshold for detecting the small flow edges
        /// @param[in] minFractionalAreaTriangles The ratio of the face area to the average area of neighboring non triangular faces
        /// @return Error code (0 Successful)
        MKERNEL_API int mkernel_remove_small_flow_edges(int meshKernelId, double smallFlowEdgesThreshold, double minFractionalAreaTriangles);

        /// @brief Triangle interpolation (ec_module)
        /// @param[in] meshGeometryDimensions Mesh dimensions
        /// @param[in] meshGeometry Mesh data
        /// @param[in] startIndex start index (not used)
        /// @param[in] samplesXCoordinate The sample x coordinates
        /// @param[in] samplesYCoordinate The sample y coordinates
        /// @param[in] samplesValue The sample values
        /// @param[in] numSamples The number of samples
        /// @param[in,out] results The interpolation results
        /// @param[in] locationType The location type (see \ref InterpolationLocation)
        /// @param[in] spherical Current projection (0 cartesian, 1 spherical)
        /// @param[in] sphericalAccurate Accurate spherical projections (0 default spherical, 1 spherical accurate)
        /// @return Error code (0 Successful)
        MKERNEL_API int triangulation(const MeshGeometryDimensions& meshGeometryDimensions,
                                      const MeshGeometry& meshGeometry,
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
        /// @param[in] meshGeometryDimensions Mesh dimensions
        /// @param[in] meshGeometry Mesh data
        /// @param[in] startIndex Mesh data start index (not used)
        /// @param[in] samplesXCoordinate The sample x coordinates
        /// @param[in] samplesYCoordinate The sample y coordinates
        /// @param[in] samplesValue The sample values
        /// @param[in] numSamples The number of samples
        /// @param[in,out] results The interpolation results
        /// @param[in] locationType The location type (see InterpolationLocation enum)
        /// @param[in] Wu1Duni A setting for 1d meshes (not used)
        /// @param[in] averagingMethod The averaging method (see Method enum)
        /// @param[in] minNumberOfSamples The minimum amount of samples (not used)
        /// @param[in] relativeSearchSize The relative search size around the location (larger increases the number of samples considered)
        /// @param[in] spherical Current projection (0 cartesian, 1 spherical)
        /// @param[in] sphericalAccurate Accurate spherical computations (0 default spherical, 1 spherical accurate)
        /// @return Error code (0 Successful)
        MKERNEL_API int averaging(const MeshGeometryDimensions& meshGeometryDimensions,
                                  const MeshGeometry& meshGeometry,
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

        /// @brief Get pointer to error message.
        /// @param[out] error_message
        /// @returns Error code
        MKERNEL_API int mkernel_get_error(const char*& error_message);

#ifdef __cplusplus
    }
#endif
} // namespace meshkernelapi
