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
#include "MeshGeometryDimensions.hpp"
#include "MeshGeometry.hpp"
#include "GeometryListNative.hpp"
#include "OrthogonalizationParametersNative.hpp"
#include "CurvilinearParametersNative.hpp"
#include "SplinesToCurvilinearParametersNative.hpp"
#include "MakeGridParametersNative.hpp"
#include "InterpolationParametersNative.hpp"
#include "SampleRefineParametersNative.hpp"

#if defined(_WIN32) 
#if !defined(GRIDGEOM_API)
#define GRIDGEOM_API __declspec(dllexport)
#endif
#else  
#define GRIDGEOM_API __attribute__((visibility("default")))
#endif

// contains all mesh instances
namespace GridGeomApi
{
#ifdef __cplusplus
    extern "C"
    {
#endif
        /// <summary>
        /// Create a new mesh state and return the generated <param name="gridStateId"/>
        /// </summary>
        /// <param name="gridStateId">Identifier for the created grid state</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_new_mesh(int& gridStateId);

        /// <summary>
        /// Deallocate mesh state
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_deallocate_state(int gridStateId);

        /// <summary>
        /// Deletes a mesh in a polygon using several options
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="disposableGeometryList">The polygon where to perform the operation</param>
        /// <param name="deletionOption">The deletion option (to be detailed)</param>
        /// <param name="invertDeletion">Inverts the deletion of selected features</param>
        /// <returns>If the operation succeeded</returns>
        GRIDGEOM_API int ggeo_delete_mesh(int gridStateId, GeometryListNative& disposableGeometryList, int deletionOption, bool invertDeletion);

        /// <summary>
        /// Synchronize provided mesh (<param name="meshGeometryDimensions"/> and <param name="meshGeometry"/>) with the grid state with <param name="gridStateId"/>
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="meshGeometryDimensions">Mesh dimensions</param>
        /// <param name="meshGeometry">Mesh data</param>
        /// <param name="isGeographic">Cartesian or spherical mesh</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_set_state(int gridStateId, MeshGeometryDimensions& meshGeometryDimensions, MeshGeometry& meshGeometry, bool isGeographic);

        /// <summary>
        /// Gets the mesh state as a <see cref="MeshGeometry"/> structure
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="meshGeometryDimensions">Mesh dimensions</param>
        /// <param name="meshGeometry">Grid data</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_get_mesh(int gridStateId, MeshGeometryDimensions& meshGeometryDimensions, MeshGeometry& meshGeometry);

        /// <summary>
        /// Gets the mesh state as a <see cref="MeshGeometry"/> structure including faces information
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="meshGeometryDimensions">Grid dimensions</param>
        /// <param name="meshGeometry">Mesh data</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_find_faces(int gridStateId, MeshGeometryDimensions& meshGeometryDimensions, MeshGeometry& meshGeometry);

        /// <summary>
        /// Orthogonalization
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="isTriangulationRequired">The option to triangulate also non triangular cells (if activated squares becomes triangles) </param>
        /// <param name="isAccountingForLandBoundariesRequired">The option to account for land boundaries</param>
        /// <param name="projectToLandBoundaryOption">The option to determine how to snap to land boundaries</param>
        /// <param name="orthogonalizationParametersNative">The structure containing the orthogonalization parameters</param>
        /// <param name="geometryListNativePolygon">The polygon where to perform the orthogonalization</param>
        /// <param name="geometryListNativeLandBoundaries">The land boundaries to account for in the orthogonalization process</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_orthogonalize(int gridStateId, int isTriangulationRequired, int isAccountingForLandBoundariesRequired, int projectToLandBoundaryOption,
            OrthogonalizationParametersNative& orthogonalizationParametersNative, GeometryListNative& geometryListNativePolygon, GeometryListNative& geometryListNativeLandBoundaries);
       
        /// <summary>
        /// Orthogonalization initialization (first function to use in interactive mode)
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="isTriangulationRequired">The option to triangulate also non triangular cells (if activated squares becomes triangles) </param>
        /// <param name="isAccountingForLandBoundariesRequired">The option to account for land boundaries</param>
        /// <param name="projectToLandBoundaryOption">The option to determine how to snap to land boundaries</param>
        /// <param name="orthogonalizationParameters">The structure containing the user defined orthogonalization parameters</param>
        /// <param name="geometryListNativePolygon">The polygon where to perform the orthogonalization</param>
        /// <param name="geometryListNativeLandBoundaries">The land boundaries to account for in the orthogonalization process</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_orthogonalize_initialize(int gridStateId,
            int isTriangulationRequired, 
            int isAccountingForLandBoundariesRequired, 
            int projectToLandBoundaryOption,
            OrthogonalizationParametersNative& orthogonalizationParametersNative, 
            GeometryListNative& geometryListNativePolygon, 
            GeometryListNative& geometryListNativeLandBoundaries);

        /// <summary>
        /// Prepare outer orthogonalization iteration (interactive mode)
        /// </summary>
        /// <param name="gridGeomId">Id of the mesh state</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_orthogonalize_prepare_outer_iteration(int gridStateId);

        /// <summary>
        /// Perform inner orthogonalization iteration (interactive mode)
        /// </summary>
        /// <param name="gridGeomId">Id of the mesh state</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_orthogonalize_inner_iteration(int gridStateId);

        /// <summary>
        /// Finalize orthogonalization outer iteration (interactive mode)
        /// </summary>
        /// <param name="gridStateId"></param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_orthogonalize_finalize_outer_iteration(int gridStateId);
         
        /// <summary>
        /// Clean up back-end orthogonalization algorithm (interactive mode)
        /// </summary>
        /// <param name="gridGeomId">Id of the mesh state</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_orthogonalize_delete(int gridStateId);

        /// <summary>
        /// Gets the orthogonality
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="geometryListIn">The orthogonality values of each edge</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_get_orthogonality(int gridStateId, GeometryListNative& geometryList);

        /// <summary>
        /// Gets the smoothness 
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="geometryListIn">The smoothness values of each edge</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_get_smoothness(int gridStateId, GeometryListNative& geometryList);

        /// <summary>
        /// Get spline intermediate points 
        /// </summary>
        /// <param name="disposableGeometryListIn">The input corner vertices of the splines</param>
        /// <param name="disposableGeometryListOut">The output spline </param>
        /// <param name="numberOfPointsBetweenVertices">The number of spline vertices between the corners points</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_get_splines(GeometryListNative& geometryListIn, GeometryListNative& geometry_list_out, int number_of_points_between_vertices);

        /// <summary>
        /// Make curvilinear grid from splines with an advancing front.
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="geometryListNative">The input splines corners</param>
        /// <param name="curvilinearParametersNative">The input parameters to generate the curvilinear grid</param> 
        /// <param name="splinesToCurvilinearParametersNative">The parameters of the advancing front algorithm</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_curvilinear_mesh_from_splines_ortho(int gridStateId, GeometryListNative& geometryListNative, CurvilinearParametersNative& curvilinearParameters, SplinesToCurvilinearParametersNative& splineToCurvilinearParameters);

        /// <summary>
        /// Generate a curvilinear grid from splines with the advancing front method. Initialization step (interactive)
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="geometryListNative">The input splines corners</param>
        /// <param name="curvilinearParametersNative">The input parameters to generate the curvilinear grid</param>
        /// <param name="splinesToCurvilinearParametersNative">The parameters of the advancing front algorithm</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_curvilinear_mesh_from_splines_ortho_initialize(int gridStateId, GeometryListNative& geometryListNative, CurvilinearParametersNative& curvilinearParametersNative, SplinesToCurvilinearParametersNative& splinesToCurvilinearParametersNative);

        /// <summary>
        /// One advancment of the front in curvilinear grid from splines (interactive)
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="layer">The layer index</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_curvilinear_mesh_from_splines_iteration(int gridStateId, int layer);

        /// <summary>
        /// Converts curvilinear grid to mesh and refreshes the state (interactive)
        /// </summary>
        /// <param name="gridStateId"></param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_curvilinear_mesh_from_splines_ortho_refresh_mesh(int gridStateId);

        /// <summary>
        /// Finalize curvilinear grid from splines algorithm
        /// </summary>
        /// <param name="gridStateId"></param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_curvilinear_mesh_from_splines_ortho_delete(int gridStateId);

        /// <summary>
        /// Make a new mesh
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="makeGridParameters">The structure containing the make grid parameters </param>
        /// <param name="geometryListNative">The polygon to account for</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_make_mesh(int gridStateId, MakeGridParametersNative& makeGridParameters, GeometryListNative& geometryListNative);

        /// <summary>
        /// Make a triangular grid in a polygon
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="geometryListNative">The polygon where to triangulate</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_make_mesh_from_polygon(int gridStateId, GeometryListNative& geometryListNative);

        /// <summary>
        /// Make a triangular mesh from samples
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="geometryListNative">The samples where to triangulate</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_make_mesh_from_samples(int gridStateId, GeometryListNative& geometryListNative);

        /// <summary>
        /// Retrives the mesh boundary polygon
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="geometryListNative">The output network boundary polygon</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_copy_mesh_boundaries_to_polygon(int gridStateId, GeometryListNative& geometryListNative);

        /// <summary>
        /// Counts the number of polygon vertices contained in the mesh boundary polygon
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="numberOfPolygonVertices">The number of polygon points</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_copy_mesh_boundaries_to_polygon_count_vertices(int gridStateId, int& numberOfPolygonVertices);

        /// <summary>
        /// Gets the refined polygon
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="geometryListIn">The input polygons</param>
        /// <param name="firstIndex">The index of the first vertex</param>
        /// <param name="secondIndex">The index of the second vertex</param>
        /// <param name="distance">The refinement distance</param>
        /// <param name="geometryListOut"></param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_refine_polygon(int gridStateId, GeometryListNative& geometryListIn, int& firstIndex, int& secondIndex, double& distance, GeometryListNative& geometryListOut);

        /// <summary>
        /// Count the number of vertices after polygon refinment
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="geometryListIn">The input polygon</param>
        /// <param name="firstIndex">The index of the first vertex</param>
        /// <param name="secondIndex">The index of the second vertex</param>
        /// <param name="distance">The refinement distance</param>
        /// <param name="numberOfPolygonVertices">The number of vertices after refinement </param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_refine_polygon_count(int gridStateId, GeometryListNative& geometryListIn, int& firstIndex, int& secondIndex, double& distance, int& numberOfPolygonVertices);

        /// <summary>
        /// Merges vertices within a distance of 0.001 m, effectively removing small edges 
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="geometryListIn">The polygon where to perform the operation</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_merge_nodes(int gridStateId, GeometryListNative& geometryListIn);

        /// <summary>
        /// Merges vertex <param name="startVertexIndex"/> to <param name="endVertexIndex"/>
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="startNode">The index of the first vertex to merge</param>
        /// <param name="endNode">The index of the second vertex to merge</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_merge_two_nodes(int gridStateId, int startNode, int endNode);

        /// <summary>
        /// Gets the selected mesh node indexes  (see how to pass arrays in https://www.mono-project.com/docs/advanced/pinvoke/#memory-management)
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="geometryListIn">The input polygons</param>
        /// <param name="numberOfMeshVertices">The number of selected nodes</param>
        /// <param name="selectedVerticesPtr">The selected vertices indexes</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_nodes_in_polygons(int gridStateId, GeometryListNative& geometryListIn, int inside, int numberOfMeshVertices, int** selectedVertices);

        /// <summary>
        /// Counts the number of selected mesh node indexes
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="geometryListIn">The input polygons</param>
        /// <param name="numberOfMeshVertices">The number of selected nodes</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_count_nodes_in_polygons(int gridStateId, GeometryListNative& geometryListIn, int inside, int& numberOfMeshVertices);

        /// <summary>
        /// Insert a new edge connecting <param name="startVertexIndex"/> and <param name="endVertexIndex"/>
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="startNode">The index of the first node to connect</param>
        /// <param name="endNode">The index of the second node to connect</param>
        /// <param name="newEdgeIndex">The index of the new edge</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_insert_edge(int gridStateId, int startNode, int endNode, int& newEdgeIndex);

        /// <summary>
        /// Inserts a new node
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="disposableGeometryList">The polygon where to perform the operation</param>
        /// <param name="vertexIndex">The index of the new mesh node</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_insert_node(int gridStateId, double xCoordinate, double yCoordinate, double zCoordinate, int& vertexIndex);

        /// <summary>
        /// Deletes a node with specified <param name="nodeIndex"/>
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="nodeIndex">The nodeIndex to delete</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_delete_node(int gridStateId, int nodeIndex);

        /// <summary>
        /// Function to move a selected node to a new position
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="geometryListIn">The new coordinate</param>
        /// <param name="nodeIndex">The node index (to be detailed)</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_move_node(int gridStateId, GeometryListNative& geometryListIn, int nodeIndex);

        /// <summary>
        /// Deletes the closest mesh edge within the search radius from the input point
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="geometryListIn">The input point coordinates</param>
        /// <param name="searchRadius">The search radius</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_delete_edge(int gridStateId, GeometryListNative& geometryListIn, double searchRadius);

        /// <summary>
        /// Deletes the closest mesh edge within the search radius from the input point
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="geometryListIn">The input point coordinates</param>
        /// <param name="searchRadius">The search radius</param>
        /// <param name="edgeIndex">The edge index</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_find_edge(int gridStateId, GeometryListNative& geometryListIn, double searchRadius, int& edgeIndex);

        /// <summary>
        /// Offset a polygon
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="geometryListIn">The coordinate of the offset point</param>
        /// <param name="innerPolygon">Compute inner/outer polygon</param>
        /// <param name="distance">The offset distance</param>
        /// <param name="geometryListOut">The offsetted polygon</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_offsetted_polygon(int gridStateId, GeometryListNative& geometryListIn, bool innerPolygon, double distance, GeometryListNative& geometryListOut);

        /// <summary>
        /// Get the number of vertices of the offsetted polygon  Count the number of vertices after polygon refinment
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="geometryListIn">The coordinate of the offset point</param>
        /// <param name="innerPolygon">Compute inner/outer polygon</param>
        /// <param name="distance">The offset distance</param>
        /// <param name="numberOfPolygonVertices">The number of vertices of the generated polygon</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_offsetted_polygon_count(int gridStateId, GeometryListNative& geometryListIn, bool innerPolygon, double distance, int& numberOfPolygonVertices);

        /// <summary>
        /// Refine a grid based on the samples contained in the geometry list
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="geometryListNative">The sample set</param>
        /// <param name="interpolationParametersNative">The interpolation parameters</param>
        /// <param name="sampleRefineParametersNative">The interpolation settings related to the samples</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_refine_mesh_based_on_samples(int gridStateId, GeometryListNative& geometryListNative, InterpolationParametersNative& interpolationParametersNative, SampleRefineParametersNative& sampleRefineParametersNative);

        /// <summary>
        /// Refine a grid based on polygon
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="geometryListNative">The closed polygon where to perform the refinement</param>
        /// <param name="interpolationParametersNative">The interpolation parameters</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_refine_mesh_based_on_polygon(int gridStateId, GeometryListNative& geometryListNative, InterpolationParametersNative& interpolationParametersNative);

        /// <summary>
        /// Finds the vertex index closest to the input point
        /// </summary>
        /// <param name="gridStateId"></param>
        /// <param name="geometryListIn"></param>
        /// <param name="searchRadius"></param>
        /// <param name="vertexIndex"></param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_get_node_index(int gridStateId, GeometryListNative& geometryListIn, double searchRadius, int& vertexIndex);

        /// <summary>
        /// Selects points in polygons
        /// </summary>
        /// <param name="gridGeomId">Id of the mesh state</param>
        /// <param name="inputPolygon">The polygon(s) used for selection</param>
        /// <param name="inputPoints">The points to select</param>
        /// <param name="selectedPoints">The selected points in the zCoordinates field (0.0 not selected, 1.0 selected)</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_points_in_polygon(int gridGeomId, GeometryListNative& inputPolygon, GeometryListNative& inputPoints, GeometryListNative& selectedPoints);

        /// <summary>
        /// Flip the edges
        /// </summary>
        /// <param name="gridStateId">Id of the mesh state</param>
        /// <param name="isTriangulationRequired">The option to triangulate also non triangular cells (if activated squares becomes triangles) </param>
        /// <param name="isAccountingForLandBoundariesRequired">The option to account for land boundaries</param>
        /// <param name="projectToLandBoundaryOption">The option to determine how to snap to land boundaries</param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_flip_edges(int gridStateId, int isTriangulationRequired, int isAccountingForLandBoundariesRequired, int projectToLandBoundaryOption);

        /// <summary>
        /// Generates curvilinear grid from splines with transfinite interpolation
        /// </summary>
        /// <param name="gridStateId"></param>
        /// <param name="geometryListNativeIn"></param>
        /// <param name="curvilinearParametersNative"></param>
        /// <returns>Error code</returns>
        GRIDGEOM_API int ggeo_curvilinear_mesh_from_splines(int gridStateId, GeometryListNative& geometryListNativeIn, CurvilinearParametersNative& curvilinearParametersNative);

#ifdef __cplusplus
    }
#endif
}