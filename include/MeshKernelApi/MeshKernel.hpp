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

        /// @brief Deallocates mesh state
        /// @param[in] meshKernelId The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_deallocate_state(int meshKernelId);

        /// @brief Deletes a mesh in a polygon using several options
        /// @param[in] meshKernelId   The id of the mesh state
        /// @param[in] polygon        The polygon where to perform the operation
        /// @param[in] deletionOption The deletion option (to be detailed)
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
        /// before passing the struct to `mkernel_get_mesh2d_data`.
        /// @param[in]  meshKernelId The id of the mesh state
        /// @param[out] mesh2d       The Mesh2D dimensions
        /// @returns Error code
        MKERNEL_API int mkernel_get_mesh2d_dimensions(int meshKernelId,
                                                      Mesh2D& mesh2d);

        /// @brief Gets the Mesh2D dimensions data
        ///
        /// This function ought to be called after `mkernel_get_mesh2d_dimensions` has been called
        /// and the pointers have been set to correctly sized memory.
        /// @param[in]     meshKernelId The id of the mesh state
        /// @param[in,out] mesh2d       The Mesh2D data
        /// @returns Error code
        MKERNEL_API int mkernel_get_mesh2d_data(int meshKernelId,
                                                Mesh2D& mesh2d);

        /// @brief Gets the curvilinear grid dimensions as a CurvilinearGrid struct (converted as set of edges and nodes)
        ///
        /// The integer parameters of the CurvilinearGrid struct are set to the corresponding dimensions
        /// The pointers are set to null, and must be set to correctly sized memory
        /// before passing the struct to `mkernel_get_curvilinear_data`.
        /// @param[in]  meshKernelId    The id of the mesh state
        /// @param[out] curvilinearGrid The CurvilinearGrid data
        /// @returns Error code
        MKERNEL_API int mkernel_get_curvilinear_dimensions(int meshKernelId,
                                                           CurvilinearGrid& curvilinearGrid);

        /// @brief Gets the curvilinear grid data as a CurvilinearGrid struct (converted as set of edges and nodes)
        ///
        /// This function ought to be called after `mkernel_get_curvilinear_dimension` has been called
        /// and the pointers have been set to correctly sized memory.
        /// @param[in]  meshKernelId    The id of the mesh state
        /// @param[out] curvilinearGrid The CurvilinearGrid data
        /// @returns Error code
        MKERNEL_API int mkernel_get_curvilinear_data(int meshKernelId,
                                                     CurvilinearGrid& curvilinearGrid);

        /// @brief Gets the Mesh1D dimensions
        ///
        /// The integer parameters of the Mesh1D struct are set to the corresponding dimensions
        /// The pointers are set to null, and must be set to correctly sized memory
        /// before passing the struct to `mkernel_get_mesh1d_data`
        /// @param[in]  meshKernelId The id of the mesh state
        /// @param[out] mesh1d       The Mesh1D dimensions
        /// @returns Error code
        MKERNEL_API int mkernel_get_mesh1d_dimensions(int meshKernelId,
                                                      Mesh1D& mesh1d);

        /// @brief Gets the Mesh1D dimensions data
        ///
        /// This function ought to be called after `mkernel_get_mesh1d_dimensions` has been called
        /// and the pointers have been set to correctly sized memory
        /// @param[in]     meshKernelId The id of the mesh state
        /// @param[in,out] mesh1d       The Mesh1D data
        /// @returns Error code
        MKERNEL_API int mkernel_get_mesh1d_data(int meshKernelId,
                                                Mesh1D& mesh1d);

        /// @brief Gets the meshkernel::Contacts size as a meshkernelapi::Contacts struct
        /// @param[in]  meshKernelId           The id of the mesh state
        /// @param[out] contacts               Contacts data
        /// @returns                           Error code
        MKERNEL_API int mkernel_count_contacts(int meshKernelId,
                                               Contacts& contacts);

        /// @brief Gets the meshkernel::Contacts data as a meshkernelapi::Contacts struct
        /// @param[in]  meshKernelId           The id of the mesh state
        /// @param[out] contacts               Contacts data
        /// @returns                           Error code
        MKERNEL_API int mkernel_get_contacts(int meshKernelId,
                                             Contacts& contacts);

        /// @brief Count the number of hanging edges
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[out] numHangingEdgesIndices
        /// @returns Error code
        MKERNEL_API int mkernel_count_hanging_edges_mesh2d(int meshKernelId, int& numHangingEdgesIndices);

        /// @brief Gets the indices of hanging edges
        /// @param[in]     meshKernelId        The id of the mesh state
        /// @param[in,out] hangingEdgesIndices Pointer to memory where the hanging edge indices will be stored
        /// @returns Error code
        MKERNEL_API int mkernel_get_hanging_edges_mesh2d(int meshKernelId, int** hangingEdgesIndices);

        /// @brief Deletes the hanging edges
        /// @param[in] meshKernelId The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_delete_hanging_edges_mesh2d(int meshKernelId);

        /// @brief Orthogonalization
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] projectToLandBoundaryOption The option to determine how to snap to land boundaries
        /// @param[in] orthogonalizationParameters The structure containing the orthogonalization parameters
        /// @param[in] polygons                    The polygon where to perform the orthogonalization
        /// @param[in] landBoundaries              The land boundaries to account for in the orthogonalization process
        /// @returns Error code
        MKERNEL_API int mkernel_compute_orthogonalization_mesh2d(int meshKernelId,
                                                                 int projectToLandBoundaryOption,
                                                                 const OrthogonalizationParameters& orthogonalizationParameters,
                                                                 const GeometryList& polygons,
                                                                 const GeometryList& landBoundaries);

        /// @brief Orthogonalization initialization (first function to use in interactive mode)
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

        /// @brief Prepares outer orthogonalization iteration (interactive mode)
        /// @param[in] meshKernelId The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_prepare_outer_iteration_orthogonalization_mesh2d(int meshKernelId);

        /// @brief Performs inner orthogonalization iteration (interactive mode)
        /// @param[in] meshKernelId The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_compute_inner_ortogonalization_iteration_mesh2d(int meshKernelId);

        /// @brief Finalizes orthogonalization outer iteration (interactive mode)
        /// @param[in] meshKernelId
        /// @returns Error code
        MKERNEL_API int mkernel_finalize_inner_ortogonalization_iteration_mesh2d(int meshKernelId);

        /// @brief Cleans up back-end orthogonalization algorithm (interactive mode)
        /// @param[in] meshKernelId The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_delete_orthogonalization_mesh2d(int meshKernelId);

        /// @brief Gets the orthogonality
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[out] geometryList The orthogonality values of each edge
        /// @returns Error code
        MKERNEL_API int mkernel_get_orthogonality_mesh2d(int meshKernelId, GeometryList& geometryList);

        /// @brief Gets the smoothness
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[out] geometryList The smoothness values of each edge
        /// @returns Error code
        MKERNEL_API int mkernel_get_smoothness_mesh2d(int meshKernelId, GeometryList& geometryList);

        /// @brief Get spline intermediate points
        /// @param[in] geometryListIn The input corner nodes of the splines
        /// @param[out] geometryListOut The output spline
        /// @param[out] numberOfPointsBetweenNodes The number of spline nodes between the corners points
        /// @returns Error code
        MKERNEL_API int mkernel_get_splines(const GeometryList& geometryListIn,
                                            GeometryList& geometryListOut,
                                            int numberOfPointsBetweenNodes);

        /// @brief Gets the coordinates of the closest existing node
        /// @param[in]  meshKernelId    Id of the grid state
        /// @param[in]  geometryListIn  Node coordinates
        /// @param[in]  searchRadius    The radius where to search for the node
        /// @param[out] geometryListOut Mesh2D node coordinates
        /// @returns Error code
        MKERNEL_API int mkernel_get_closest_node_mesh2d(int meshKernelId,
                                                        const GeometryList& geometryListIn,
                                                        double searchRadius,
                                                        GeometryList& geometryListOut);

        /// @brief Makes a triangular grid in a polygon
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] geometryList The polygon where to triangulate
        /// @returns Error code
        MKERNEL_API int mkernel_make_mesh_from_polygon_mesh2d(int meshKernelId, const GeometryList& geometryList);

        /// @brief Makes a triangular mesh from samples
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] geometryList The samples where to triangulate
        /// @returns Error code
        MKERNEL_API int mkernel_make_mesh_from_samples_mesh2d(int meshKernelId, const GeometryList& geometryList);

        /// @brief Retrieves the mesh boundary polygon
        /// @param[in]  meshKernelId The id of the mesh state
        /// @param[out] geometryList The output network boundary polygon
        /// @returns Error code
        MKERNEL_API int mkernel_get_mesh_boundaries_to_polygon_mesh2d(int meshKernelId, GeometryList& geometryList);

        /// @brief Counts the number of polygon nodes contained in the mesh boundary polygon
        /// @param[in]  meshKernelId         The id of the mesh state
        /// @param[out] numberOfPolygonNodes The number of polygon points
        /// @returns Error code
        MKERNEL_API int mkernel_count_mesh_boundaries_to_polygon_mesh2d(int meshKernelId,
                                                                        int& numberOfPolygonNodes);

        /// @brief Gets the refined polygon
        /// @param[in]  meshKernelId   The id of the mesh state
        /// @param[in]  geometryListIn The input polygons
        /// @param[in]  firstIndex     The index of the first node
        /// @param[in]  secondIndex    The index of the second node
        /// @param[in]  distance       The refinement distance
        /// @param[out] geometryListOut
        /// @returns Error code
        MKERNEL_API int mkernel_refine_polygon(int meshKernelId,
                                               const GeometryList& geometryListIn,
                                               int firstIndex,
                                               int secondIndex,
                                               double distance,
                                               GeometryList& geometryListOut);

        /// @brief Counts the number of nodes after polygon refinement
        /// @param[in] meshKernelId   The id of the mesh state
        /// @param[in] geometryListIn The input polygon
        /// @param[in] firstIndex     The index of the first node
        /// @param[in] secondIndex    The index of the second node
        /// @param[in] distance       The refinement distance
        /// @param[out] numberOfPolygonNodes The number of nodes after refinement
        /// @returns Error code
        MKERNEL_API int mkernel_count_refine_polygon(int meshKernelId,
                                                     const GeometryList& geometryListIn,
                                                     int firstIndex,
                                                     int secondIndex,
                                                     double distance,
                                                     int& numberOfPolygonNodes);

        /// @brief Merges nodes within a distance of 0.001 m, effectively removing small edges
        /// @param[in] meshKernelId   The id of the mesh state
        /// @param[in] geometryListIn The polygon where to perform the operation
        /// @returns Error code
        MKERNEL_API int mkernel_merge_nodes_mesh2d(int meshKernelId, const GeometryList& geometryListIn);

        /// @brief Merges node \p startNode to \p endNode
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] startNode    The index of the first node to merge
        /// @param[in] endNode      The index of the second node to merge
        /// @returns Error code
        MKERNEL_API int mkernel_merge_two_nodes_mesh2d(int meshKernelId, int startNode, int endNode);

        /// @brief Gets the selected mesh node indices
        /// @param[in]  meshKernelId   The id of the mesh state
        /// @param[in]  geometryListIn The input polygons
        /// @param[in]  inside         Count nodes indices inside (1) or outside (0) the polygon
        /// @param[out] selectedNodes  The selected nodes indices
        /// @returns Error code
        MKERNEL_API int mkernel_get_nodes_in_polygons(int meshKernelId,
                                                      const GeometryList& geometryListIn,
                                                      int inside,
                                                      int** selectedNodes);

        /// @brief Counts the number of selected mesh node indices
        /// @param[in]  meshKernelId      The id of the mesh state
        /// @param[in]  geometryListIn    The input polygons
        /// @param[in]  inside            Count nodes inside (1) or outside (0) the polygon
        /// @param[out] numberOfMeshNodes The number of selected nodes
        /// @returns Error code
        MKERNEL_API int mkernel_count_nodes_in_polygons(int meshKernelId,
                                                        const GeometryList& geometryListIn,
                                                        int inside,
                                                        int& numberOfMeshNodes);

        /// @brief Insert a new edge connecting \p startNode and \p endNode
        /// @param[in]  meshKernelId The id of the mesh state
        /// @param[in]  startNode    The index of the first node to connect
        /// @param[in]  endNode      The index of the second node to connect
        /// @param[out] newEdgeIndex The index of the new edge
        /// @returns Error code
        MKERNEL_API int mkernel_insert_edge_mesh2d(int meshKernelId, int startNode, int endNode, int& newEdgeIndex);

        /// @brief Inserts a new node
        /// @param[in]  meshKernelId The id of the mesh state
        /// @param[in]  xCoordinate  X-coordinate of the new node
        /// @param[in]  yCoordinate  Y-coordinate of the new node
        /// @param[out] nodeIndex    The index of the new mesh node
        /// @returns Error code
        MKERNEL_API int mkernel_insert_node_mesh2d(int meshKernelId,
                                                   double xCoordinate,
                                                   double yCoordinate,
                                                   int& nodeIndex);

        /// @brief Deletes a node with specified \p nodeIndex
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] nodeIndex    The nodeIndex to delete
        /// @returns Error code
        MKERNEL_API int mkernel_delete_node_mesh2d(int meshKernelId, int nodeIndex);

        /// @brief Moves a selected node to a new position
        /// @param[in] meshKernelId   The id of the mesh state
        /// @param[in] geometryListIn The new coordinate
        /// @param[in] nodeIndex      The node index (to be detailed)
        /// @returns Error code
        MKERNEL_API int mkernel_move_node_mesh2d(int meshKernelId, const GeometryList& geometryListIn, int nodeIndex);

        /// @brief Deletes the closest mesh edge within the search radius from the input point
        /// @param[in] meshKernelId   The id of the mesh state
        /// @param[in] geometryListIn The input point coordinates
        /// @returns Error code
        MKERNEL_API int mkernel_delete_edge_mesh2d(int meshKernelId, const GeometryList& geometryListIn);

        /// @brief Deletes the closest mesh edge within the search radius from the input point
        /// @param[in] meshKernelId   The id of the mesh state
        /// @param[in] geometryListIn The input point coordinates
        /// @param[out] edgeIndex     The edge index
        /// @returns Error code
        MKERNEL_API int mkernel_find_edge_mesh2d(int meshKernelId, const GeometryList& geometryListIn, int& edgeIndex);

        /// @brief Offsets a polygon
        /// @param[in] meshKernelId   The id of the mesh state
        /// @param[in] geometryListIn The polygon to be offsetted
        /// @param[in] innerPolygon   Compute inner (true) or outer (false) polygon
        /// @param[in] distance       The offset distance
        /// @param[out] geometryListOut The offsetted polygon
        /// @returns Error code
        MKERNEL_API int mkernel_get_offsetted_polygon(int meshKernelId,
                                                      const GeometryList& geometryListIn,
                                                      bool innerPolygon,
                                                      double distance,
                                                      GeometryList& geometryListOut);

        /// @brief Gets the number of nodes of the offsetted polygon  Count the number of nodes after polygon refinement
        /// @param[in] meshKernelId   The id of the mesh state
        /// @param[in] geometryListIn The polygon to be offsetted
        /// @param[in] innerPolygon   Compute inner (true) or outer (false) polygon
        /// @param[in] distance       The offset distance
        /// @param[out] numberOfPolygonNodes The number of nodes of the generated polygon
        /// @returns Error code
        MKERNEL_API int mkernel_count_offsetted_polygon(int meshKernelId,
                                                        const GeometryList& geometryListIn,
                                                        bool innerPolygon,
                                                        double distance,
                                                        int& numberOfPolygonNodes);

        /// @brief Refines a grid based on the samples contained in the geometry list
        /// @param[in] meshKernelId            The id of the mesh state
        /// @param[in] geometryList            The sample set
        /// @param[in] interpolationParameters The interpolation parameters
        /// @param[in] sampleRefineParameters  The interpolation settings related to the samples
        /// @returns Error code
        MKERNEL_API int mkernel_refine_based_on_samples_mesh2d(int meshKernelId,
                                                               const GeometryList& geometryList,
                                                               const InterpolationParameters& interpolationParameters,
                                                               const SampleRefineParameters& sampleRefineParameters);

        /// @brief Refines a grid based on polygon
        /// @param[in] meshKernelId            The id of the mesh state
        /// @param[in] geometryList            The closed polygon where to perform the refinement
        /// @param[in] interpolationParameters The interpolation parameters
        /// @returns Error code
        MKERNEL_API int mkernel_refine_based_on_polygon_mesh2d(int meshKernelId,
                                                               const GeometryList& geometryList,
                                                               const InterpolationParameters& interpolationParameters);

        /// @brief Finds the node index closest to the input point
        /// @param[in] meshKernelId   The id of the mesh state
        /// @param[in] geometryListIn The input point from where starting the search
        /// @param[in] searchRadius   The search radius to use for the search
        /// @param[out] nodeIndex     The index of the found node
        /// @returns Error code
        MKERNEL_API int mkernel_get_node_index_mesh2d(int meshKernelId,
                                                      const GeometryList& geometryListIn,
                                                      double searchRadius,
                                                      int& nodeIndex);

        /// @brief Selects points in polygons
        /// @param[in]  meshKernelId   The id of the mesh state
        /// @param[in]  inputPolygon   The polygon(s) used for selection
        /// @param[in]  inputPoints    The points to select
        /// @param[out] selectedPoints The selected points in the zCoordinates field (0.0 not selected, 1.0 selected)
        /// @returns Error code
        MKERNEL_API int mkernel_get_points_in_polygon(int meshKernelId,
                                                      const GeometryList& inputPolygon,
                                                      const GeometryList& inputPoints,
                                                      GeometryList& selectedPoints);

        /// @brief Flips the edges
        /// @param[in] meshKernelId                The id of the mesh state
        /// @param[in] isTriangulationRequired     The option to triangulate also non triangular cells (if activated squares becomes triangles)
        /// @param[in] projectToLandBoundaryOption The option to determine how to snap to land boundaries
        /// @returns Error code
        MKERNEL_API int mkernel_flip_edges_mesh2d(int meshKernelId,
                                                  int isTriangulationRequired,
                                                  int projectToLandBoundaryOption);

        /// @brief Gets the number of obtuse triangles (those having one edge longer than the sum of the other two)
        /// @param[in]  meshKernelId       The id of the mesh state
        /// @param[out] numObtuseTriangles The number of obtuse triangles
        /// @return
        MKERNEL_API int mkernel_count_obtuse_triangles_mesh2d(int meshKernelId, int& numObtuseTriangles);

        /// @brief Gets the obtuse triangle mass centers (those having one edge longer than the sum of the other two)
        /// @param[in]  meshKernelId  The id of the mesh state
        /// @param[out] result        The obtuse triangles mass centers
        /// @return Error code (0 Successful)
        MKERNEL_API int mkernel_get_obtuse_triangles_mass_centers_mesh2d(int meshKernelId, GeometryList& result);

        /// @brief Counts the small flow edges (flow edges are the edges connecting the face circumcenters)
        /// @param[in] meshKernelId            The id of the mesh state
        /// @param[in] smallFlowEdgesThreshold The configurable threshold for detecting the small flow edges
        /// @param[out] numSmallFlowEdges      The number of the small flow edges
        /// @return
        MKERNEL_API int mkernel_count_small_flow_edge_centers_mesh2d(int meshKernelId,
                                                                     double smallFlowEdgesThreshold,
                                                                     int& numSmallFlowEdges);

        /// @brief Gets the small flow edges (flow edges are the edges connecting the face circumcenters)
        /// @param[in] meshKernelId            The id of the mesh state
        /// @param[in] smallFlowEdgesThreshold The configurable threshold for detecting the small flow edges
        /// @param[out] result                 The center points of the small flow edges
        /// @return Error code (0 Successful)
        MKERNEL_API int mkernel_get_small_flow_edge_centers_mesh2d(int meshKernelId,
                                                                   double smallFlowEdgesThreshold,
                                                                   GeometryList& result);

        /// @brief Deletes the small flow edges (flow edges are the edges connecting the face circumcenters)
        /// @param[in] meshKernelId               The id of the mesh state
        /// @param[in] smallFlowEdgesThreshold    The configurable threshold for detecting the small flow edges
        /// @param[in] minFractionalAreaTriangles The ratio of the face area to the average area of neighboring non triangular faces
        /// @return Error code (0 Successful)
        MKERNEL_API int mkernel_delete_small_flow_edges_mesh2d(int meshKernelId,
                                                               double smallFlowEdgesThreshold,
                                                               double minFractionalAreaTriangles);

        /// @brief Computes 1d-2d contacts, where every single 1d node is connected to one 2d face circumcenter (ggeo_make1D2Dinternalnetlinks_dll)
        ///
        /// \see meshkernel::Contacts::ComputeSingleContacts
        /// @param[in]  meshKernelId  The id of the mesh state
        /// @param[in]  oneDNodeMask  The mask to apply to 1d nodes (1 = connect node, 0 = do not generate contacts)
        /// @param[in]  polygons      The polygons selecting the area where the 1d-2d contacts will be generated.
        /// @return                   Error code (0 Successful)
        MKERNEL_API int mkernel_compute_single_contacts(int meshKernelId,
                                                        const int* oneDNodeMask,
                                                        const GeometryList& polygons);

        /// @brief Computes 1d-2d contacts, where a single 1d node is connected to multiple 2d face circumcenters (ggeo_make1D2Dembeddedlinks_dll)
        ///
        /// \see meshkernel::Contacts::ComputeMultipleContacts
        /// @param[in]  meshKernelId  The id of the mesh state
        /// @param[in]  oneDNodeMask  The mask to apply to 1d nodes (1 = connect node, 0 = do not generate contacts)
        /// @return                   Error code (0 Successful)
        MKERNEL_API int mkernel_compute_multiple_contacts(int meshKernelId,
                                                          const int* oneDNodeMask);

        /// @brief Computes 1d-2d contacts, where a 2d face per polygon is connected to the closest 1d node (ggeo_make1D2Droofgutterpipes_dll)
        ///
        /// \see meshkernel::Contacts::ComputeContactsWithPolygons
        /// @param[in]  meshKernelId  The id of the mesh state
        /// @param[in]  oneDNodeMask  The mask to apply to 1d nodes (1 = connect node, 0 = do not generate contacts)
        /// @param[in]  polygons      The polygons to connect
        /// @return                   Error code (0 Successful)
        MKERNEL_API int mkernel_compute_with_polygons_contacts(int meshKernelId,
                                                               const int* oneDNodeMask,
                                                               const GeometryList& polygons);

        /// @brief Computes 1d-2d contacts, where 1d nodes are connected to the 2d faces mass centers containing the input point (ggeo_make1D2Dstreetinletpipes_dll)
        ///
        /// \see meshkernel::Contacts::ComputeContactsWithPoints
        /// @param[in]  meshKernelId  The id of the mesh state
        /// @param[in]  oneDNodeMask  The mask to apply to 1d nodes (1 = connect node, 0 = do not generate contacts)
        /// @param[in]  points        The points selecting the faces to connect
        /// @return                   Error code (0 Successful)
        MKERNEL_API int mkernel_compute_with_points_contacts(int meshKernelId,
                                                             const int* oneDNodeMask,
                                                             const GeometryList& points);

        /// @brief Computes 1d-2d contacts, where 1d nodes are connected to the closest 2d faces at the boundary (ggeo_make1D2DRiverLinks_dll)
        ///
        /// \see meshkernel::Contacts::ComputeBoundaryContacts
        /// @param[in]  meshKernelId The id of the mesh state.
        /// @param[in]  oneDNodeMask The mask to apply to 1d nodes (1 = connect node, 0 = do not generate contacts)
        /// @param[in]  polygons     The points selecting the faces to connect.
        /// @param[in]  searchRadius The radius used for searching neighboring faces, if equal to doubleMissingValue, the search radius will be calculated internally.
        /// @return                  Error code (0 Successful)
        MKERNEL_API int mkernel_compute_boundary_contacts(int meshKernelId,
                                                          const int* oneDNodeMask,
                                                          const GeometryList& polygons,
                                                          double searchRadius);

        /// @brief Curvilinear grid refinement
        ///
        /// \p geometryListFirstPoint and \p geometryListSecondPoint must lie on the same gridline
        /// @param[in] meshKernelId            The id of the mesh state.
        /// @param[in] geometryListFirstPoint  The geometry list containing the first node of the segment defining the refinement zone.
        /// @param[in] geometryListSecondPoint The geometry list containing the second node of the segment defining the refinement zone.
        /// @param[in] refinement              The number of refinement lines between \p geometryListFirstPoint and \p geometryListSecondPoint
        /// @return Error code (0 Successful)
        MKERNEL_API int mkernel_refine_curvilinear(int meshKernelId,
                                                   const GeometryList& geometryListFirstPoint,
                                                   const GeometryList& geometryListSecondPoint,
                                                   int refinement);

        /// @brief Curvilinear grid derefinement
        ///
        /// \p geometryListFirstPoint and \p geometryListSecondPoint must lie on the same gridline
        /// @param meshKernelId            The id of the mesh state.
        /// @param geometryListFirstPoint  The geometry list containing the first node of the segment defining the derefinement zone.
        /// @param geometryListSecondPoint The geometry list containing the second node of the segment defining the derefinement zone.
        /// @return Error code (0 Successful)
        MKERNEL_API int mkernel_derefine_curvilinear(int meshKernelId,
                                                     const GeometryList& geometryListFirstPoint,
                                                     const GeometryList& geometryListSecondPoint);

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
        MKERNEL_API int mkernel_compute_orthogonal_curvilinear(int meshKernelId,
                                                               const GeometryList& geometryList,
                                                               const CurvilinearParameters& curvilinearParameters,
                                                               const SplinesToCurvilinearParameters& splinesToCurvilinearParameters);

        /// @brief Generates a curvilinear grid from splines with the advancing front method. Initialization step (interactive)
        /// @param[in] meshKernelId                    The id of the mesh state
        /// @param[in] geometryList                    The input splines corners
        /// @param[in] curvilinearParameters           The input parameters to generate the curvilinear grid
        /// @param[in] splinesToCurvilinearParameters  The parameters of the advancing front algorithm
        /// @returns Error code
        MKERNEL_API int mkernel_initialize_orthogonal_curvilinear(int meshKernelId,
                                                                  const GeometryList& geometryList,
                                                                  const CurvilinearParameters& curvilinearParameters,
                                                                  const SplinesToCurvilinearParameters& splinesToCurvilinearParameters);

        /// @brief One advancement of the front in curvilinear grid from splines (interactive)
        /// @param[in] meshKernelId The id of the mesh state
        /// @param[in] layer        The layer index
        /// @returns Error code
        MKERNEL_API int mkernel_iterate_orthogonal_curvilinear(int meshKernelId, int layer);

        /// @brief Converts curvilinear grid to mesh and refreshes the state (interactive)
        /// @param[in] meshKernelId The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_refresh_orthogonal_curvilinear(int meshKernelId);

        /// @brief Finalizes curvilinear grid from splines algorithm
        /// @param[in] meshKernelId The id of the mesh state
        /// @returns Error code
        MKERNEL_API int mkernel_delete_orthogonal_curvilinear(int meshKernelId);

        /// @brief Makes a new mesh
        /// @param[in] meshKernelId       The id of the mesh state
        /// @param[in] makeGridParameters The structure containing the make grid parameters
        /// @param[in] geometryList       The polygon to account for
        /// @returns Error code
        MKERNEL_API int mkernel_make_uniform_curvilinear(int meshKernelId,
                                                         const MakeMeshParameters& makeGridParameters,
                                                         const GeometryList& geometryList);

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
        /// @return Error code (0 Successful)
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
        /// @param[in] mesh2d             Mesh2D data
        /// @param[in] startIndex         Mesh2D data start index (not used)
        /// @param[in] samplesXCoordinate The sample x-coordinates
        /// @param[in] samplesYCoordinate The sample y-coordinates
        /// @param[in] samplesValue       The sample values
        /// @param[in] numSamples         The number of samples
        /// @param[out] results           The interpolation results
        /// @param[in] locationType       The location type (see MeshLocations enum)
        /// @param[in] Wu1Duni            A setting for 1d meshes (not used)
        /// @param[in] averagingMethod    The averaging method (see Method enum)
        /// @param[in] minNumberOfSamples The minimum amount of samples (not used)
        /// @param[in] relativeSearchSize The relative search size around the location (larger increases the number of samples considered)
        /// @param[in] spherical          Current projection (0 cartesian, 1 spherical)
        /// @param[in] sphericalAccurate  Accurate spherical computations (0 default spherical, 1 spherical accurate)
        /// @return Error code (0 Successful)
        MKERNEL_API int averaging(const Mesh2D& mesh2d,
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
        /// @param[out] type         The entity type (node, edge or face, see MeshLocations)
        /// @returns Error code
        MKERNEL_API int mkernel_get_geometry_error(int& invalidIndex, int& type);

#ifdef __cplusplus
    }
#endif
} // namespace meshkernelapi
