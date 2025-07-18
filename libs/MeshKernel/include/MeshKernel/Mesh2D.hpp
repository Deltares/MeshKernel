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
#include <array>
#include <ranges>
#include <utility>
#include <vector>

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Definitions.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/Polygon.hpp>
#include <MeshKernel/UndoActions/CompoundUndoAction.hpp>
#include <MeshKernel/UndoActions/SphericalCoordinatesOffsetAction.hpp>
#include <MeshKernel/UndoActions/UndoAction.hpp>

/// \namespace meshkernel
/// @brief Contains the logic of the C++ static library
namespace meshkernel
{
    // Forward declarations
    class CurvilinearGrid;
    class Polygons;
    class GeometryList;

    /// @brief A class derived from Mesh, which describes unstructures 2d meshes.
    ///
    /// When communicating with the client only unstructured meshes are used.
    /// Some algorithms generate curvilinear grids, but these are converted to a mesh
    /// instance when communicating with the client.
    class Mesh2D final : public Mesh
    {
    public:
        using Mesh::CommitAction;
        using Mesh::RestoreAction;

        /// Enumerator describing the different options to delete a mesh
        enum DeleteMeshOptions
        {
            InsideNotIntersected = 0,
            InsideAndIntersected = 1,
            FacesWithIncludedCircumcenters = 2
        };

        /// Enumerator for different properties on a 2D mesh
        enum class Property
        {
            Orthogonality = 0,
            EdgeLength = 1
        };

        /// @brief Default destructor
        ~Mesh2D() override = default;

        /// @brief Default constructor
        Mesh2D();

        /// @brief Construct a mesh2d using only the projection
        /// @param[in] projection The projection to use
        explicit Mesh2D(Projection projection);

        /// @brief Construct a mesh2d starting from the edges and nodes
        /// @param[in] edges The input edges
        /// @param[in] nodes The input nodes
        /// @param[in] projection The projection to use
        Mesh2D(const std::vector<Edge>& edges,
               const std::vector<Point>& nodes,
               Projection projection);

        /// @brief Construct a mesh2d from face nodes and num face nodes
        /// @param[in] edges The input edges
        /// @param[in] nodes The input nodes
        /// @param[in] faceNodes The input face nodes
        /// @param[in] numFaceNodes For each face, the number of nodes
        /// @param[in] projection The mesh projection
        Mesh2D(const std::vector<Edge>& edges,
               const std::vector<Point>& nodes,
               const std::vector<std::vector<UInt>>& faceNodes,
               const std::vector<std::uint8_t>& numFaceNodes,
               Projection projection);

        /// @brief Create triangular grid from nodes (triangulatesamplestonetwork)
        /// @param[in] nodes Input nodes
        /// @param[in] polygons Selection polygon
        /// @param[in] projection The projection to use
        Mesh2D(const std::vector<Point>& nodes, const Polygons& polygons, Projection projection);

        /// @brief Perform complete administration
        void Administrate(CompoundUndoAction* undoAction = nullptr) override;

        /// @brief Compute face circumcenters
        void ComputeCircumcentersMassCentersAndFaceAreas(bool computeMassCenters = false);

        /// @brief Constructs the face nodes mapping, face mass centers and areas
        void FindFaces();

        /// @brief Find remaining face information given the face nodes mapping
        /// @param[in] faceNodes The input face nodes
        /// @param[in] numFaceNodes For each face, the number of nodes
        void FindFacesGivenFaceNodesMapping(const std::vector<std::vector<UInt>>& faceNodes,
                                            const std::vector<std::uint8_t>& numFaceNodes);

        /// @brief Offset the x coordinates if m_projection is spherical
        /// @param[in] minx
        /// @param[in] maxx
        [[nodiscard]] std::unique_ptr<SphericalCoordinatesOffsetAction> OffsetSphericalCoordinates(double minx, double maxx);

        /// @brief Apply the coordinate offset action
        void CommitAction(const SphericalCoordinatesOffsetAction& undoAction);

        /// @brief Undo the coordinate offset action
        ///
        /// Restore mesh to state before coordinate offset action was applied
        void RestoreAction(const SphericalCoordinatesOffsetAction& undoAction);

        /// @brief For a face create a closed polygon and fill local mapping caches (get_cellpolygon)
        /// @param[in]  faceIndex              The face index
        /// @param[out] polygonNodesCache      The node cache array filled with the nodes values
        /// @param[out] localNodeIndicesCache  The consecutive node index in polygonNodesCache (0, 1, 2,...)
        /// @param[out] globalEdgeIndicesCache The edge cache array filled with the global edge indices
        void ComputeFaceClosedPolygonWithLocalMappings(UInt faceIndex,
                                                       std::vector<Point>& polygonNodesCache,
                                                       std::vector<UInt>& localNodeIndicesCache,
                                                       std::vector<UInt>& globalEdgeIndicesCache) const;

        /// @brief For a closed polygon, compute the circumcenter of a face (getcircumcenter)
        /// @param[in,out] polygon       Cache storing the face nodes
        /// @param[in]     edgesNumFaces For meshes, the number of faces sharing the edges
        /// @returns       The computed circumcenter
        [[nodiscard]] Point ComputeFaceCircumenter(std::vector<Point>& polygon,
                                                   const std::vector<UInt>& edgesNumFaces) const;

        /// @brief Gets the mass centers of obtuse triangles
        /// @returns The center of obtuse triangles
        [[nodiscard]] std::vector<Point> GetObtuseTrianglesCenters();

        /// @brief Gets the edges crossing the small flow edges
        /// @param[in] smallFlowEdgesThreshold The configurable threshold for detecting the small flow edges
        /// @returns The indices of the edges crossing small flow edges
        [[nodiscard]] std::vector<UInt> GetEdgesCrossingSmallFlowEdges(double smallFlowEdgesThreshold);

        /// @brief Gets the flow edges centers from the crossing edges
        /// @param[in] edges The crossing edges indices
        /// @returns The centers of the flow edges
        [[nodiscard]] std::vector<Point> GetFlowEdgesCenters(const std::vector<UInt>& edges) const;

        /// @brief Deletes small flow edges (removesmallflowlinks, part 1)
        ///
        /// An unstructured mesh can be used to calculate water flow. This involves
        /// a pressure gradient between the circumcenters of neighbouring faces.
        /// That procedure is numerically unreliable when the distance between face
        /// circumcenters (flow edges) becomes too small. Let's consider the following figure
        /// \image html coincide_circumcenter.svg  "Coincide circumcenter"
        /// The algorithm works as follow:
        ///
        /// -   Any degenerated triangle (e.g. those having a coinciding node) is
        ///     removed by collapsing the second and third node into the first one.
        ///
        /// -   The edges crossing small flow edges are found. The flow edge length
        ///     is computed from the face circumcenters and compared to an estimated
        ///     cut off distance. The cutoff distance is computed using the face
        ///     areas as follow:
        ///
        ///     \f$\textrm{cutOffDistance} = \textrm{threshold} \cdot 0.5 \cdot (\sqrt{\textrm{Area}_I}+\sqrt{\textrm{Area}_{II}})\f$
        ///
        /// -   All small flow edges are flagged with invalid indices and removed
        ///     from the mesh. Removal occors in the \ref Mesh2D::Administrate method.
        /// @param[in] smallFlowEdgesThreshold The configurable threshold for detecting the small flow edges
        [[nodiscard]] std::unique_ptr<meshkernel::UndoAction> DeleteSmallFlowEdges(double smallFlowEdgesThreshold);

        /// @brief Deletes small triangles at the boundaries (removesmallflowlinks, part 2)
        ///
        /// This algorithm removes triangles having the following properties:
        /// - The are at mesh boundary.
        ///
        /// - One or more neighboring faces are non-triangles.
        ///
        /// - The ratio of the face area to the average area of neighboring non
        ///   triangles is less than a minimum ratio (defaults to 0.2).
        ///
        /// - The absolute cosine of one internal angle is less than 0.2.
        ///
        /// These triangles having the above properties are merged by collapsing the
        /// face nodes to the node having the minimum absolute cosine (e.g. the node
        /// where the internal angle is closer to 90 degrees).
        /// @param[in] minFractionalAreaTriangles Small triangles at the boundaries will be eliminated.
        /// This threshold is the ration of the face area to the average area of neighboring faces.
        [[nodiscard]] std::unique_ptr<UndoAction> DeleteSmallTrianglesAtBoundaries(double minFractionalAreaTriangles);

        /// @brief Computes node neighbours
        void ComputeNodeNeighbours(std::vector<std::vector<UInt>>& nodesNodes, UInt& maxNumNeighbours) const;

        /// @brief Gets the aspect ratios (the ratios edges lengths to flow edges lengths)
        /// @param[in,out] aspectRatios The aspect ratios (passed as reference to avoid re-allocation)
        void ComputeAspectRatios(std::vector<double>& aspectRatios) const;

        ///  @brief Classifies the nodes (makenetnodescoding)
        void ClassifyNodes();

        /// @brief Get the node type
        MeshNodeType GetNodeType(const UInt nodeId) const { return m_nodesTypes[nodeId]; }

        /// @brief Get the node type
        void GetNodeTypes(std::vector<MeshNodeType>& nodeTypes) const { nodeTypes = m_nodesTypes; }

        /// @brief Deletes coinciding triangles
        [[nodiscard]] std::unique_ptr<UndoAction> DeleteDegeneratedTriangles();

        /// @brief Transform non-triangular faces to triangular faces
        [[nodiscard]] std::unique_ptr<UndoAction> TriangulateFaces();

        /// @brief Transform non-triangular faces inside a polygon to triangular faces
        [[nodiscard]] std::unique_ptr<UndoAction> TriangulateFaces(const Polygons& polygon);

        /// @brief Make a dual face around the node, enlarged by a factor
        /// @param[in] edgeCentres Centre point of each of the edges.
        /// @param[in] node The node index
        /// @param[in] enlargementFactor The factor by which the dual face is enlarged
        /// @param[out] dualFace The dual face to be calculated
        void MakeDualFace(const std::span<const Point> edgeCentres,
                          UInt node,
                          double enlargementFactor,
                          std::vector<Point>& dualFace) const;

        /// @brief Sorts the faces around a node, sorted in counter clock wise order
        /// @param[in] node The node index
        /// @return The face indexses
        [[nodiscard]] std::vector<UInt> SortedFacesAroundNode(UInt node) const;

        /// @brief Convert all mesh boundaries to a vector of polygon nodes, including holes (copynetboundstopol)
        /// @param[in] polygon The polygon where the operation is performed
        /// @return The resulting polygon mesh boundary
        [[nodiscard]] std::vector<Point> ComputeBoundaryPolygons(const std::vector<Point>& polygon);

        /// @brief Constructs a polygon from the meshboundary, by walking through the mesh
        /// @param[in] polygon The input polygon
        /// @param[in,out] isVisited the visited mesh nodes
        /// @param[in,out] currentNode the current node
        /// @param[out] meshBoundaryPolygon The resulting polygon points
        void WalkBoundaryFromNode(const Polygon& polygon,
                                  std::vector<bool>& isVisited,
                                  UInt& currentNode,
                                  std::vector<Point>& meshBoundaryPolygon) const;

        /// @brief Gets the hanging edges
        /// @return A vector with the indices of the hanging edges
        [[nodiscard]] std::vector<UInt> GetHangingEdges() const;

        /// @brief Deletes the hanging edges
        [[nodiscard]] std::unique_ptr<UndoAction> DeleteHangingEdges();

        /// @brief For a collection of points, compute the face indices including them.
        /// @param[in] points The input point vector.
        /// @return The face indices including the points.
        [[nodiscard]] std::vector<UInt> PointFaceIndices(const std::vector<Point>& points);

        /// @brief Deletes a mesh in a polygon, using several options (delnet)
        /// @param[in] polygon        The polygon where to perform the operation
        ///                           If this Polygons instance contains multiple polygons, the first one will be taken.
        /// @param[in] deletionOption The deletion option
        /// @param[in] invertDeletion Inverts the selected node to delete (instead of outside the polygon, inside the polygon)
        [[nodiscard]] std::unique_ptr<UndoAction> DeleteMesh(const Polygons& polygon, DeleteMeshOptions deletionOption, bool invertDeletion);

        /// @brief  This method generates a mask indicating which locations are within the specified  range of the given metric.
        ///
        /// @param[in] location The location representing the location where to filter the object.
        /// @param[in] property The property by which to filter locations.
        /// @param[in] minValue The minimum value of the metric for filtering.
        /// @param[in] maxValue The maximum value of the metric for filtering.
        /// @ return A vector of boolean values. Each element corresponds to a location and is `true` if the location's metric is within the specified range, and `false` otherwise.
        [[nodiscard]] std::vector<bool> FilterBasedOnMetric(Location location, Property property, double minValue, double maxValue) const;

        /// @brief Inquire if a segment is crossing a face
        /// @param[in] firstPoint The first point of the segment
        /// @param[in] secondPoint The second point of the segment
        /// @return A tuple with the intersectedFace face index and intersected  edge index
        [[nodiscard]] std::tuple<UInt, UInt> IsSegmentCrossingABoundaryEdge(const Point& firstPoint, const Point& secondPoint) const;

        /// @brief Masks the edges of all faces entirely included in all polygons
        /// @param[in] polygons The selection polygon
        /// @param[in] invertSelection Invert selection
        /// @param[in] includeIntersected Included the edges intersected by the polygon
        /// @return The edge mask
        [[nodiscard]] std::vector<int> MaskEdgesOfFacesInPolygon(const Polygons& polygons, bool invertSelection, bool includeIntersected) const;

        /// @brief From the edge mask compute the node mask
        /// @param[in] edgeMask The edge mask
        /// @return The node mask
        [[nodiscard]] std::vector<int> NodeMaskFromEdgeMask(std::vector<int> const& edgeMask) const;

        /// @brief Mask all nodes included in all polygons
        /// @param[in] polygons The input polygon
        /// @param[in] inside   Inside or outside option
        /// @return The node mask
        [[nodiscard]] std::vector<int> NodeMaskFromPolygon(const Polygons& polygons, bool inside) const;

        /// @brief Find edge on the opposite side of the element
        /// @note Currently only valid of quadrilateral elements.
        /// Will throw exception NotImplementedError for non-quadrilateral element shapes.
        UInt FindOppositeEdge(const UInt faceId, const UInt edgeId) const;

        /// @brief Get the next face adjacent to the edge on the opposite side.
        /// @param [in] faceId The starting face
        /// @param [in] edgeId The starting edge
        /// @return Id of neighbour face along the edge
        UInt NextFace(const UInt faceId, const UInt edgeId) const;

        /// @brief Merges mesh connectivity.
        ///
        /// Only merges the mesh connectivity graphs and updates indices.
        /// @note Does not do any administration on the node, edges or elements,
        /// it may be required to call Administrate after merging
        static std::unique_ptr<Mesh2D> Merge(const Mesh2D& mesh1, const Mesh2D& mesh2);

        /// @brief Merges mesh node and edge connectivity into a single mesh.
        static std::unique_ptr<Mesh2D> Merge(const std::span<const Point>& mesh1Nodes,
                                             const std::span<const Edge>& mesh1Edges,
                                             const std::span<const Point>& mesh2Nodes,
                                             const std::span<const Edge>& mesh2Edges,
                                             const Projection projection);

        /// @brief Get the mesh bounding box
        ///
        /// @return The mesh bounding box
        [[nodiscard]] BoundingBox GetBoundingBox() const;

        /// @brief Get the bounding boxes of the mesh edges
        ///
        /// @return The mesh edges bounding boxes
        [[nodiscard]] std::vector<BoundingBox> GetEdgesBoundingBoxes() const;

        /// @brief Find all faces that have the given node as a vertex.
        ///
        /// @param [in] nodeIndex Index of the node
        /// @param [out] sharedFaces On exit will contain only indices of faces that contain nodeIndex as a node.
        void FindFacesConnectedToNode(UInt nodeIndex, std::vector<UInt>& sharedFaces) const;

        /// @brief Get indices of all nodes that are connected directly to a give node along connected edges
        ///
        /// @param [in] nodeIndex Index of the node
        /// @param [out] connectedNodes
        void GetConnectingNodes(UInt nodeIndex, std::vector<UInt>& connectedNodes) const;

        /// @brief Find all unique nodes.
        ///
        /// @param [in] nodeIndex Index of the node
        /// @param [in] sharedFaces List of faces that share the nodeIndex as a common node
        /// @param [in, out] connectedNodes List of nodes that are in the patch of shared faces
        /// @param [out] faceNodeMapping Mapping from node index to the position in connectedNodes list.
        void FindNodesSharedByFaces(UInt nodeIndex, const std::vector<UInt>& sharedFaces, std::vector<UInt>& connectedNodes, std::vector<std::vector<UInt>>& faceNodeMapping) const;

        /// @brief Determine if the node is at the start or end of the edge.
        ///
        /// Returns 0 when the node is at the start of the edge, 1 when it is at the end
        /// and the null value when the edge is not connected to the node.
        UInt IsStartOrEnd(const UInt edgeId, const UInt nodeId) const;

        /// @brief Determine if the element lies on the left or right side of the edge
        ///
        /// Returns 0 when the element is on the left and 1 when it is on the right.
        /// If one or other edge is not connected to the element then a null value will be returned.
        UInt IsLeftOrRight(const UInt elementId, const UInt edgeId) const;

        /// @brief Find the id of the element that is common to both edges.
        ///
        /// If no such element can be found then the null value will be returned.
        UInt FindCommonFace(const UInt edge1, const UInt edge2) const;

    private:
        // orthogonalization
        static constexpr double m_minimumEdgeLength = 1e-4;               ///< Minimum edge length
        static constexpr double m_curvilinearToOrthogonalRatio = 0.5;     ///< Ratio determining curvilinear-like(0.0) to pure(1.0) orthogonalization
        static constexpr double m_minimumCellArea = 1e-12;                ///< Minimum cell area
        static constexpr UInt m_maximumNumberOfHangingNodesAlongEdge = 5; ///< The maximum number of hanging nodes along a single element edge

        /// @brief Bounded array for storing hanging node indices.
        using HangingNodeIndexArray = std::array<UInt, m_maximumNumberOfHangingNodesAlongEdge>;

        /// @brief Compute the average area of the neighbouring faces of a triangle element
        void ComputeAverageAreOfNeighbouringFaces(const UInt faceId, UInt& numNonBoundaryFaces, double& averageOtherFacesArea) const;

        /// @brief Find the smallest angle for the three corners of the triangle
        void FindSmallestCornerAngle(const UInt faceId,
                                     double& minCosPhiSmallTriangle,
                                     UInt& nodeToPreserve,
                                     UInt& firstNodeToMerge,
                                     UInt& secondNodeToMerge,
                                     UInt& thirdEdgeSmallTriangle) const;

        /// @brief Delete the small triangle
        void DeleteSmallTriangle(const UInt nodeToPreserve,
                                 const UInt firstNodeToMerge,
                                 const UInt secondNodeToMerge,
                                 bool& nodesMerged,
                                 CompoundUndoAction& undoAction);

        /// @brief Find nodes to be deleted.
        void FindNodesToDelete(const Polygons& polygon,
                               const bool invertDeletion,
                               std::vector<bool>& isNodeInsidePolygon,
                               std::vector<bool>& deleteNode) const;

        /// @brief Delete node and edges
        void DeletedMeshNodesAndEdges(const std::function<bool(UInt)>& excludedFace,
                                      std::vector<bool>& deleteNode,
                                      CompoundUndoAction& deleteMeshAction);

        /// @brief Find nodes contained within the polygon
        std::vector<int> ComputeNodeMask(const Polygons& polygons) const;

        /// @brief Find edges that are contained within the polygon
        std::vector<int> ComputeEdgeMask(const std::vector<int>& nodeMask,
                                         bool includeIntersected) const;

        /// @brief Remove edges that are intersecting
        void RemoveIntersected(const std::vector<int>& edgeMask,
                               std::vector<int>& secondEdgeMask) const;

        /// @brief Invert the selection of edges
        void InvertSelection(const std::vector<int>& edgeMask,
                             std::vector<int>& secondEdgeMask) const;

        /// @brief Find the mesh faces that lie entirely within the polygon.
        std::vector<bool> FindFacesEntirelyInsidePolygon(const std::vector<bool>& isNodeInsidePolygon) const;

        /// @brief Deletes the mesh faces inside a polygon
        /// @param[in] polygon        The polygon where to perform the operation
        ///                           If this Polygons instance contains multiple polygons, the first one will be taken.
        /// @param[in] invertDeletion Inverts the selected node to delete (instead of outside the polygon, inside the polygon)
        [[nodiscard]] std::unique_ptr<UndoAction> DeleteMeshFaces(const Polygons& polygon, bool invertDeletion);

        /// @brief Find cells recursive, works with an arbitrary number of edges
        /// @param[in] startNode The starting node
        /// @param[in] node The current node
        /// @param[in] previousEdge The previously visited edge
        /// @param[in] numClosingEdges The number of edges closing a face (3 for triangles, 4 for quads, etc)
        /// @param[in,out] edges The vector storing the current edges forming a face
        /// @param[in,out] nodes The vector storing the current nodes forming a face
        /// @param[in,out] sortedEdges The caching array used for sorting the edges, used to inquire if an edge has been already visited
        /// @param[in,out] sortedNodes The caching array used for sorting the nodes, used to inquire if a node has been already visited
        /// @param[in,out] nodalValues The nodal values building a closed polygon
        void FindFacesRecursive(UInt startNode,
                                UInt node,
                                UInt previousEdge,
                                UInt numClosingEdges,
                                std::vector<UInt>& edges,
                                std::vector<UInt>& nodes,
                                std::vector<UInt>& sortedEdges,
                                std::vector<UInt>& sortedNodes,
                                std::vector<Point>& nodalValues);

        /// @brief Checks if a triangle has an acute angle (checktriangle)
        /// @param[in] faceNodes The face nodes composing the triangles
        /// @param[in] nodes The node coordinates
        /// @returns If triangle has an acute triangle
        [[nodiscard]] bool HasTriangleNoAcuteAngles(const std::vector<UInt>& faceNodes, const std::vector<Point>& nodes) const;

        /// @brief Determine if there are duplicate node id's on the node array
        ///
        /// The parameter sortedNodes, is a temporary array, that reduces the need to re-allocate any memory locally to this function
        bool HasDuplicateNodes(const UInt numClosingEdges, const std::vector<UInt>& node, std::vector<UInt>& sortedNodes) const;

        /// @brief Determine if there are duplicate edge-facw id's on the edges array
        ///
        /// The parameter sortedEdgesFaces, is a temporary array, that reduces the need to re-allocate any memory locally to this function
        bool HasDuplicateEdgeFaces(const UInt numClosingEdges, const std::vector<UInt>& edges, std::vector<UInt>& sortedEdgesFaces) const;

        /// @brief Resizes and initializes face vectors
        void ResizeAndInitializeFaceVectors();

        /// @brief Perform complete administration
        /// @param[in] faceNodes The input face nodes
        /// @param[in] numFaceNodes For each face, the number of nodes
        void DoAdministrationGivenFaceNodesMapping(const std::vector<std::vector<UInt>>& faceNodes,
                                                   const std::vector<std::uint8_t>& numFaceNodes);

        /// @brief Perform complete administration
        /// @param[in,out] undoAction if not null then collect any undo actions generated during the administration.
        void DoAdministration(CompoundUndoAction* undoAction = nullptr);

        /// @brief Initialise the node type array for nodes that lie on the boundary
        void InitialiseBoundaryNodeClassification(std::vector<int>& intNodeType) const;

        /// @brief Classify a single node
        MeshNodeType ClassifyNode(const UInt nodeId) const;

        // /// @brief Count the number of valid edges in list
        // UInt CountNumberOfValidEdges(const std::vector<UInt>& edgesNumFaces, const UInt numNodes) const;

        /// @brief Compute mid point and normal of polygon segment
        void ComputeMidPointsAndNormals(const std::vector<Point>& polygon,
                                        const std::vector<UInt>& edgesNumFaces,
                                        const UInt numNodes,
                                        std::array<Point, constants::geometric::maximumNumberOfNodesPerFace>& middlePoints,
                                        std::array<Point, constants::geometric::maximumNumberOfNodesPerFace>& normals,
                                        UInt& pointCount) const;

        /// @brief Compute circumcentre of face
        Point ComputeCircumCentre(const Point& centerOfMass,
                                  const UInt pointCount,
                                  const std::array<Point, constants::geometric::maximumNumberOfNodesPerFace>& middlePoints,
                                  const std::array<Point, constants::geometric::maximumNumberOfNodesPerFace>& normals) const;

        /// @brief Compute edge and average flow length
        void ComputeAverageFlowEdgesLength(std::vector<double>& edgesLength,
                                           std::vector<double>& averageFlowEdgesLength) const;

        /// @brief Compute average edge length and aspect ratios
        void ComputeAverageEdgeLength(const std::vector<double>& edgesLength,
                                      const std::vector<double>& averageFlowEdgesLength,
                                      std::vector<bool>& curvilinearGridIndicator,
                                      std::vector<std::array<double, 2>>& averageEdgesLength,
                                      std::vector<double>& aspectRatios) const;

        /// @brief Set the node indices of edges with no attached faces to invalid id.
        /// @returns The number of edges that have been invalidated
        UInt InvalidateEdgesWithNoFace();

        std::vector<MeshNodeType> m_nodesTypes; ///< The node types (nb)
    };

} // namespace meshkernel
