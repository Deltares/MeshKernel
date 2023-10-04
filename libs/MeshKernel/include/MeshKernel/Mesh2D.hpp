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

#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/Polygon.hpp>

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
    class Mesh2D : public Mesh
    {
    public:
        /// Enumerator describing the different options to delete a mesh
        enum DeleteMeshOptions
        {
            AllNodesInside = 0,
            FacesWithIncludedCircumcenters = 1,
            FacesCompletelyIncluded = 2
        };

        /// Enumerator describing the different node types
        enum class NodeTypes
        {
            internalNode,
            onRing,
            cornerNode,
            hangingNode,
            other
        };

        /// @brief Default constructor
        Mesh2D() = default;

        /// @brief Construct a mesh2d using only the projection
        /// @param[in] projection The projection to use
        Mesh2D(Projection projection);

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
               const std::vector<UInt>& numFaceNodes,
               Projection projection);

        /// @brief Create triangular grid from nodes (triangulatesamplestonetwork)
        /// @param[in] nodes Input nodes
        /// @param[in] polygons Selection polygon
        /// @param[in] projection The projection to use
        Mesh2D(const std::vector<Point>& nodes, const Polygons& polygons, Projection projection);

        /// @brief Perform complete administration
        void Administrate() override;

        /// @brief Compute face circumcenters
        void ComputeCircumcentersMassCentersAndFaceAreas(bool computeMassCenters = false);

        /// @brief Find faces: constructs the m_facesNodes mapping, face mass centers and areas (findcells)
        void FindFaces();

        /// @brief Offset the x coordinates if m_projection is spherical
        /// @param[in] minx
        /// @param[in] miny
        void OffsetSphericalCoordinates(double minx, double miny);

        /// @brief For a face create a closed polygon and fill local mapping caches (get_cellpolygon)
        /// @param[in]  faceIndex              The face index
        /// @param[out] polygonNodesCache      The node cache array filled with the nodes values
        /// @param[out] localNodeIndicesCache  The consecutive node index in polygonNodesCache (0, 1, 2,...)
        /// @param[out] globalEdgeIndicesCache The edge cache array filled with the global edge indices
        void ComputeFaceClosedPolygonWithLocalMappings(UInt faceIndex,
                                                       std::vector<Point>& polygonNodesCache,
                                                       std::vector<UInt>& localNodeIndicesCache,
                                                       std::vector<UInt>& globalEdgeIndicesCache) const;

        /// @brief For a face create a closed polygon
        /// @param[in]     faceIndex         The face index
        /// @param[in,out] polygonNodesCache The cache array to be filled with the nodes values
        void ComputeFaceClosedPolygon(UInt faceIndex, std::vector<Point>& polygonNodesCache) const;

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
        void DeleteSmallFlowEdges(double smallFlowEdgesThreshold);

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
        void DeleteSmallTrianglesAtBoundaries(double minFractionalAreaTriangles);

        /// @brief Computes m_nodesNodes, see class members
        void ComputeNodeNeighbours();

        /// @brief Get the orthogonality values, the inner product of edges and segments connecting the face circumcenters
        /// @return The edge orthogonality
        [[nodiscard]] std::vector<double> GetOrthogonality();

        /// @brief Gets the smoothness values, ratios of the face areas
        /// @return The smoothness at the edges
        [[nodiscard]] std::vector<double> GetSmoothness();

        /// @brief Gets the aspect ratios (the ratios edges lengths to flow edges lengths)
        /// @param[in,out] aspectRatios The aspect ratios (passed as reference to avoid re-allocation)
        void ComputeAspectRatios(std::vector<double>& aspectRatios);

        ///  @brief Classifies the nodes (makenetnodescoding)
        void ClassifyNodes();

        /// @brief Deletes coinciding triangles
        void DeleteDegeneratedTriangles();

        /// @brief Transform non-triangular faces in triangular faces
        void TriangulateFaces();

        /// @brief Make a dual face around the node, enlarged by a factor
        /// @param[in] node The node index
        /// @param[in] enlargementFactor The factor by which the dual face is enlarged
        /// @param[out] dualFace The dual face to be calculated
        void MakeDualFace(UInt node, double enlargementFactor, std::vector<Point>& dualFace);

        /// @brief Sorts the faces around a node, sorted in counter clock wise order
        /// @param[in] node The node index
        /// @return The face indexses
        [[nodiscard]] std::vector<UInt> SortedFacesAroundNode(UInt node) const;

        /// @brief Convert all mesh boundaries to a vector of polygon nodes, including holes (copynetboundstopol)
        /// @param[in] polygon The polygon where the operation is performed
        /// @return The resulting polygon mesh boundary
        [[nodiscard]] std::vector<Point> MeshBoundaryToPolygon(const std::vector<Point>& polygon);

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
        void DeleteHangingEdges();

        /// @brief For a collection of points, compute the face indices including them.
        /// @param[in] points The input point vector.
        /// @return The face indices including the points.
        [[nodiscard]] std::vector<UInt> PointFaceIndices(const std::vector<Point>& points);

        /// @brief Deletes a mesh in a polygon, using several options (delnet)
        /// @param[in] polygon        The polygon where to perform the operation
        ///                           If this Polygons instance contains multiple polygons, the first one will be taken.
        /// @param[in] deletionOption The deletion option
        /// @param[in] invertDeletion Inverts the selected node to delete (instead of outside the polygon, inside the polygon)
        void DeleteMesh(const Polygons& polygon, int deletionOption, bool invertDeletion);

        /// @brief Inquire if a segment is crossing a face
        /// @param[in] firstPoint The first point of the segment
        /// @param[in] secondPoint The second point of the segment
        /// @return A tuple with the intersectedFace face index and intersected  edge index
        [[nodiscard]] std::tuple<UInt, UInt> IsSegmentCrossingABoundaryEdge(const Point& firstPoint, const Point& secondPoint) const;

        /// @brief Gets the edges and faces intersected by a polyline, with additional information on the intersections
        /// @param[in] polyLine An input polyline, defined as a series of points
        /// @return A tuple containing a vector of EdgeMeshPolylineIntersections and FaceMeshPolylineIntersections
        [[nodiscard]] std::tuple<std::vector<EdgeMeshPolylineIntersection>,
                                 std::vector<FaceMeshPolylineIntersection>>
        GetPolylineIntersections(const std::vector<Point>& polyLine);

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

        UInt m_maxNumNeighbours = 0; ///< Maximum number of neighbours

        /// @brief Merges mesh connectivity.
        ///
        /// Only merges the mesh connectivity graphs and updates indices.
        /// @note Does not do any administration on the node, edges or elements,
        /// it may be required to call Administrate after merging
        static Mesh2D Merge(const Mesh2D& mesh1, const Mesh2D& mesh2);

    private:
        // orthogonalization
        static constexpr double m_minimumEdgeLength = 1e-4;               ///< Minimum edge length
        static constexpr double m_curvilinearToOrthogonalRatio = 0.5;     ///< Ratio determining curvilinear-like(0.0) to pure(1.0) orthogonalization
        static constexpr double m_minimumCellArea = 1e-12;                ///< Minimum cell area
        static constexpr double m_weightCircumCenter = 1.0;               ///< Weight circum center
        static constexpr UInt m_maximumNumberOfHangingNodesAlongEdge = 5; ///< The maximum number of hanging nodes along a single element edge

        /// @brief Bounded array for storing hanging node indices.
        using HangingNodeIndexArray = std::array<UInt, m_maximumNumberOfHangingNodesAlongEdge>;

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

        /// @brief Resizes and initializes face vectors
        void ResizeAndInitializeFaceVectors()
        {
            // face administration
            m_edgesNumFaces.resize(m_edges.size());
            std::ranges::fill(m_edgesNumFaces, 0);

            m_edgesFaces.resize(m_edges.size());
            std::ranges::fill(m_edgesFaces, std::array{constants::missing::uintValue, constants::missing::uintValue});

            m_facesMassCenters.clear();
            m_faceArea.clear();
            m_facesNodes.clear();
            m_facesEdges.clear();
            m_facesCircumcenters.clear();
            m_numFacesNodes.clear();

            m_facesMassCenters.reserve(GetNumNodes());
            m_faceArea.reserve(GetNumNodes());
            m_facesNodes.reserve(GetNumNodes());
            m_facesEdges.reserve(GetNumNodes());
            m_facesCircumcenters.reserve(GetNumNodes());
            m_numFacesNodes.reserve(GetNumNodes());
        }
    };

} // namespace meshkernel
