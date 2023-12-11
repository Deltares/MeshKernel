//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2023.
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

#include <vector>

#include "MeshKernel/BoundingBox.hpp"
#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Point.hpp"

namespace meshkernel
{

    /// @brief Forward declaration of the LandBoundary
    class LandBoundary;

    /// @brief A closed polygon.
    ///
    /// A polygon consists of at least 3 distinct points, and 1 point to close the polygon.
    class Polygon
    {
    public:
        /// @brief Default constructor.
        Polygon() = default;

        /// @brief Default copy constructor.
        Polygon(const Polygon& copy) = default;

        /// @brief Default move constructor.
        Polygon(Polygon&& copy) = default;

        /// @brief Constructor
        Polygon(const std::vector<Point>& points,
                Projection projection);

        /// @brief Constructor
        /// @note Points are moved to the polygon
        Polygon(std::vector<Point>&& points,
                Projection projection);

        ///  @brief Reset the polygon
        void Reset(const std::vector<Point>& points,
                   Projection projection);

        /// @brief Copy assignment operator
        Polygon& operator=(const Polygon& copy);

        /// @brief Move assignment operator
        Polygon& operator=(Polygon&& copy);

        /// @brief Return the number of points in the polygon
        UInt Size() const;

        /// @brief Return vector of nodes of the polygon
        const std::vector<Point>& Nodes() const;

        /// @brief Return single point at position
        Point& Node(const size_t i);

        /// @brief Return single point at position
        const Point& Node(const size_t i) const;

        /// @brief Get the bounding box of the polygon
        const BoundingBox& GetBoundingBox() const;

        /// @brief Determine if the polygon contains the point
        bool Contains(const Point& point) const;

        /// @brief Snap the section of the polygon defined by start- and end-index to the land boundary
        ///
        /// @note The bounding box may be changed
        void SnapToLandBoundary(const size_t startIndex, const size_t endIndex, const LandBoundary& boundary);

        /// @brief Refine the polygon
        /// @return The points for the refined polygon
        std::vector<Point> Refine(size_t startIndex, size_t endIndex, double refinementDistance) const;

        /// @brief Refine the polygon
        /// @return The points for the refined polygon
        std::vector<Point> LinearRefine(size_t startIndex, size_t endIndex) const;

        /// @brief Compute the area of the polygon, its centre of mass and the direction
        std::tuple<double, Point, TraversalDirection> FaceAreaAndCenterOfMass() const;

        /// @brief Compute the area of the polygon, its centre of mass and the direction
        static std::tuple<double, Point, TraversalDirection> FaceAreaAndCenterOfMass(const std::vector<Point>& polygon, const Projection projection);

        /// @brief Compute the area of the polygon, its centre of mass and the direction
        ///
        /// This version uses an indirect indexing of the set of nodes, this can be used to
        /// reduce the need to copy nodes to a separate array.
        static std::tuple<double, Point, TraversalDirection> FaceAreaAndCenterOfMass(const std::vector<Point>& nodes,
                                                                                     const std::vector<UInt>& nodeIndices,
                                                                                     const Projection projection,
                                                                                     bool isClosed);

        /// @brief Compute the perimiter length of the closed polygon
        double PerimeterLength() const;

        /// @brief Computes the edge lengths of the polygon
        /// @return edgeLengths The length of each polygon edge
        std::vector<double> EdgeLengths() const;

        /// @brief Compute the poygon offset.
        std::vector<Point> ComputeOffset(double displacement, const bool innerAndOuter) const;

        /// @brief Get the projection used.
        Projection GetProjection() const;

    private:
        /// @brief Refines the segment between two polygon nodes, starting with the node specified by the iterator up to,
        /// but not including, the next node.
        /// @param refinedPolygon     [in,out] a buffer of points into which the refined points are written
        /// @param nodeIterator       [in] position in the original, unrefined polygon that contains the first point of
        ///                           the refinement
        /// @param refinementDistance [in] the distance between two refined nodes
        /// @param projection         [in] the projection used for computing the length of the segment to be refined
        static void RefineSegment(std::vector<meshkernel::Point>& refinedPolygon,
                                  const std::vector<meshkernel::Point>::const_iterator& nodeIterator,
                                  const double refinementDistance,
                                  const meshkernel::Projection projection);

        /// @brief Compute the average length of a segment
        static void computeAverageLengths(const std::vector<double>& cumulativeDistances, std::vector<double>& averageDistances);

        /// @brief Smooth the cumulative distance from the start of the polyline
        static void smoothCumulativeDistance(const std::vector<double>& averageDistances, std::vector<double>& cumulativeDistances);

        /// @brief Smooth average length of polyline segment
        static void smoothAverageLengths(const std::vector<double>& cumulativeDistances,
                                         const double firstDistance,
                                         const double lastDistance,
                                         std::vector<double>& averageLengths);

        /// @brief Interpolate at the point on a polyline.
        static meshkernel::Point interpolatePointOnPolyline(const std::vector<meshkernel::Point>& points,
                                                            const std::vector<double>& cumulativeDistances,
                                                            const double pointDistance);

        /// @brief Check polygon has a valid state and initialise it.
        void Initialise();

        /// @brief Determine if the polygon contains the point for Cartesian coordinate system
        ///
        /// Also for spherical coordinates
        bool ContainsCartesian(const Point& point) const;

        /// @brief Determine if the polygon contains the point for accurate spherical coordinate system
        bool ContainsSphericalAccurate(const Point& point) const;

        /// @brief The point sequence making up the corners of the polygon
        std::vector<Point> m_nodes;

        /// @brief The current projection
        Projection m_projection = Projection::cartesian;

        /// @brief The bounding box containing the polygon
        BoundingBox m_boundingBox;
    };

} // namespace meshkernel

inline meshkernel::UInt meshkernel::Polygon::Size() const
{
    return static_cast<UInt>(m_nodes.size());
}

inline const std::vector<meshkernel::Point>& meshkernel::Polygon::Nodes() const
{
    return m_nodes;
}

inline meshkernel::Point& meshkernel::Polygon::Node(const size_t i)
{
    return m_nodes[i];
}

inline const meshkernel::Point& meshkernel::Polygon::Node(const size_t i) const
{
    return m_nodes[i];
}

inline const meshkernel::BoundingBox& meshkernel::Polygon::GetBoundingBox() const
{
    return m_boundingBox;
}

inline meshkernel::Projection meshkernel::Polygon::GetProjection() const
{
    return m_projection;
}
