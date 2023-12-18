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

#include <iostream>
#include <numeric>

#include "MeshKernel/BoundingBox.hpp"
#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Utilities/LinearAlgebra.hpp"
#include "MeshKernel/Utilities/RTreeBase.hpp"
#include "MeshKernel/Vector.hpp"

namespace meshkernel
{
    /// @brief Resizes and fills a two dimensional vector
    /// @tparam T The type of the vector elements
    /// @param[in] v The input two dimensional vector
    /// @param[in] firstDimension The first new dimension
    /// @param[in] secondDimension The second new dimension
    /// @param[in] fill Whatever fill or not fill the vector with missing values
    /// @param[in] fillValue The fill value
    template <typename T>
    void ResizeAndFill2DVector(std::vector<std::vector<T>>& v, UInt const& firstDimension, UInt const& secondDimension, bool fill = false, const T& fillValue = {})
    {
        v.resize(firstDimension);
        for (auto& e : v)
        {
            e.resize(secondDimension);
            if (fill)
            {
                std::fill(e.begin(), e.end(), fillValue);
            }
        }
    }

    /// @brief Resizes and fills a three dimensional vector
    /// @tparam T The type of the vector elements
    /// @param[in] v The input three dimensional vector
    /// @param[in] firstDimension The first new dimension
    /// @param[in] secondDimension The second new dimension
    /// @param[in] thirdDim The third new dimension
    /// @param[in] fill Whatever fill or not fill the vector with missing values
    /// @param[in] fillValue The fill value
    template <typename T>
    void ResizeAndFill3DVector(std::vector<std::vector<std::vector<T>>>& v, UInt const& firstDimension, UInt const& secondDimension, UInt const& thirdDim, bool fill = false, const T& fillValue = {})
    {
        v.resize(firstDimension);
        for (auto& e : v)
        {
            e.resize(secondDimension);
            for (auto& ee : e)
            {
                ee.resize(thirdDim);
                if (fill)
                {
                    std::fill(ee.begin(), ee.end(), fillValue);
                }
            }
        }
    }

    /// @brief Defines generic dot product for one dimension
    /// @tparam T Requires * operator
    /// @param[in] dx1 First component
    /// @param[in] dx2 Second component
    /// @return The dot product
    template <typename T>
    [[nodiscard]] T DotProduct(const T& dx1, const T& dx2)
    {
        return dx1 * dx2;
    }

    /// @brief Defines generic dot product of infinite dimensions
    /// @tparam T Requires * operator
    /// @tparam ...Args Parameter pack, requires a multiple of 2 number of parameters
    /// @param[in] dx1 First component
    /// @param[in] dx2 Second component
    /// @param[in] ...args Parameter pack
    /// @return The dot product
    template <typename T, typename... Args>
    [[nodiscard]] T DotProduct(const T& dx1, const T& dx2, Args&... args)
    {
        return dx1 * dx2 + DotProduct(args...);
    }

    /// @brief Defines vector product for cartesian 3D-space
    /// @param[in] a The first cartesian 3D point
    /// @param[in] b The second cartesian 3D point
    /// @return The vector product
    [[nodiscard]] Cartesian3DPoint VectorProduct(const Cartesian3DPoint& a, const Cartesian3DPoint& b);

    /// @brief Defines inner product in cartesian 3D-space
    /// @param[in] a The first cartesian 3D point
    /// @param[in] b The second cartesian 3D point
    /// @return The resulting inner product
    [[nodiscard]] double InnerProduct(const Cartesian3DPoint& a, const Cartesian3DPoint& b);

    /// @brief Find index of a certain element
    /// @param[in] vec The vector to search in
    /// @param[in] el The element to search for
    /// @returns The index of element
    template <typename T>
    [[nodiscard]] T FindIndex(const std::vector<T>& vec, T el)
    {
        for (UInt n = 0; n < vec.size(); n++)
        {
            if (vec[n] == el)
            {
                return static_cast<T>(n);
            }
        }

        return 0;
    }

    /// @brief Find all start-end positions in a vector separated by a separator
    /// @param[in] vec The vector with separator
    /// @param[in] start The start of the range to search for
    /// @param[in] end The end of the range to search for
    /// @param[in] separator The value of the separator
    /// @returns Indices of elements
    std::vector<std::pair<UInt, UInt>> FindIndices(const std::vector<Point>& vec, size_t start, size_t end, double separator);

    /// @brief Determine if the number of invalid points in the point array.
    UInt InvalidPointCount(const std::vector<Point>& points);

    /// @brief Determine if the number of invalid points in the section of the point array.
    /// @param[in] points The array of points
    /// @param[in] start The start of the range in which to search
    /// @param[in] end The end of the range in which to search
    UInt InvalidPointCount(const std::vector<Point>& points, size_t start, size_t end);

    /// @brief Sort a vector and return the sorted indices
    /// @param[in] v The vector to sort
    /// @returns The indices of elements
    template <typename T>
    [[nodiscard]] std::vector<UInt> SortedIndices(const std::vector<T>& v)
    {
        std::vector<UInt> indices(v.size());
        iota(indices.begin(), indices.end(), 0);
        std::ranges::stable_sort(indices.begin(), indices.end(), [&v](UInt i1, UInt i2)
                                 { return v[i1] < v[i2]; });
        return indices;
    }

    /// @brief Reorder vector accordingly to a specific order
    /// @param[in] v The vector to reorder
    /// @param[in] order The order to use
    /// @returns The reordered vector
    template <typename T>
    auto ReorderVector(const std::vector<T>& v, const std::vector<UInt>& order)
    {
        std::vector<T> ordered;
        ordered.reserve(v.size());
        for (const auto& value : order)
        {
            ordered.emplace_back(v[value]);
        }
        return ordered;
    }

    /// @brief Algorithm performing the zero's search using the golden section algorithm's
    /// @param[in] func Function to search for a root
    /// @param[in] min The minimum value of the interval
    /// @param[in] max The maximum value of the interval
    /// @returns The value where the function is approximately 0
    template <typename F>
    [[nodiscard]] double FindFunctionRootWithGoldenSectionSearch(F func, double min, double max)
    {
        // golden distance factors
        const double c = 0.38196602;
        const double r = 0.61803399;
        const double tolerance = 1e-5;

        double left = min;
        double middle = (min + max) * 0.5;
        double right = max;

        double x0 = left;
        double x1 = middle - c * (middle - left);
        double x2 = middle;
        double x3 = right;
        if (std::abs(right - middle) > std::abs(middle - left))
        {
            x1 = middle;
            x2 = middle + c * (right - left);
        }

        double f1 = func(x1);
        double f2 = func(x2);

        while (std::abs(x3 - x0) > tolerance * std::max(std::abs(x1) + std::abs(x2), 1e-8)) // tolerance * std::max(std::abs(x1) + std::abs(x2), 1e-10)
        {
            if (f2 < f1)
            {
                x0 = x1;
                x1 = x2;
                x2 = r * x1 + c * x3;

                f1 = f2;
                f2 = func(x2);
            }
            else
            {
                x3 = x2;
                x2 = x1;
                x1 = r * x2 + c * x0;

                f2 = f1;
                f1 = func(x1);
            }
        }

        return f1 < f2 ? x1 : x2;
    }

    /// @brief Get the next forward index.
    /// @param[in] currentIndex The current index.
    /// @param[in] size The size of the vector.
    /// @returns The next forward index.
    [[nodiscard]] UInt NextCircularForwardIndex(UInt currentIndex, UInt size);

    /// @brief Get the next backward index.
    /// @param[in] currentIndex The current index.
    /// @param[in] size The size of the vector.
    /// @returns The next backward index.
    [[nodiscard]] UInt NextCircularBackwardIndex(UInt currentIndex, UInt size);

    /// @brief Determines if a point is close to the poles (latitude close to 90 degrees).
    /// @param[in] point The current point.
    /// @returns If the point is on the pole.
    [[nodiscard]] bool IsPointOnPole(const Point& point);

    /// @brief Transforms 2D point in spherical coordinates to 3D cartesian coordinates.
    /// @param[in] sphericalPoint The current spherical point (2 coordinates).
    /// @returns The converted cartesian 3d point.
    [[nodiscard]] Cartesian3DPoint SphericalToCartesian3D(const Point& sphericalPoint);

    /// @brief Transforms 3D cartesian coordinates to 2D point in spherical coordinates
    /// @param[in] cartesianPoint The 3d cartesian point
    /// @param[in] referenceLongitude The reference longitude
    /// @returns The spherical coordinate
    [[nodiscard]] Point Cartesian3DToSpherical(const Cartesian3DPoint& cartesianPoint, double referenceLongitude);

    /// @brief Computes the cross product between two segments (duitpl)
    /// @param[in] firstSegmentFirstPoint   The first point of the first segment
    /// @param[in] firstSegmentSecondPoint  The second point of the first segment
    /// @param[in] secondSegmentFirstPoint   The second point of the second segment
    /// @param[in] secondSegmentSecondPoint The second point of the second segment
    /// @param[in] projection               The coordinate system projection
    /// @return The cross product value
    [[nodiscard]] double crossProduct(const Point& firstSegmentFirstPoint, const Point& firstSegmentSecondPoint, const Point& secondSegmentFirstPoint, const Point& secondSegmentSecondPoint, const Projection& projection);

    /// @brief Checks if a point is in polygonNodes using the winding number method
    /// @param[in] point The point to check
    /// @param[in] polygonNodes A series of closed polygons
    /// @param[in] startNode The start index in polygonNodes
    /// @param[in] endNode  The end index in polygonNodes
    /// @param[in] projection The coordinate system projection.
    /// @param[in] polygonCenter A coordinate needed in case of sphericalAccurate projection
    /// @returns If point is inside the designed polygon
    [[nodiscard]] bool IsPointInPolygonNodes(const Point& point,
                                             const std::vector<Point>& polygonNodes,
                                             const Projection& projection,
                                             Point polygonCenter = {constants::missing::doubleValue, constants::missing::doubleValue},
                                             UInt startNode = constants::missing::uintValue,
                                             UInt endNode = constants::missing::uintValue);

    /// @brief Computes three base components
    void ComputeThreeBaseComponents(const Point& point, std::array<double, 3>& exxp, std::array<double, 3>& eyyp, std::array<double, 3>& ezzp);

    /// @brief Computes two base components
    void ComputeTwoBaseComponents(const Point& point, double (&elambda)[3], double (&ephi)[3]);

    /// @brief Gets dx for the given projection
    /// @param[in] firstPoint
    /// @param[in] secondPoint
    /// @param[in] projection The coordinate system projection.
    [[nodiscard]] double GetDx(const Point& firstPoint, const Point& secondPoint, const Projection& projection);

    /// @brief Get the delta (dx, dy) for the given projection
    /// @param[in] firstPoint
    /// @param[in] secondPoint
    /// @param[in] projection The coordinate system projection.
    Vector GetDelta(const Point& firstPoint, const Point& secondPoint, const Projection& projection);

    /// @brief Get the normal to the line described by the two points.
    ///
    /// The vector is normalised.
    /// @param[in] start Point at the start of the line
    /// @param[in] end Point at the end of the line
    /// @param[in] projection The coordinate system projection.
    Vector ComputeNormalToline(const Point& start, const Point& end, const Projection projection);

    /// @brief Gets dy for the given projection
    /// @param[in] firstPoint
    /// @param[in] secondPoint
    /// @param[in] projection The coordinate system projection.
    [[nodiscard]] double GetDy(const Point& firstPoint, const Point& secondPoint, const Projection& projection);

    /// @brief Outer product of two segments (dprodout)
    [[nodiscard]] double OuterProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment,
                                                 const Point& firstPointSecondSegment, const Point& secondPointSecondSegment, const Projection& projection);

    /// @brief Computes the middle point.
    /// @param[in] firstPoint  The first point of the segment.
    /// @param[in] secondPoint The second point of the segment.
    /// @param[in] projection  The coordinate system projection.
    /// @return The middle point.
    [[nodiscard]] Point ComputeMiddlePoint(const Point& firstPoint, const Point& secondPoint, const Projection& projection);

    /// @brief Computes the middle point (account for poles, latitudes close to 90 degrees)
    /// @param[in] firstPoint  The first point of the segment.
    /// @param[in] secondPoint The second point of the segment.
    /// @param[in] projection  The coordinate system projection.
    /// @return The middle point.
    [[nodiscard]] Point ComputeMiddlePointAccountingForPoles(const Point& firstPoint, const Point& secondPoint, const Projection& projection);

    /// @brief Normalized vector of a segment in direction 1 -> 2 with the insidePoint orientation
    /// @param[in] firstPoint  The first point of the segment.
    /// @param[in] secondPoint The second point of the segment.
    /// @param[in] insidePoint The inside point of the segment
    /// @param[in] projection  The coordinate system projection.
    /// @return The Normal vector
    [[nodiscard]] Point NormalVector(const Point& firstPoint, const Point& secondPoint, const Point& insidePoint, const Projection& projection);

    /// @brief Transforms vector with components in global spherical coordinate directions(xglob, yglob)
    ///       to local coordinate directions(xloc, yloc) around reference point(xref, yref)
    void TransformGlobalVectorToLocal(const Point& reference, const Point& globalCoordinates, const Point& globalComponents, const Projection& projection, Point& localComponents);

    /// @brief Computes the normal vector outside (normalout)
    ///
    /// \see NormalVectorInside
    ///
    Point NormalVectorOutside(const Point& firstPoint, const Point& secondPoint, const Projection& projection);

    /// @brief Computes the normal vector to a line 1-2, which is *outward* w.r.t.
    ///         an 'inside' point 3.
    ///
    /// Similar to NormalVectorOutside, except that the normal vector may be flipped based on the 'inside' point.
    void NormalVectorInside(const Point& firstPoint, const Point& secondPoint, const Point& insidePoint, Point& normal, bool& flippedNormal, const Projection& projection);

    /// @brief Moves a point by adding an increment vector to it.
    /// @param[in] normal         The increment direction.
    /// @param[in] increment      The increment to use.
    /// @param[in] referencePoint The reference point containing the reference latitude to use.
    /// @param[in] projection     The coordinate system projection.
    /// @param[in,out] point The point to be incremented.
    void AddIncrementToPoint(const Point& normal, double increment, const Point& referencePoint, const Projection& projection, Point& point);

    /// @brief For a given polygon the function may shift the input coordinates
    /// @param[in,out] polygon    The input polygon.
    /// @note To be called for spherical coordinate system only
    void TranslateSphericalCoordinates(std::vector<Point>& polygon);

    /// @brief For a given polygon compute a reference point
    /// @param[in] polygon    The input polygon.
    /// @param[in]     projection The coordinate system projection.
    /// @return The reference point
    [[nodiscard]] Point ReferencePoint(const std::vector<Point>& polygon, const Projection& projection);

    /// @brief For a given polygon compute a reference point
    /// @param[in] nodes    The input nodes
    /// @param[in] polygonIndices  The polygon node indices.
    /// @param[in] projection The coordinate system projection.
    /// @return The reference point
    [[nodiscard]] Point ReferencePoint(const std::vector<Point>& nodes,
                                       const std::vector<UInt>& polygonIndices,
                                       const Projection& projection);

    /// @brief Computes the squared distance between two points
    ///        This is faster than ComputeDistance because it does not take the square root
    /// @param[in] firstPoint  The first point.
    /// @param[in] secondPoint The second point.
    /// @param[in] projection  The coordinate system projection.
    /// @return The squared distance
    [[nodiscard]] double ComputeSquaredDistance(const Point& firstPoint, const Point& secondPoint, const Projection& projection);

    /// @brief Computes the  distance between two points (dbdistance)
    /// @param[in] firstPoint  The first point.
    /// @param[in] secondPoint The second point.
    /// @param[in] projection  The coordinate system projection.
    /// @return The  distance
    [[nodiscard]] double ComputeDistance(const Point& firstPoint, const Point& secondPoint, const Projection& projection);

    /// @brief Computes the perpendicular distance of a point from a segment firstNode - secondNode (dlinedis3)
    /// @param[in] point      The point to consider in the distance calculation.
    /// @param[in] firstNode  The first point of the segment.
    /// @param[in] secondNode The second point of the segment.
    /// @param[in] projection The coordinate system projection.
    /// @return The normal distance from the segment, the intersection of the normal projection on the segment, the distance from the first node, expressed as ratio of the segment length
    [[maybe_unused]] std::tuple<double, Point, double> DistanceFromLine(const Point& point, const Point& firstNode, const Point& secondNode, const Projection& projection);

    /// @brief Inner product of two segments (dprodin)
    /// @param[in] firstPointFirstSegment   The first point of the first segment
    /// @param[in] secondPointFirstSegment  The second point of the first segment
    /// @param[in] firstPointSecondSegment  The first point of the second segment
    /// @param[in] secondPointSecondSegment The second point of the second segment
    /// @param[in] projection               The coordinate system projection
    /// @return The resulting inner product
    [[nodiscard]] double InnerProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment, const Point& firstPointSecondSegment, const Point& secondPointSecondSegment, const Projection& projection);

    /// @brief The normalized inner product of two segments (dcosphi)
    /// @param[in] firstPointFirstSegment   The first point of the first segment
    /// @param[in] secondPointFirstSegment  The second point of the first segment
    /// @param[in] firstPointSecondSegment  The first point of the second segment
    /// @param[in] secondPointSecondSegment The second point of the second segment
    /// @param[in] projection               The coordinate system projection
    /// @return The resulting normalized inner product
    [[nodiscard]] double NormalizedInnerProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment, const Point& firstPointSecondSegment, const Point& secondPointSecondSegment, const Projection& projection);

    /// @brief Computes the circumcenter of a triangle
    /// @param[in] firstNode  The first triangle node
    /// @param[in] secondNode The second triangle node
    /// @param[in] thirdNode  The third triangle node
    /// @param[in] projection The coordinate system projection
    /// @return The resulting circumcenter
    [[nodiscard]] Point CircumcenterOfTriangle(const Point& firstNode, const Point& secondNode, const Point& thirdNode, const Projection& projection);

    /// @brief Determines if two segments are crossing (cross, cross3D)
    /// @param[in]  firstSegmentFirstPoint   The first point of the first segment
    /// @param[in]  firstSegmentSecondPoint  The second point of the first segment
    /// @param[in]  secondSegmentFirstPoint   The first point of the second segment
    /// @param[in]  secondSegmentSecondPoint The second point of the second segment
    /// @param[in]  adimensionalCrossProduct Whether to compute the dimensionless cross product
    /// @param[in]  projection               The coordinate system projection
    /// @return A tuple with:
    ///  If the two segments are crossing
    ///  The intersection point
    ///  The cross product of the intersection
    ///  The distance of the intersection from the first node of the first segment, expressed as a ratio of the segment length
    ///  The distance of the intersection from the first node of the second segment, expressed as a ratio of the segment length
    [[nodiscard]] std::tuple<bool, Point, double, double, double> AreSegmentsCrossing(const Point& firstSegmentFirstPoint,
                                                                                      const Point& firstSegmentSecondPoint,
                                                                                      const Point& secondSegmentFirstPoint,
                                                                                      const Point& secondSegmentSecondPoint,
                                                                                      bool adimensionalCrossProduct,
                                                                                      const Projection& projection);

    /// @brief Computes the coordinate of a point on a spline, given the dimensionless distance from the first corner point (splint)
    /// @param[in] coordinates                 The spline node coordinates
    /// @param[in] coordinatesDerivatives      The spline nodal derivatives
    /// @param[in] pointAdimensionalCoordinate The adimensinal coordinate where to perform the interpolation
    /// @returns The interpolated point
    template <typename T>
    [[nodiscard]] T ComputePointOnSplineAtAdimensionalDistance(const std::vector<T>& coordinates,
                                                               const std::vector<T>& coordinatesDerivatives,
                                                               double pointAdimensionalCoordinate)
    {
        T pointCoordinate{};
        if (pointAdimensionalCoordinate < 0)
        {
            return pointCoordinate;
        }

        const double eps = 1e-5;
        const double splFac = 1.0;
        const auto coordinate = std::floor(pointAdimensionalCoordinate);
        if (pointAdimensionalCoordinate - coordinate < eps)
        {
            return pointCoordinate = coordinates[static_cast<UInt>(coordinate)];
        }

        const UInt low = static_cast<UInt>(coordinate);
        const UInt high = low + 1;
        const double a = high - pointAdimensionalCoordinate;
        const double b = pointAdimensionalCoordinate - low;

        pointCoordinate = coordinates[low] * a + coordinates[high] * b +
                          (coordinatesDerivatives[low] * (pow(a, 3) - a) + coordinatesDerivatives[high] * (pow(b, 3) - b)) / 6.0 * splFac;

        return pointCoordinate;
    }

    /// @brief Computes dimensionless distances of a vector of points such as the first entry has distance 0 and the last entry has distance 1.
    /// @param[in] v          The vector of points
    /// @param[in] projection The projection to use.
    /// @returns The dimensionless distances and the dimensional total distance (used for normalization)
    std::tuple<std::vector<double>, double> ComputeAdimensionalDistancesFromPointSerie(const std::vector<Point>& v, const Projection& projection);

    /// @brief Computes the sign of a type
    /// @tparam    T   A signed type
    /// @param[in] val the value to use for computing a sign
    /// @returns -1 for negatives and +1 for positives
    template <typename T>
    [[nodiscard]] int sgn(T val)
    {
        return (T(0) < val ? 1 : 0) - (val < T(0) ? 1 : 0);
    }

    /// @brief Computes the transfinite discretization inside the area defined by 4 sides, each one discretized with a series of points (tranfn2).
    /// @param[in] leftDiscretization   The first side of the area.
    /// @param[in] rightDiscretization  The second side of the area.
    /// @param[in] bottomDiscretization The third side of the area.
    /// @param[in] upperDiscretization  The fourth side of the area.
    /// @param[in] projection           The projection to use.
    /// @param[in] numM                 The number of columns to generate (horizontal direction).
    /// @param[in] numN                 The number of rows to generate (vertical direction).
    /// @returns The resulting dicretization (expressed as number of points).
    [[nodiscard]] lin_alg::Matrix<Point> DiscretizeTransfinite(const std::vector<Point>& leftDiscretization,
                                                               const std::vector<Point>& rightDiscretization,
                                                               const std::vector<Point>& bottomDiscretization,
                                                               const std::vector<Point>& upperDiscretization,
                                                               const Projection& projection,
                                                               UInt numM,
                                                               UInt numN);

    /// @brief Computes the edge centers
    /// @param[in] nodes The vector of edge nodes.
    /// @param[in] edges The vector of edge indices.
    /// @return The vector containing the edge centers.
    [[nodiscard]] std::vector<Point> ComputeEdgeCenters(const std::vector<Point>& nodes, const std::vector<Edge>& edges);

    /// @brief Given a triangles with values on each node, computes the interpolated value inside the triangle, using linear interpolation.
    /// @param[in] interpolationPoint The point where to interpolate.
    /// @param[in] polygon            The polygon containing the triangle nodes.
    /// @param[in] values             The values at each node.
    /// @param[in] projection         The projection to use.
    /// @return The interpolated value.
    [[nodiscard]] double LinearInterpolationInTriangle(const Point& interpolationPoint, const std::vector<Point>& polygon, const std::vector<double>& values, const Projection& projection);

    /// @brief Given a series of point computes the average coordinate
    /// @param[in] points The point series.
    /// @param[in] projection The projection to use.
    /// @return The average coordinate.
    [[nodiscard]] Point ComputeAverageCoordinate(const std::vector<Point>& points, const Projection& projection);

    /// @brief Cartesian projection of a point on a segment defined by other two points
    /// @param firstNode The first node of the segment
    /// @param secondNode The second node of the segment
    /// @param pointToProject The point to project
    /// @return The projected point, the distance of the projection from the first point (expressed as a ratio),
    /// a boolean indicating if the projected point is within the first and second point.
    std::tuple<Point, double, bool> OrthogonalProjectionOnSegment(Point const& firstNode,
                                                                  Point const& secondNode,
                                                                  Point const& pointToProject);

    /// @brief Calculates the absolute difference between to `Index` numbers.
    ///
    /// @param[in] number_1 The first number
    /// @param[in] number_2 The second number
    UInt AbsoluteDifference(UInt number_1, UInt number_2);

    /// @brief Computes the discretization points along a polyline
    /// @param polyline A polyline described by its nodes
    /// @param chainages The chainages used for dicretizing the current polyline
    /// @param projection The projection to use
    /// @return The discretized polyline
    [[nodiscard]] std::vector<Point> ComputePolyLineDiscretization(std::vector<Point> const& polyline, std::vector<double>& chainages, Projection projection);

    /// @brief Computes the chainages of each polyline node
    /// @param polyline A polyline described by its nodes
    /// @param projection The projection to use
    /// @return A vector containing the chainage volau of the polyline nodes
    [[nodiscard]] std::vector<double> ComputePolyLineNodalChainages(std::vector<Point> const& polyline, Projection projection);

    /// @brief Computes the lengths of each polyline segment
    /// @param polyline A polyline described by its nodes
    /// @param projection The projection to use
    /// @return A vector containing the lengths of each polyline segment
    [[nodiscard]] std::vector<double> ComputePolyLineEdgesLengths(std::vector<Point> const& polyline, Projection projection);

    /// @brief Compute the matrix norm
    /// @brief x [in] To be detailed
    /// @brief y [in] To be detailed
    /// @brief matCoefficients [in] To be detailed
    /// @returns The computed matrix norm
    [[nodiscard]] double MatrixNorm(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& matCoefficients);

    /// @brief Print the (simplified) graph in a form that can be loaded into matlab/octave.
    ///
    /// Only nodes and node connectivity need be printed to visualise the graph.
    void Print(const std::vector<Point>& nodes, const std::vector<Edge>& edges, std::ostream& out = std::cout);

    /// @brief Increment a valid value by an increment
    inline void IncrementValidValue(UInt& value, const UInt increment)
    {
        if (value != constants::missing::uintValue) [[likely]]
        {
            value += increment;
        }
    }

} // namespace meshkernel
