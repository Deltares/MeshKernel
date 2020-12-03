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

#include <algorithm>
#include <cmath>
#include <numeric>

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/SpatialTrees.hpp>

namespace meshkernel
{
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
    [[nodiscard]] Cartesian3DPoint VectorProduct(Cartesian3DPoint a, Cartesian3DPoint b);

    /// @brief Defines inner product in cartesian 3D-space
    /// @param[in] a The first cartesian 3D point
    /// @param[in] b The second cartesian 3D point
    /// @return The resulting inner product
    [[nodiscard]] double InnerProduct(Cartesian3DPoint a, Cartesian3DPoint b);

    /// @brief Find index of a certain element
    /// @param[in] vec The vector to search in
    /// @param[in] el The element to search for
    /// @returns The index of element
    template <typename T>
    [[nodiscard]] int FindIndex(const std::vector<T>& vec, T el)
    {
        int index = 0;
        for (int n = 0; n < vec.size(); n++)
        {
            if (vec[n] == el)
            {
                index = n;
                break;
            }
        }
        return index;
    }

    /// @brief Find all start-end positions in a vector separated by a separator
    /// @param[in] vec The vector with separator
    /// @param[in] start The start of the range to search for
    /// @param[in] end The end of the range to search for
    /// @param[in] separator The value of the separator
    /// @returns Indices of elements
    std::vector<std::vector<size_t>> FindIndices(const std::vector<Point>& vec,
                                                 size_t start,
                                                 size_t end,
                                                 double separator);

    /// @brief Sort a vector and return the sorted indices
    /// @param[in] v The vector to sort
    /// @returns The indices of elements
    template <typename T>
    [[nodiscard]] std::vector<int> SortedIndices(const std::vector<T>& v)
    {
        std::vector<int> indices(v.size());
        iota(indices.begin(), indices.end(), 0);
        std::stable_sort(indices.begin(), indices.end(), [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });
        return indices;
    }

    /// @brief Reorder vector accordingly to a specific order
    /// @param[in] v The vector to reorder
    /// @param[in] order The order to use
    /// @returns The reordered vector
    template <typename T>
    auto ReorderVector(const std::vector<T>& v, const std::vector<int>& order)
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
        //golden distance factors
        const double c = 0.38196602;
        const double r = 0.61803399;
        const double tolerance = 0.00001;

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
            x2 = middle + c * (middle - left);
        }

        double f1 = func(x1);
        double f2 = func(x2);

        while (std::abs(x3 - x0) > tolerance * std::max(std::abs(x1) + std::abs(x2), 1e-8))
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

    /// @brief Get the next forward index
    /// @param[in] currentIndex The current index
    /// @param[in] size The size of the vector
    /// @returns The next forward index
    [[nodiscard]] int NextCircularForwardIndex(int currentIndex, int size);

    /// @brief Get the next backward index
    /// @param[in] currentIndex The current index
    /// @param[in] size The size of the vector
    /// @returns The next backward index
    [[nodiscard]] int NextCircularBackwardIndex(int currentIndex, int size);

    /// @brief Determines if point is close to the poles (latitude close to 90 degrees)
    /// @param[in] point The current point
    /// @returns If the point is on the pole
    [[nodiscard]] bool IsPointOnPole(const Point& point);

    /// @brief Transforms 2D point in spherical coordinates to 3D cartesian coordinates
    /// @param[in] sphericalPoint The current spherical point (2 coordinates)
    /// @param[out] cartesianPoint The converted cartesian 3d point
    Cartesian3DPoint SphericalToCartesian3D(const Point& sphericalPoint);

    /// @brief Transforms 3D cartesian coordinates to 2D point in spherical coordinates
    /// @param[in] cartesianPoint
    /// @param[in] referenceLongitude
    /// @param[out] sphericalPoint
    void Cartesian3DToSpherical(const Cartesian3DPoint& cartesianPoint, double referenceLongitude, Point& sphericalPoint);

    /// @brief Tests if a point is Left|On|Right of an infinite line.
    /// @param[in] leftPoint
    /// @param[in] rightPoint
    /// @param[in] point
    /// @returns
    ///          - >0 for point left of the line through leftPoint and rightPoint
    ///          - =0 for point  on the line
    ///          - <0 for point  right of the line
    [[nodiscard]] double IsLeft(const Point& leftPoint, const Point& rightPoint, const Point& point);

    /// @brief Checks if a point is in polygonNodes using the winding number method
    /// @param[in] point
    /// @param[in] polygonNodes A closed polygonNodes consisting of a vector of numberOfPolygonPoints + 1 in counter clockwise order
    /// @param[in] startNode
    /// @param[in] endNode
    /// @param[in] projection
    /// @param[in] polygonCenter
    /// @returns If point is in polygon nodes
    [[nodiscard]] bool IsPointInPolygonNodes(Point point,
                                             const std::vector<Point>& polygonNodes,
                                             int startNode,
                                             int endNode,
                                             Projection projection,
                                             Point polygonCenter = {doubleMissingValue, doubleMissingValue});

    /// @brief Computes three base components
    void ComputeThreeBaseComponents(const Point& point, std::array<double, 3>& exxp, std::array<double, 3>& eyyp, std::array<double, 3>& ezzp);

    /// @brief Computes two base components
    void ComputeTwoBaseComponents(const Point& point, double (&elambda)[3], double (&ephi)[3]);

    /// @brief Gets dx for the given projection
    /// @param[in] firstPoint
    /// @param[in] secondPoint
    /// @param[in] projection
    [[nodiscard]] double GetDx(const Point& firstPoint, const Point& secondPoint, const Projection& projection);

    /// @brief Gets dy for the given projection
    /// @param[in] firstPoint
    /// @param[in] secondPoint
    /// @param[in] projection
    [[nodiscard]] double GetDy(const Point& firstPoint, const Point& secondPoint, const Projection& projection);

    /// @brief Outer product of two segments (dprodout)
    [[nodiscard]] double OuterProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment,
                                                 const Point& firstPointSecondSegment, const Point& secondPointSecondSegment, const Projection& projection);

    /// @brief Compute the middle point
    void MiddlePoint(const Point& firstPoint, const Point& secondPoint, Point& result, const Projection& projection);

    /// @brief Compute the middle point
    void ComputeMiddlePoint(const Point& firstPoint, const Point& secondPoint, const Projection& projection, Point& centre);

    /// @brief Normalized vector in direction 1 -> 2, in the orientation of (xu,yu)
    void NormalVector(const Point& firstPoint, const Point& secondPoint, const Point& insidePoint, Point& result, const Projection& projection);

    /// @brief Transforms vector with components in global spherical coordinate directions(xglob, yglob)
    ///       to local coordinate directions(xloc, yloc) around reference point(xref, yref)
    void TransformGlobalVectorToLocal(const Point& reference, const Point& globalCoordinates, const Point& globalComponents, Projection projection, Point& localComponents);

    /// @brief Computes the normal vector outside
    ///
    /// \see NormalVectorInside
    void NormalVectorOutside(const Point& firstPoint, const Point& secondPoint, Point& result, const Projection& projection);

    /// @brief Computes the normal vector to a line 1-2, which is *outward* w.r.t.
    ///         an 'inside' point 3.
    ///
    /// Similar to NormalVectorOutside, except that the normal vector may be flipped based on the 'inside' point.
    void NormalVectorInside(const Point& firstPoint, const Point& secondPoint, const Point& insidePoint, Point& normal, bool& flippedNormal, Projection projection);

    /// @brief Moves a point by adding a vector to it
    /// @param[in,out] point The point to be moved
    /// @param[in] normal The direction the point will be moved
    /// @param[in] increment The length of the movement
    void Add(Point& point, const Point& normal, double increment, double xf, const Projection& projection);

    void ReferencePoint(std::vector<Point>& polygon, const int numPoints, double& minX, double& minY, const Projection& projection);

    /// @brief Computes the squared distance between two points
    /// @param firstPoint The first point
    /// @param secondPoint The second point
    /// @param projection The coordinate system projection
    /// @return The squared distance
    [[nodiscard]] double ComputeSquaredDistance(Point firstPoint, Point secondPoint, Projection projection);

    /// @brief Computes the  distance between two points (dbdistance)
    /// @param firstPoint The first point
    /// @param secondPoint The second point
    /// @param projection The coordinate system projection
    /// @return The  distance
    [[nodiscard]] double ComputeDistance(Point firstPoint, Point secondPoint, Projection projection);

    /// @brief Computes the perpendicular distance of a point from a segment firstNode - secondNode (dlinedis3)
    /// @param point The point to consider in the distance calculation
    /// @param firstNode The first point of the segment
    /// @param secondNode The second point of the segment
    /// @param projection The coordinate system projection
    /// @param normalPoint The intersection of the normal projection with the segment
    /// @param ratio the distance from the first node, expressed as a ratio to the segment length
    /// @return The normal distance from the segment
    [[nodiscard]] double DistanceFromLine(Point point, Point firstNode, Point secondNode, Projection projection, Point& normalPoint, double& ratio);

    /// dprodin inner product of two segments
    [[nodiscard]] double InnerProductTwoSegments(Point firstPointFirstSegment, Point secondPointFirstSegment, Point firstPointSecondSegment, Point secondPointSecondSegment, Projection projection);

    // dcosphi
    [[nodiscard]] double NormalizedInnerProductTwoSegments(Point firstPointFirstSegment, Point secondPointFirstSegment, Point firstPointSecondSegment, Point secondPointSecondSegment, Projection projection);

    /// @brief
    /// @param p1
    /// @param p2
    /// @param p3
    /// @param projection
    /// @return
    Point CircumcenterOfTriangle(const Point& p1, const Point& p2, const Point& p3, Projection projection);

    /// @brief (cross, cross3D)
    /// @param firstSegmentFistPoint
    /// @param firstSegmentSecondPoint
    /// @param secondSegmentFistPoint
    /// @param secondSegmentSecondPoint
    /// @param adimensionalCrossProduct
    /// @param intersectionPoint
    /// @param crossProduct
    /// @param ratioFirstSegment
    /// @param ratioSecondSegment
    /// @param projection
    /// @return
    [[nodiscard]] bool AreLinesCrossing(const Point& firstSegmentFistPoint,
                                        const Point& firstSegmentSecondPoint,
                                        const Point& secondSegmentFistPoint,
                                        const Point& secondSegmentSecondPoint,
                                        bool adimensionalCrossProduct,
                                        Point& intersectionPoint,
                                        double& crossProduct,
                                        double& ratioFirstSegment,
                                        double& ratioSecondSegment,
                                        const Projection& projection);

    /// @brief Computes the area of a polygon, its center of mass, and the orientation of the edges
    /// @param polygon The input vector containing the nodes of the polygon.
    /// @param numberOfPolygonPoints The number of points to account for in the input vector.
    /// @param projection The projection to use.
    /// @param area The resulting area.
    /// @param centerOfMass The resulting center of mass.
    /// @param isCounterClockWise The orientation of the edges.
    void FaceAreaAndCenterOfMass(std::vector<Point>& polygon, size_t numberOfPolygonPoints, Projection projection, double& area, Point& centerOfMass, bool& isCounterClockWise);

    /// @brief Interpolate spline points
    /// @param[in] coordinates
    /// @param[in] coordinatesDerivatives
    /// @param[in] pointAdimensionalCoordinate
    /// @param[out] pointCoordinate
    /// @returns False if an error occurs
    template <typename T>
    [[nodiscard]] bool InterpolateSplinePoint(const std::vector<T>& coordinates,
                                              const std::vector<T>& coordinatesDerivatives,
                                              double pointAdimensionalCoordinate,
                                              T& pointCoordinate)
    {
        if (pointAdimensionalCoordinate < 0)
        {
            return false;
        }

        const double eps = 1e-5;
        const double splFac = 1.0;
        const auto intCoordinate = int(std::floor(pointAdimensionalCoordinate));
        if (pointAdimensionalCoordinate - intCoordinate < eps)
        {
            pointCoordinate = coordinates[intCoordinate];
            return true;
        }

        int low = intCoordinate;
        int high = low + 1;
        double a = high - pointAdimensionalCoordinate;
        double b = pointAdimensionalCoordinate - low;

        pointCoordinate = coordinates[low] * a + coordinates[high] * b +
                          (coordinatesDerivatives[low] * (pow(a, 3) - a) + coordinatesDerivatives[high] * (pow(b, 3) - b)) / 6.0 * splFac;

        return true;
    }

    template <class T>
    void SwapVectorElements(std::vector<T>& v, int numElements)
    {
        if (numElements > v.size())
        {
            return;
        }

        for (int i = 0; i < numElements / 2; i++)
        {
            const auto a = v[i];
            v[i] = v[i + 1];
            v[i + 1] = a;
        }
    }

    void ComputeAdimensionalDistancesFromPointSerie(const std::vector<Point>& v, Projection projection, std::vector<double>& result, double& totalDistance);

    // get the sign
    template <typename T>
    [[nodiscard]] int sgn(T val)
    {
        return (T(0) < val ? 1 : 0) - (val < T(0) ? 1 : 0);
    }

    //(DUITPL)
    [[nodiscard]] int TwoSegmentsSign(const Point& p1, const Point& p2, const Point& p3, const Point& p4, Projection projection);

    //(TRANFN2)
    std::vector<std::vector<Point>> InterpolateTransfinite(const std::vector<Point>& sideOne,
                                                           const std::vector<Point>& sideTwo,
                                                           const std::vector<Point>& sideThree,
                                                           const std::vector<Point>& sideFour,
                                                           Projection projection,
                                                           int numM,
                                                           int numN);

    [[nodiscard]] std::vector<Point> ComputeEdgeCenters(int numEdges, const std::vector<Point>& nodes, const std::vector<Edge>& edges);

    double LinearInterpolationInTriangle(Point interpolationPoint, const std::vector<Point>& polygon, const std::vector<double>& values, Projection projection);

    /// @brief Given a vector of coordinates, get the lowest upper and right points
    /// @tparam T Requires IsCoordinate<T>
    /// @param[in] values The values
    /// @param[out] lowerLeft The lower left corner
    /// @param[out] upperRight The upper right corner
    template <typename T>
    void GetBoundingBox(const std::vector<T>& values, Point& lowerLeft, Point& upperRight)
    {

        double minx = std::numeric_limits<double>::max();
        double maxx = std::numeric_limits<double>::lowest();
        double miny = std::numeric_limits<double>::max();
        double maxy = std::numeric_limits<double>::lowest();
        for (int n = 0; n < values.size(); n++)
        {
            bool isInvalid = IsEqual(values[n].x, doubleMissingValue) ||
                             IsEqual(values[n].y, doubleMissingValue);

            if (isInvalid)
            {
                continue;
            }

            minx = std::min(minx, values[n].x);
            maxx = std::max(maxx, values[n].x);
            miny = std::min(miny, values[n].y);
            maxy = std::max(maxy, values[n].y);
        }

        lowerLeft = {minx, miny};
        upperRight = {maxx, maxy};
    }

    /// @brief Checks if value is in bounding box
    /// @tparam T Requires IsCoordinate<T>
    /// @param[in] point The point to inquire
    /// @param[in] lowerLeft The lower left corner of the bounding box
    /// @param[in] upperRight The upper right corner of the bounding box
    /// @returns If the point is in the bounding box
    template <typename T>
    bool IsValueInBoundingBox(T point, Point lowerLeft, Point upperRight)
    {

        return point.x >= lowerLeft.x && point.x <= upperRight.x &&
               point.y >= lowerLeft.y && point.y <= upperRight.y;
    }

    [[nodiscard]] Point ComputeAverageCoordinate(const std::vector<Point>& points, int numPoints, Projection projection);

} // namespace meshkernel
