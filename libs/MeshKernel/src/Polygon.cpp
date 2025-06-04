//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2025.
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

#include <cmath>
#include <limits>

#include "MeshKernel/Cartesian3DPoint.hpp"
#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/LandBoundary.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/Polygon.hpp"
#include "MeshKernel/Vector.hpp"

meshkernel::Polygon::Polygon(const std::vector<Point>& points,
                             Projection projection) : m_nodes(points), m_projection(projection)
{
    Initialise();
}

meshkernel::Polygon::Polygon(std::vector<Point>&& points,
                             Projection projection) : m_nodes(points), m_projection(projection)
{
    Initialise();
}

void meshkernel::Polygon::Initialise()
{
    if (0 < m_nodes.size() && m_nodes.size() < constants::geometric::numNodesInTriangle + 1)
    {
        throw ConstraintError("Insufficient nodes in the polygon: {}, require at least 3 (+1, making 4, to close)",
                              m_nodes.size());
    }

    if (m_nodes.size() > 0 && m_nodes[0] != m_nodes[m_nodes.size() - 1])
    {
        throw ConstraintError("Polygon is not closed");
    }

    if (InvalidPointCount(m_nodes) > 0)
    {
        throw ConstraintError("Polygon nodes contains invalid nodes");
    }

    m_boundingBox.Reset(m_nodes);
}

meshkernel::Polygon& meshkernel::Polygon::operator=(const Polygon& copy)
{
    if (this != &copy)
    {
        m_nodes = copy.m_nodes;
        m_projection = copy.m_projection;
        m_boundingBox = copy.m_boundingBox;
    }

    return *this;
}

meshkernel::Polygon& meshkernel::Polygon::operator=(Polygon&& copy)
{
    if (this != &copy)
    {
        m_nodes = std::move(copy.m_nodes);
        m_projection = copy.m_projection;
        m_boundingBox = copy.m_boundingBox;
    }

    return *this;
}

void meshkernel::Polygon::Reset(const std::vector<Point>& points,
                                Projection projection)
{
    m_projection = projection;
    m_nodes = points;
    Initialise();
}

bool meshkernel::Polygon::ContainsCartesian(const Point& point) const
{
    int windingNumber = 0;

    for (size_t n = 0; n < m_nodes.size() - 1; n++)
    {
        // TODO always Cartesian
        // So Dx and Dy can be simplified (no branching)
        // Then for 2 or more points, return multiple cross product values
        const auto crossProductValue = crossProduct(m_nodes[n], m_nodes[n + 1], m_nodes[n], point, Projection::cartesian);

        if (IsEqual(crossProductValue, 0.0))
        {
            // check if is on the line or outside
            const double deltaXSegment = GetDeltaXCartesian(m_nodes[n], m_nodes[n + 1]);
            const double lambdaX = std::abs(deltaXSegment) > 0.0 ? GetDeltaXCartesian(m_nodes[n], point) / deltaXSegment : 0.0;

            const double deltaYSegment = GetDeltaYCartesian(m_nodes[n], m_nodes[n + 1]);
            const double lambdaY = std::abs(deltaYSegment) > 0.0 ? GetDeltaYCartesian(m_nodes[n], point) / deltaYSegment : 0.0;

            if (lambdaX >= 0.0 && lambdaX <= 1.0 && lambdaY >= 0.0 && lambdaY <= 1.0)
            {
                return true;
            }
        }

        if (m_nodes[n].y <= point.y) // an upward crossing
        {
            if (m_nodes[n + 1].y > point.y && crossProductValue > 0.0)

            {
                ++windingNumber; // have  a valid up intersect
            }
        }
        else
        {
            if (m_nodes[n + 1].y <= point.y && crossProductValue < 0.0) // a downward crossing
            {
                --windingNumber; // have  a valid down intersect
            }
        }
    }

    // If winding number is not zero then the point is contained within the polygon
    return windingNumber != 0;
}

bool meshkernel::Polygon::ContainsSphericalAccurate(const Point& point) const
{
    // get 3D polygon coordinates
    std::vector<Cartesian3DPoint> cartesian3DPoints;
    cartesian3DPoints.reserve(Size());

    for (UInt i = 0; i < m_nodes.size(); ++i)
    {
        cartesian3DPoints.emplace_back(SphericalToCartesian3D(m_nodes[i]));
    }

    // enlarge around polygon
    const double enlargementFactor = 1.000001;

    // TODO set to centre?
    Point polygonCenter;
    const Cartesian3DPoint polygonCenterCartesian3D{SphericalToCartesian3D(polygonCenter)};
    for (UInt i = 0; i < m_nodes.size(); ++i)
    {
        cartesian3DPoints[i].x = polygonCenterCartesian3D.x + enlargementFactor * (cartesian3DPoints[i].x - polygonCenterCartesian3D.x);
        cartesian3DPoints[i].y = polygonCenterCartesian3D.y + enlargementFactor * (cartesian3DPoints[i].y - polygonCenterCartesian3D.y);
        cartesian3DPoints[i].z = polygonCenterCartesian3D.z + enlargementFactor * (cartesian3DPoints[i].z - polygonCenterCartesian3D.z);
    }

    // convert point
    const Cartesian3DPoint pointCartesian3D{SphericalToCartesian3D(point)};

    // get test direction: e_lambda
    const double lambda = point.x * constants::conversion::degToRad;
    const Cartesian3DPoint ee{-std::sin(lambda), std::cos(lambda), 0.0};
    int inside = 0;

    // loop over the polygon nodes
    for (UInt i = 0; i < m_nodes.size() - 1; ++i)
    {
        const auto nextNode = NextCircularForwardIndex(i, static_cast<UInt>(m_nodes.size()));
        const auto xiXxip1 = VectorProduct(cartesian3DPoints[i], cartesian3DPoints[nextNode]);
        const auto xpXe = VectorProduct(pointCartesian3D, ee);

        const double D = InnerProduct(xiXxip1, ee);
        double zeta = 0.0;
        double xi = 0.0;
        double eta = 0.0;

        if (std::abs(D) > 0.0)
        {
            xi = -InnerProduct(xpXe, cartesian3DPoints[nextNode]) / D;
            eta = InnerProduct(xpXe, cartesian3DPoints[i]) / D;
            zeta = -InnerProduct(xiXxip1, pointCartesian3D) / D;
        }

        if (IsEqual(zeta, 0.0))
        {
            return true;
        }

        if (xi >= 0.0 && eta >= 0.0 && zeta >= 0.0)
        {
            inside = 1 - inside;
        }
    }

    return inside == 1;
}

bool meshkernel::Polygon::Contains(const Point& pnt) const
{
    if (!pnt.IsValid())
    {
        throw ConstraintError("Point is not valid");
    }

    if (m_nodes.empty())
    {
        return true;
    }

    if (m_nodes.size() < constants::geometric::numNodesInTriangle)
    {
        return false;
    }

    if (!m_boundingBox.Contains(pnt))
    {
        return false;
    }

    if (m_projection == Projection::cartesian || m_projection == Projection::spherical)
    {
        return ContainsCartesian(pnt);
    }
    // projection = Projection::sphericalAccurate
    return ContainsSphericalAccurate(pnt);
}

// TODO does this need start and end points, probably not all the polygon is to be snapped
void meshkernel::Polygon::SnapToLandBoundary(const size_t startIndex, const size_t endIndex, const LandBoundary& landBoundary)
{

    if (startIndex >= m_nodes.size())
    {
        throw ConstraintError("The start index is not valid: {}.", startIndex);
    }

    if (endIndex >= m_nodes.size())
    {
        throw ConstraintError("The end index is not valid: {}.", endIndex);
    }

    const auto numNodes = Size();
    for (auto i = 0u; i < numNodes; i++)
    {
        // Adjust currentIndex for wrap-around
        const auto currentIndex = (startIndex + i) % numNodes;

        if (m_nodes[currentIndex].IsValid())
        {
            m_nodes[currentIndex] = landBoundary.FindNearestPoint(m_nodes[currentIndex], m_projection);
        }

        // Break the loop if we have reached the lastIndex
        if (currentIndex == endIndex)
        {
            break;
        }
    }

    if (m_projection == Projection::spherical)
    {
        // TODO Should this be called for spherical accurate too?
        // TODO WHere did I find this? Is it correct?
        TranslateSphericalCoordinates(m_nodes);
    }

    // Now update the bounding box
    m_boundingBox.Reset(m_nodes);
}

std::vector<double> meshkernel::Polygon::EdgeLengths() const
{
    std::vector<double> edgeLengths;
    edgeLengths.reserve(m_nodes.size());

    for (size_t p = 0; p < m_nodes.size(); ++p)
    {
        size_t firstNode = p;
        size_t secondNode = p + 1;

        if (secondNode == m_nodes.size())
        {
            secondNode = 0;
        }

        edgeLengths.emplace_back(ComputeDistance(m_nodes[firstNode], m_nodes[secondNode], m_projection));
    }

    return edgeLengths;
}

void meshkernel::Polygon::RefineSegment(std::vector<meshkernel::Point>& refinedPolygon,
                                        const std::vector<meshkernel::Point>::const_iterator& nodeIterator,
                                        const double refinementDistance,
                                        const meshkernel::Projection projection)
{
    // Line segment starting point.
    const auto& n0 = *nodeIterator;
    const auto& n1 = *std::next(nodeIterator);

    refinedPolygon.push_back(n0);

    const double segmentLength = ComputeDistance(n0, n1, projection);
    long int n = std::lround(segmentLength / refinementDistance);

    for (long int i = 1; i < n; ++i)
    {
        double lambda = static_cast<double>(i) / static_cast<double>(n);
        refinedPolygon.push_back((1.0 - lambda) * n0 + lambda * n1);
    }
}

void meshkernel::Polygon::computeAverageLengths(const std::vector<double>& cumulativeDistances, std::vector<double>& averageDistances)
{
    averageDistances[0] = cumulativeDistances[1] - cumulativeDistances.front();
    averageDistances.back() = cumulativeDistances.back() - cumulativeDistances[cumulativeDistances.size() - 2];
    for (meshkernel::UInt i = 1; i < averageDistances.size() - 1; ++i)
    {
        averageDistances[i] = 0.5 * (cumulativeDistances[i + 1] - cumulativeDistances[i - 1]);
    }
}

void meshkernel::Polygon::smoothCumulativeDistance(const std::vector<double>& averageDistances, std::vector<double>& cumulativeDistances)
{
    double disc = std::accumulate(averageDistances.begin(),
                                  averageDistances.end(), 0.0) -
                  0.5 * (averageDistances.front() + averageDistances.back());
    double dfac = cumulativeDistances.back() / disc;

    double cumDistance = 0.0;
    for (meshkernel::UInt i = 1; i < cumulativeDistances.size() - 1; ++i)
    {
        cumDistance += 0.5 * (averageDistances[i - 1] + averageDistances[i]);
        cumulativeDistances[i] = cumDistance * dfac;
    }
}

void meshkernel::Polygon::smoothAverageLengths(const std::vector<double>& cumulativeDistances,
                                               const double firstDistance,
                                               const double lastDistance,
                                               std::vector<double>& averageLengths)
{
    if (cumulativeDistances.size() <= 1)
    {
        return;
    }

    for (meshkernel::UInt i = 0; i < cumulativeDistances.size(); ++i)
    {
        const double factor = cumulativeDistances[i] / cumulativeDistances.back();
        averageLengths[i] = (1.0 - factor) * firstDistance + factor * lastDistance;
    }
}

meshkernel::Point meshkernel::Polygon::interpolatePointOnPolyline(const std::vector<meshkernel::Point>& points,
                                                                  const std::vector<double>& cumulativeDistances,
                                                                  const double pointDistance)
{
    if (cumulativeDistances.empty())
    {
        throw meshkernel::ConstraintError("Distances vector is empty!");
    }

    meshkernel::UInt intervalIndex = meshkernel::constants::missing::uintValue;
    for (meshkernel::UInt i = 0; i < cumulativeDistances.size(); ++i)
    {
        if (cumulativeDistances[i] > pointDistance)
        {
            intervalIndex = i;
            break;
        }
    }

    if (intervalIndex == meshkernel::constants::missing::uintValue)
    {
        intervalIndex = static_cast<meshkernel::UInt>(cumulativeDistances.size()) - static_cast<meshkernel::UInt>(1);
    }

    const double dt = cumulativeDistances[intervalIndex] - cumulativeDistances[intervalIndex - 1];
    double ti = 0.0;

    if (dt != 0.0)
    {
        ti = (pointDistance - cumulativeDistances[intervalIndex - 1]) / dt;
    }

    return (1.0 - ti) * points[intervalIndex - 1] + ti * points[intervalIndex];
}

void meshkernel::Polygon::ComputeResampledNodes(const size_t numberOfNewNodes, const std::vector<double>& segmentLengths, const std::vector<size_t>& nodeIndices, std::vector<Point>& refinedPolygon) const
{
    double delta = segmentLengths.back() / static_cast<double>(numberOfNewNodes);
    double distanceAlongPolygon = 0.0;
    size_t segmentEndIndex = 0;

    for (size_t i = 1; i < numberOfNewNodes; ++i)
    {
        distanceAlongPolygon += delta;

        while (distanceAlongPolygon > segmentLengths[segmentEndIndex])
        {
            if (IsEqual(distanceAlongPolygon, segmentLengths[segmentEndIndex]))
            {
                break;
            }

            ++segmentEndIndex;

            if (segmentEndIndex == segmentLengths.size())
            {
                throw meshkernel::ConstraintError("Inconsistency between distanceAlongPolygon and segmentLengths");
            }
        }

        double lambda = (distanceAlongPolygon - segmentLengths[segmentEndIndex - 1]) / (segmentLengths[segmentEndIndex] - segmentLengths[segmentEndIndex - 1]);

        size_t nodeStartIndex = nodeIndices[segmentEndIndex - 1];
        size_t nodeEndIndex = nodeIndices[segmentEndIndex];

        Point newPoint = (1.0 - lambda) * m_nodes[nodeStartIndex] + lambda * m_nodes[nodeEndIndex];

        refinedPolygon.push_back(newPoint);
    }
}

std::vector<meshkernel::Point> meshkernel::Polygon::Refine(const UInt startIndex, const UInt endIndex, const double refinementDistance) const
{
    if (startIndex == endIndex)
    {
        return m_nodes;
    }

    if (startIndex >= m_nodes.size() || endIndex >= m_nodes.size())
    {
        throw ConstraintError("The indices are not valid: {}, {}.", startIndex, endIndex);
    }

    std::vector<Point> refinedPolygon;
    const auto iStart = static_cast<std::vector<Point>::difference_type>(startIndex);
    const auto iEnd = static_cast<std::vector<Point>::difference_type>(endIndex);
    const auto from = std::next(m_nodes.begin(), iStart);
    const auto to = std::next(m_nodes.begin(), iEnd);

    if (startIndex < endIndex)
    {
        std::vector<double> segmentLengths(endIndex - startIndex + 1);
        std::vector<size_t> nodeIndices(endIndex + m_nodes.size() - startIndex);
        segmentLengths[0] = 0.0;
        nodeIndices[0] = startIndex;

        std::iota(nodeIndices.begin(), nodeIndices.end(), startIndex);

        for (size_t i = startIndex + 1; i <= endIndex; ++i)
        {
            segmentLengths[i - startIndex] = segmentLengths[i - 1 - startIndex] + ComputeDistance(m_nodes[i - 1], m_nodes[i], m_projection);
        }

        const double length = segmentLengths.back();
        const size_t numberOfNewNodes = static_cast<size_t>(std::ceil(length / refinementDistance));
        const size_t sizeEstimate = m_nodes.size() + numberOfNewNodes;

        refinedPolygon.reserve(sizeEstimate);

        // copy the first segments, up to startIndex
        std::copy(m_nodes.begin(), from + 1, std::back_inserter(refinedPolygon));

        ComputeResampledNodes(numberOfNewNodes, segmentLengths, nodeIndices, refinedPolygon);

        // copy the nodes starting with endIndex until the end
        std::copy(to, m_nodes.end(), std::back_inserter(refinedPolygon));
    }
    else
    {
        std::vector<double> segmentLengths(endIndex + m_nodes.size() - startIndex);
        std::vector<size_t> nodeIndices(endIndex + m_nodes.size() - startIndex);

        segmentLengths[0] = 0.0;
        nodeIndices[0] = startIndex;
        size_t count = 1;

        Point segmentStartPoint = m_nodes[startIndex];

        for (size_t i = startIndex + 1; i < m_nodes.size() - 1; ++i)
        {
            segmentLengths[count] = segmentLengths[count - 1] + ComputeDistance(segmentStartPoint, m_nodes[i], m_projection);
            segmentStartPoint = m_nodes[i];
            nodeIndices[count] = i;
            ++count;
        }

        for (size_t i = 0; i <= endIndex; ++i)
        {
            segmentLengths[count] = segmentLengths[count - 1] + ComputeDistance(segmentStartPoint, m_nodes[i], m_projection);
            segmentStartPoint = m_nodes[i];
            nodeIndices[count] = i;
            ++count;
        }

        const double length = segmentLengths[segmentLengths.size() - 1];

        const size_t numberOfNewNodes = static_cast<size_t>(std::ceil(length / refinementDistance));
        const size_t sizeEstimate = m_nodes.size() + numberOfNewNodes;

        refinedPolygon.reserve(sizeEstimate);

        // copy the first segments, up to startIndex
        std::copy(to, from + 1, std::back_inserter(refinedPolygon));

        ComputeResampledNodes(numberOfNewNodes, segmentLengths, nodeIndices, refinedPolygon);

        // Close the polygon
        refinedPolygon.push_back(*to);
    }

    return refinedPolygon;
}

void meshkernel::Polygon::GetPolygonNodes(const UInt startIndex,
                                          const UInt endIndex,
                                          std::vector<Point>& polygonNodes) const
{
    if (startIndex < endIndex)
    {
        for (UInt i = 0; i < polygonNodes.size(); ++i)
        {
            auto polygonNodeIndex = i + startIndex;
            if (polygonNodeIndex >= m_nodes.size())
            {
                polygonNodeIndex = polygonNodeIndex - static_cast<UInt>(m_nodes.size());
            }
            polygonNodes[i] = m_nodes[polygonNodeIndex];
        }
    }
    else
    {
        UInt count = 0;

        for (UInt i = startIndex; i < m_nodes.size(); ++i)
        {
            polygonNodes[count] = m_nodes[i];
            ++count;
        }

        // Do not include the start/end point twice.
        for (UInt i = 1; i <= endIndex; ++i)
        {
            polygonNodes[count] = m_nodes[i];
            ++count;
        }
    }
}

std::vector<double> meshkernel::Polygon::ComputeCumulativeDistances(const std::vector<meshkernel::Point>& polygonNodes) const
{
    std::vector<double> cumulativeDistances(polygonNodes.size(), 0.0);
    cumulativeDistances[0] = 0.0;

    for (UInt i = 1; i < polygonNodes.size(); ++i)
    {
        cumulativeDistances[i] = cumulativeDistances[i - 1] + ComputeDistance(polygonNodes[i], polygonNodes[i - 1], m_projection);
    }

    return cumulativeDistances;
}

std::tuple<meshkernel::UInt, meshkernel::UInt> meshkernel::Polygon::FindMinMaxRatioIndex(const std::vector<double>& averageLengths,
                                                                                         const std::vector<double>& actualAverageLengths) const
{
    double minRatio = std::numeric_limits<double>::max();
    UInt minRatioIndex = constants::missing::uintValue;
    double maxRatio = std::numeric_limits<double>::lowest();
    UInt maxRatioIndex = constants::missing::uintValue;

    for (UInt i = 0; i < averageLengths.size() - 1; ++i)
    {
        const double currentRatio = actualAverageLengths[i] / averageLengths[i];

        if (i > 0 && currentRatio < minRatio)
        {
            minRatioIndex = i;
            minRatio = currentRatio;
        }

        if (currentRatio > maxRatio)
        {
            maxRatioIndex = i;
            maxRatio = currentRatio;
        }
    }

    return {minRatioIndex, maxRatioIndex};
}

std::vector<meshkernel::Point> meshkernel::Polygon::LinearRefine(const UInt startIndex, const UInt endIndex) const
{

    if (startIndex >= m_nodes.size())
    {
        throw ConstraintError("The start index is greater than the number of points in the polygon: {} >= {}.",
                              startIndex,
                              m_nodes.size());
    }

    if (endIndex >= m_nodes.size())
    {
        throw ConstraintError("The end index is greater than the number of points in the polygon: {} >= {}.",
                              endIndex,
                              m_nodes.size());
    }

    if (startIndex == endIndex)
    {
        return m_nodes;
    }

    if (m_nodes.size() < 4)
    {
        return m_nodes;
    }

    UInt numPolygonNodes;
    if (endIndex > startIndex)
    {
        numPolygonNodes = static_cast<UInt>(endIndex) - static_cast<UInt>(startIndex) + 1u;
    }
    else
    {
        numPolygonNodes = static_cast<UInt>(m_nodes.size()) - static_cast<UInt>(startIndex) + static_cast<UInt>(endIndex);
    }

    std::vector<Point> polygonNodes(numPolygonNodes);
    std::vector<double> averageLengths(numPolygonNodes, 0.0);
    std::vector<double> actualAverageLengths(numPolygonNodes, 0.0);
    std::vector<Point> result(numPolygonNodes);

    GetPolygonNodes(startIndex, endIndex, polygonNodes);

    std::vector<double> cumulativeDistances(ComputeCumulativeDistances(polygonNodes));
    std::vector<double> initialCumulativeDistances(cumulativeDistances);

    computeAverageLengths(cumulativeDistances, averageLengths);
    const double firstLength = averageLengths.front();
    const double lastLength = averageLengths.back();
    double cumulativeDistanceTarget = 0.0;
    const UInt maxInnerIter = 20;

    // Get an estimate of number of nodes to be added
    const UInt estimateOfNodesToBeAdded = static_cast<UInt>(2.0 * cumulativeDistances.back() / (firstLength + lastLength));
    const UInt maxOuterIter = static_cast<UInt>(m_nodes.size()) + static_cast<UInt>(5) * estimateOfNodesToBeAdded;

    UInt outerIter = 0;

    while (outerIter < maxOuterIter)
    {
        for (UInt iter = 1; iter <= maxInnerIter; ++iter)
        {
            for (UInt i = 0; i < cumulativeDistances.size(); ++i)
            {
                result[i] = interpolatePointOnPolyline(polygonNodes, initialCumulativeDistances, cumulativeDistances[i]);
            }

            computeAverageLengths(cumulativeDistances, actualAverageLengths);
            smoothAverageLengths(cumulativeDistances, firstLength, lastLength, averageLengths);

            smoothCumulativeDistance(averageLengths, cumulativeDistances);
            cumulativeDistanceTarget = std::accumulate(averageLengths.begin(), averageLengths.end(), 0.0) - 0.5 * (averageLengths.front() +
                                                                                                                   averageLengths.back());
        }

        auto [minRatioIndex, maxRatioIndex] = FindMinMaxRatioIndex(averageLengths, actualAverageLengths);

        if (minRatioIndex != constants::missing::uintValue && cumulativeDistanceTarget - 1.5 * averageLengths[minRatioIndex] > initialCumulativeDistances.back())
        {
            --numPolygonNodes;

            for (UInt i = minRatioIndex; i < numPolygonNodes; ++i)
            {
                cumulativeDistances[i] = cumulativeDistances[i + 1];
            }

            averageLengths.resize(numPolygonNodes);
            actualAverageLengths.resize(numPolygonNodes);
            cumulativeDistances.resize(numPolygonNodes);
            result.resize(numPolygonNodes);
        }
        else if (cumulativeDistanceTarget + 0.5 * actualAverageLengths[maxRatioIndex] < initialCumulativeDistances.back())
        {
            ++numPolygonNodes;
            averageLengths.resize(numPolygonNodes);
            actualAverageLengths.resize(numPolygonNodes);
            cumulativeDistances.resize(numPolygonNodes);
            result.resize(numPolygonNodes);

            const UInt lowerLimit = (maxRatioIndex == constants::missing::uintValue ? 0u : maxRatioIndex) + 2u;

            for (UInt i = static_cast<UInt>(cumulativeDistances.size()) - 1u; i >= lowerLimit; --i)
            {
                cumulativeDistances[i] = cumulativeDistances[i - 1];
            }

            cumulativeDistances[maxRatioIndex + 1] = 0.5 * (cumulativeDistances[maxRatioIndex] + cumulativeDistances[maxRatioIndex + 2]);
        }
        else
        {
            break;
        }

        outerIter++;
    }

    if (endIndex > startIndex)
    {
        result.insert(result.begin(), m_nodes.begin(), m_nodes.begin() + startIndex);
        result.insert(result.end(), m_nodes.begin() + endIndex + 1u, m_nodes.end());
    }
    else
    {
        result.insert(result.begin(), m_nodes.begin() + endIndex, m_nodes.begin() + startIndex);
    }

    return result;
}

std::tuple<double, meshkernel::Point, meshkernel::TraversalDirection> meshkernel::Polygon::FaceAreaAndCenterOfMass(const std::vector<Point>& polygon, const Projection projection)
{

    if (polygon.size() < constants::geometric::numNodesInTriangle)
    {
        throw std::invalid_argument("FaceAreaAndCenterOfMass: The polygon has less than 3 unique nodes.");
    }

    double area = 0.0;

    const double minArea = 1e-8;
    const auto numberOfPointsOpenedPolygon = static_cast<UInt>(polygon.size()) - 1;

    const double updateStepSize = 0.1;

    Point centreOfMass(0.0, 0.0);

    for (UInt n = 0; n < numberOfPointsOpenedPolygon; ++n)
    {
        centreOfMass += polygon[n];
    }

    centreOfMass *= 1.0 / static_cast<double>(numberOfPointsOpenedPolygon);
    // Will be non-unity for spherical coordinates only
    const double xTransformation = projection == Projection::cartesian ? 1.0 : 1.0 / std::cos(centreOfMass.y * constants::conversion::degToRad);
    const double circumcentreTolerance = constants::geometric::circumcentreTolerance * (projection == Projection::cartesian ? 1.0 : 1.0 / (constants::geometric::earth_radius * constants::conversion::degToRad));

    if (numberOfPointsOpenedPolygon == constants::geometric::numNodesInTriangle)
    {

        Point midPoint1 = 0.5 * (polygon[0] + polygon[1]);
        Point midPoint2 = 0.5 * (polygon[1] + polygon[2]);
        Point midPoint3 = 0.5 * (polygon[2] + polygon[0]);

        Vector edgeVector1 = static_cast<Vector>(NormalVector(polygon[0], polygon[1], midPoint1, projection));
        Vector edgeVector2 = static_cast<Vector>(NormalVector(polygon[1], polygon[2], midPoint2, projection));
        Vector edgeVector3 = static_cast<Vector>(NormalVector(polygon[2], polygon[0], midPoint3, projection));

        edgeVector1.normalise();
        edgeVector2.normalise();
        edgeVector3.normalise();

        Vector edgeVectorSum = edgeVector1 + edgeVector2 + edgeVector3;

        double edgeVectorSumLength = edgeVectorSum.length();

        for (UInt i = 1; i <= constants::numeric::MaximumNumberOfCircumcentreIterations; ++i)
        {
            Vector delta1 = GetDelta(midPoint1, centreOfMass, projection);
            Vector delta2 = GetDelta(midPoint2, centreOfMass, projection);
            Vector delta3 = GetDelta(midPoint3, centreOfMass, projection);

            double ds = dot(delta1, edgeVector1) + dot(delta2, edgeVector2) + dot(delta3, edgeVector3);

            if (projection != Projection::cartesian)
            {
                ds *= constants::conversion::radToDeg * constants::geometric::inverse_earth_radius;
            }

            centreOfMass.x -= updateStepSize * ds * edgeVectorSum.x() * xTransformation;
            centreOfMass.y -= updateStepSize * ds * edgeVectorSum.y();

            if (ds * edgeVectorSumLength < circumcentreTolerance || i == constants::numeric::MaximumNumberOfCircumcentreIterations)
            {
                break;
            }
        }
    }
    else if (numberOfPointsOpenedPolygon == constants::geometric::numNodesInQuadrilateral)
    {

        Point midPoint1 = 0.5 * (polygon[0] + polygon[1]);
        Point midPoint2 = 0.5 * (polygon[1] + polygon[2]);
        Point midPoint3 = 0.5 * (polygon[2] + polygon[3]);
        Point midPoint4 = 0.5 * (polygon[3] + polygon[0]);

        Vector edgeVector1 = static_cast<Vector>(NormalVector(polygon[0], polygon[1], midPoint1, projection));
        Vector edgeVector2 = static_cast<Vector>(NormalVector(polygon[1], polygon[2], midPoint2, projection));
        Vector edgeVector3 = static_cast<Vector>(NormalVector(polygon[2], polygon[3], midPoint3, projection));
        Vector edgeVector4 = static_cast<Vector>(NormalVector(polygon[3], polygon[0], midPoint4, projection));

        edgeVector1.normalise();
        edgeVector2.normalise();
        edgeVector3.normalise();
        edgeVector4.normalise();

        Vector edgeVectorSum = edgeVector1 + edgeVector2 + edgeVector3 + edgeVector4;

        if (projection != Projection::cartesian)
        {
            edgeVectorSum.x() *= constants::conversion::radToDeg * constants::geometric::inverse_earth_radius;
        }

        double edgeVectorSumLength = edgeVectorSum.length();

        for (UInt i = 1; i <= constants::numeric::MaximumNumberOfCircumcentreIterations; ++i)
        {
            Vector delta1 = GetDelta(midPoint1, centreOfMass, projection);
            Vector delta2 = GetDelta(midPoint2, centreOfMass, projection);
            Vector delta3 = GetDelta(midPoint3, centreOfMass, projection);
            Vector delta4 = GetDelta(midPoint4, centreOfMass, projection);

            double ds = dot(delta1, edgeVector1) + dot(delta2, edgeVector2) + dot(delta3, edgeVector3) + dot(delta4, edgeVector4);

            if (projection != Projection::cartesian)
            {
                ds *= constants::conversion::radToDeg * constants::geometric::inverse_earth_radius;
            }

            centreOfMass.x -= updateStepSize * ds * edgeVectorSum.x() * xTransformation;
            centreOfMass.y -= updateStepSize * ds * edgeVectorSum.y();

            if (ds * edgeVectorSumLength < circumcentreTolerance)
            {
                break;
            }
        }
    }
    else
    {
        for (UInt j = 1; j <= constants::numeric::MaximumNumberOfCircumcentreIterations; ++j)
        {
            Vector edgeVectorSum(0.0, 0.0);
            double ds = 0.0;

            for (UInt i = 0; i < numberOfPointsOpenedPolygon; ++i)
            {
                const auto nextNode = NextCircularForwardIndex(i, numberOfPointsOpenedPolygon);

                Point midPoint = 0.5 * (polygon[i] + polygon[nextNode]);
                Vector edgeVector = static_cast<Vector>(NormalVector(polygon[i], polygon[nextNode], midPoint, projection));
                Vector delta = GetDelta(midPoint, centreOfMass, projection);

                edgeVector.normalise();
                ds += dot(delta, edgeVector);
                edgeVectorSum += edgeVector;
            }

            if (projection != Projection::cartesian)
            {
                ds *= constants::conversion::radToDeg * constants::geometric::inverse_earth_radius;
            }

            centreOfMass.x -= updateStepSize * ds * edgeVectorSum.x() * xTransformation;
            centreOfMass.y -= updateStepSize * ds * edgeVectorSum.y();

            if (j > 1 && ds * edgeVectorSum.length() < circumcentreTolerance)
            {
                break;
            }
        }
    }

    area = ComputeArea(polygon, projection);
    TraversalDirection direction = area > 0.0 ? TraversalDirection::AntiClockwise : TraversalDirection::Clockwise;

    area = std::abs(area) < minArea ? minArea : area;

    return {std::abs(area), centreOfMass, direction};
}

std::tuple<double, meshkernel::Point, meshkernel::TraversalDirection> meshkernel::Polygon::FaceAreaAndCenterOfMass(const std::vector<Point>& nodes,
                                                                                                                   const std::vector<UInt>& nodeIndices,
                                                                                                                   const Projection projection,
                                                                                                                   const bool isClosed)
{

    if (nodeIndices.size() < constants::geometric::numNodesInTriangle)
    {
        throw std::invalid_argument("FaceAreaAndCenterOfMass: The polygon has less than 3 unique nodes.");
    }

    const auto numberOfPointsOpenedPolygon = static_cast<UInt>(nodeIndices.size()) - (isClosed ? 1 : 0);

    std::vector<Point> polygonPoints(numberOfPointsOpenedPolygon + 1);

    for (size_t i = 0; i < numberOfPointsOpenedPolygon; ++i)
    {
        polygonPoints[i] = nodes[nodeIndices[i]];
    }

    polygonPoints[numberOfPointsOpenedPolygon] = nodes[nodeIndices[0]];

    return FaceAreaAndCenterOfMass(polygonPoints, projection);
}

std::tuple<double, meshkernel::Point, meshkernel::TraversalDirection> meshkernel::Polygon::FaceAreaAndCenterOfMass() const
{
    return FaceAreaAndCenterOfMass(m_nodes, m_projection);
}

double meshkernel::Polygon::PerimeterLength() const
{
    if (m_nodes.size() <= 1)
    {
        return 0.0;
    }

    double perimeter = 0.0;

    for (size_t i = 1; i < m_nodes.size(); ++i)
    {
        perimeter += ComputeDistance(m_nodes[i - 1], m_nodes[i], m_projection);
    }

    return perimeter;
}

std::vector<meshkernel::Point> meshkernel::Polygon::ComputeOffset(double distance, const bool innerAndOuter) const
{
    std::vector<Vector> normalVectors(m_nodes.size());

    Vector normal;
    Vector firstNormal;
    // Initialise the previous normal with the normal from the last segment of the polygon
    // The last segment consists of the second to last point and the first point.
    // Second to last because the polygon is closed, so the last point is equal to the first.
    Vector previousNormal = ComputeNormalToline(m_nodes[m_nodes.size() - 2], m_nodes[0], m_projection);

    for (size_t i = 0; i < m_nodes.size(); ++i)
    {
        if (i < m_nodes.size() - 1)
        {
            normal = ComputeNormalToline(m_nodes[i], m_nodes[i + 1], m_projection);
        }
        else
        {
            normal = firstNormal;
        }

        if (i == 0)
        {
            firstNormal = normal;
        }

        const double factor = 1.0 / (1.0 + dot(previousNormal, normal));
        normalVectors[i] = factor * (previousNormal + normal);

        previousNormal = normal;
    }

    std::vector<Point> offsetPoints(m_nodes.size() + (innerAndOuter ? m_nodes.size() + 1 : 0), Point());

    // negative sign introduced because normal vector pointing inward
    distance = -distance;

    // TODO should this be for spherical accurate too.
    // The perhaps should either:
    // 1. a function to determine if m_projection is a spherical kind (either spherical or spherical-accurate, or any other spherical type)
    // 2. IFF projection will only include the current 3 types (check it is not cartesian) (should really add an uninitialised value)
    if (m_projection == Projection::spherical)
    {
        distance = distance / (constants::geometric::earth_radius * constants::conversion::degToRad);
    }

    for (UInt i = 0; i < m_nodes.size(); ++i)
    {
        Vector normal = distance * normalVectors[i];

        // TODO should this be for spherical accurate too.
        if (m_projection == Projection::spherical)
        {
            // TODO is it worth having a special function for this, it must be elsewhere in the code too.
            normal.x() /= std::cos((m_nodes[i].y + 0.5 * normal.y()) * constants::conversion::degToRad);
        }

        offsetPoints[i] = m_nodes[i] + normal;

        if (innerAndOuter)
        {
            offsetPoints[i + m_nodes.size() + 1] = m_nodes[i] - normal;
        }
    }

    return offsetPoints;
}

std::tuple<double, double> meshkernel::Polygon::SegmentLengthExtrema() const
{

    if (m_nodes.size() <= 1)
    {
        return {constants::missing::doubleValue, constants::missing::doubleValue};
    }

    double minimumSegmentLength = std::numeric_limits<double>::max();
    double maximumSegmentLength = 0.0;

    for (size_t i = 1; i < m_nodes.size(); ++i)
    {
        double segmentLength = ComputeDistance(m_nodes[i - 1], m_nodes[i], m_projection);

        minimumSegmentLength = std::min(minimumSegmentLength, segmentLength);
        maximumSegmentLength = std::max(maximumSegmentLength, segmentLength);
    }

    return {minimumSegmentLength, maximumSegmentLength};
}
