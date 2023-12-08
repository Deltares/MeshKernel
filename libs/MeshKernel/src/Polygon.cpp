#include <cmath>

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
            // point on the line
            return true;
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

    if (startIndex > endIndex || endIndex >= m_nodes.size())
    {
        throw ConstraintError("The indices are not valid: {}, {}.", startIndex, endIndex);
    }

    // snap polygon section to land boundary
    for (size_t i = startIndex; i <= endIndex; ++i)
    {
        if (m_nodes[i].IsValid())
        {
            m_nodes[i] = landBoundary.FindNearestPoint(m_nodes[i], m_projection);
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

namespace
{
    /// @brief Refines the segment between two polygon nodes, starting with the node specified by the iterator up to,
    /// but not including, the next node.
    /// @param refinedPolygon     [in,out] a buffer of points into which the refined points are written
    /// @param nodeIterator       [in] position in the original, unrefined polygon that contains the first point of
    ///                           the refinement
    /// @param refinementDistance [in] the distance between two refined nodes
    /// @param projection         [in] the projection used for computing the length of the segment to be refined
    void RefineSegment(std::vector<meshkernel::Point>& refinedPolygon,
                       const std::vector<meshkernel::Point>::const_iterator& nodeIterator,
                       const double refinementDistance,
                       const meshkernel::Projection projection)
    {
        // Line segment starting point.
        const auto& n0 = *nodeIterator;
        const auto& n1 = *std::next(nodeIterator);

        refinedPolygon.push_back(n0);

        const double segmentLength = ComputeDistance(n0, n1, projection);

        int n = std::lround(segmentLength / refinementDistance);
        const double refinedLength = n * refinementDistance;
        if (refinedLength > segmentLength ||
            meshkernel::IsEqual(refinedLength, segmentLength, meshkernel::constants::geometric::refinementTolerance))
        {
            --n;
        }

        // Refined segment step size.
        const meshkernel::Point delta = (n1 - n0) * refinementDistance / segmentLength;
        for (auto i = 1; i <= n; ++i)
        {
            refinedPolygon.push_back(n0 + i * delta);
        }
    }

    void computeAverageLengths(const std::vector<double>& cumulativeDistances, std::vector<double>& averageDistances)
    {
        averageDistances[0] = cumulativeDistances[1] - cumulativeDistances.front();
        averageDistances.back() = cumulativeDistances.back() - cumulativeDistances[cumulativeDistances.size() - 2];
        for (meshkernel::UInt i = 1; i < averageDistances.size() - 1; ++i)
        {
            averageDistances[i] = 0.5 * (cumulativeDistances[i + 1] - cumulativeDistances[i - 1]);
        }
    }

    void smoothCumulativeDistance(const std::vector<double>& averageDistances, std::vector<double>& cumulativeDistances)
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

    void smoothAverageLengths(const std::vector<double>& cumulativeDistances,
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

    meshkernel::Point interpolatePointOnPolyline(const std::vector<meshkernel::Point>& points,
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
} // namespace

std::vector<meshkernel::Point> meshkernel::Polygon::Refine(const size_t startIndex, const size_t endIndex, const double refinementDistance) const
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

    auto computeDistance = [&projection = std::as_const(m_projection)](auto l, auto p)
    {
        return l + ComputeDistance(p, *std::next(&p), projection);
    };

    if (startIndex < endIndex)
    {
        // estimate the size of the resulting polygon
        const double length = std::accumulate(from, to, 0.0, computeDistance);
        const size_t sizeEstimate = m_nodes.size() + std::lround(std::ceil(length / refinementDistance));
        refinedPolygon.reserve(sizeEstimate);

        // copy the first segments, up to startIndex
        std::copy(m_nodes.begin(), from, std::back_inserter(refinedPolygon));

        // refine edges starting with startIndex
        for (auto it = from; it != to; ++it)
        {
            RefineSegment(refinedPolygon, it, refinementDistance, m_projection);
        }

        // copy the nodes starting with endIndex until the end
        std::copy(to, m_nodes.end(), std::back_inserter(refinedPolygon));
    }
    else
    {
        // estimate the size of the resulting polygon
        const auto last = std::prev(m_nodes.end());
        const double length =
            std::accumulate(m_nodes.begin(), from, 0.0, computeDistance) +
            std::accumulate(to, last, 0.0, computeDistance);

        const size_t sizeEstimate = m_nodes.size() + static_cast<size_t>(length / refinementDistance);
        refinedPolygon.reserve(sizeEstimate);

        // refine edges from the start to endIndex
        for (auto it = m_nodes.begin(); it != to; ++it)
            RefineSegment(refinedPolygon, it, refinementDistance, m_projection);

        // copy from endIndex up to startIndex
        std::copy(to, from, std::back_inserter(refinedPolygon));

        // refine edges from startIndex to the end of the polygon
        for (auto it = from; it != last; ++it)
            RefineSegment(refinedPolygon, it, refinementDistance, m_projection);
        refinedPolygon.push_back(m_nodes.back());
    }
    return refinedPolygon;
}

std::vector<meshkernel::Point> meshkernel::Polygon::LinearRefine(const size_t startIndex, const size_t endIndex) const
{
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
    for (UInt i = 0; i < polygonNodes.size(); ++i)
    {
        auto polygonNodeIndex = i + startIndex;
        if (polygonNodeIndex >= m_nodes.size())
        {
            polygonNodeIndex = polygonNodeIndex - m_nodes.size();
        }
        polygonNodes[i] = m_nodes[polygonNodeIndex];
    }

    std::vector<double> cumulativeDistances(numPolygonNodes, 0.0);
    cumulativeDistances[0] = 0.0;
    for (UInt i = 1; i < polygonNodes.size(); ++i)
    {
        cumulativeDistances[i] = cumulativeDistances[i - 1] + ComputeDistance(polygonNodes[i], polygonNodes[i - 1], m_projection);
    }
    std::vector<double> initialCumulativeDistances(cumulativeDistances);

    computeAverageLengths(cumulativeDistances, averageLengths);
    const double firstLength = averageLengths.front();
    const double lastLength = averageLengths.back();
    double cumulativeDistanceTarget = 0.0;
    const UInt maxInnerIter = 20;
    const UInt maxOuterIter = numPolygonNodes;
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

        double minRatio = 1.0e9;
        UInt minRatioIndex = constants::missing::uintValue;
        double maxRatio = -1.0e9;
        UInt maxRatioIndex = constants::missing::uintValue;
        double minLength = 1.0e30;

        for (UInt i = 0; i < averageLengths.size() - 1; ++i)
        {
            minLength = std::min(averageLengths[i], minLength);
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

    Point centreOfMass(0.0, 0.0);
    double area = 0.0;

    const double minArea = 1e-8;
    const Point reference = ReferencePoint(polygon, projection);
    const auto numberOfPointsOpenedPolygon = static_cast<UInt>(polygon.size()) - 1;

    if (numberOfPointsOpenedPolygon == constants::geometric::numNodesInTriangle)
    {
        Vector delta1 = GetDelta(reference, polygon[0], projection);
        Vector delta2 = GetDelta(reference, polygon[1], projection);
        Vector delta3 = GetDelta(reference, polygon[2], projection);

        Vector middle1 = 0.5 * (delta1 + delta2);
        Vector middle2 = 0.5 * (delta2 + delta3);
        Vector middle3 = 0.5 * (delta3 + delta1);

        delta1 = GetDelta(polygon[0], polygon[1], projection);
        delta2 = GetDelta(polygon[1], polygon[2], projection);
        delta3 = GetDelta(polygon[2], polygon[0], projection);

        double xds1 = delta1.y() * middle1.x() - delta1.x() * middle1.y();
        double xds2 = delta2.y() * middle2.x() - delta2.x() * middle2.y();
        double xds3 = delta3.y() * middle3.x() - delta3.x() * middle3.y();

        area = 0.5 * (xds1 + xds2 + xds3);

        centreOfMass += xds1 * middle1;
        centreOfMass += xds2 * middle2;
        centreOfMass += xds3 * middle3;
    }
    else if (numberOfPointsOpenedPolygon == constants::geometric::numNodesInQuadrilateral)
    {
        Vector delta1 = GetDelta(reference, polygon[0], projection);
        Vector delta2 = GetDelta(reference, polygon[1], projection);
        Vector delta3 = GetDelta(reference, polygon[2], projection);
        Vector delta4 = GetDelta(reference, polygon[3], projection);

        Vector middle1 = 0.5 * (delta1 + delta2);
        Vector middle2 = 0.5 * (delta2 + delta3);
        Vector middle3 = 0.5 * (delta3 + delta4);
        Vector middle4 = 0.5 * (delta4 + delta1);

        delta1 = GetDelta(polygon[0], polygon[1], projection);
        delta2 = GetDelta(polygon[1], polygon[2], projection);
        delta3 = GetDelta(polygon[2], polygon[3], projection);
        delta4 = GetDelta(polygon[3], polygon[0], projection);

        double xds1 = delta1.y() * middle1.x() - delta1.x() * middle1.y();
        double xds2 = delta2.y() * middle2.x() - delta2.x() * middle2.y();
        double xds3 = delta3.y() * middle3.x() - delta3.x() * middle3.y();
        double xds4 = delta4.y() * middle4.x() - delta4.x() * middle4.y();

        area = 0.5 * (xds1 + xds2 + xds3 + xds4);

        centreOfMass += xds1 * middle1;
        centreOfMass += xds2 * middle2;
        centreOfMass += xds3 * middle3;
        centreOfMass += xds4 * middle4;
    }
    else
    {

        for (UInt n = 0; n < numberOfPointsOpenedPolygon; ++n)
        {
            const auto nextNode = NextCircularForwardIndex(n, numberOfPointsOpenedPolygon);

            Vector delta = GetDelta(reference, polygon[n], projection);
            Vector deltaNext = GetDelta(reference, polygon[nextNode], projection);
            Vector middle = 0.5 * (delta + deltaNext);
            delta = GetDelta(polygon[n], polygon[nextNode], projection);

            // Rotate by 3pi/2
            Vector normal(delta.y(), -delta.x());
            double xds = dot(normal, middle);
            area += 0.5 * xds;

            centreOfMass += xds * middle;
        }
    }

    TraversalDirection direction = area > 0.0 ? TraversalDirection::AntiClockwise : TraversalDirection::Clockwise;

    area = std::abs(area) < minArea ? minArea : area;
    centreOfMass *= 1.0 / (3.0 * area);

    // TODO SHould this also apply to spheciral accurate?
    if (projection == Projection::spherical)
    {
        centreOfMass.y /= (constants::geometric::earth_radius * constants::conversion::degToRad);
        centreOfMass.x /= (constants::geometric::earth_radius * constants::conversion::degToRad * std::cos((centreOfMass.y + reference.y) * constants::conversion::degToRad));
    }

    centreOfMass += reference;

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

    Point centreOfMass(0.0, 0.0);
    double area = 0.0;

    const double minArea = 1e-8;
    const Point reference = ReferencePoint(nodes, nodeIndices, projection);
    const auto numberOfPointsOpenedPolygon = static_cast<UInt>(nodeIndices.size()) - (isClosed ? 1 : 0);

    if (numberOfPointsOpenedPolygon == constants::geometric::numNodesInTriangle)
    {
        Vector delta1 = GetDelta(reference, nodes[nodeIndices[0]], projection);
        Vector delta2 = GetDelta(reference, nodes[nodeIndices[1]], projection);
        Vector delta3 = GetDelta(reference, nodes[nodeIndices[2]], projection);

        Vector middle1 = 0.5 * (delta1 + delta2);
        Vector middle2 = 0.5 * (delta2 + delta3);
        Vector middle3 = 0.5 * (delta3 + delta1);

        delta1 = GetDelta(nodes[nodeIndices[0]], nodes[nodeIndices[1]], projection);
        delta2 = GetDelta(nodes[nodeIndices[1]], nodes[nodeIndices[2]], projection);
        delta3 = GetDelta(nodes[nodeIndices[2]], nodes[nodeIndices[0]], projection);

        double xds1 = delta1.y() * middle1.x() - delta1.x() * middle1.y();
        double xds2 = delta2.y() * middle2.x() - delta2.x() * middle2.y();
        double xds3 = delta3.y() * middle3.x() - delta3.x() * middle3.y();

        area = 0.5 * (xds1 + xds2 + xds3);

        centreOfMass.x = xds1 * middle1.x();
        centreOfMass.y = xds1 * middle1.y();
        centreOfMass += xds2 * middle2;
        centreOfMass += xds3 * middle3;
    }
    else if (numberOfPointsOpenedPolygon == constants::geometric::numNodesInQuadrilateral)
    {
        Vector delta1 = GetDelta(reference, nodes[nodeIndices[0]], projection);
        Vector delta2 = GetDelta(reference, nodes[nodeIndices[1]], projection);
        Vector delta3 = GetDelta(reference, nodes[nodeIndices[2]], projection);
        Vector delta4 = GetDelta(reference, nodes[nodeIndices[3]], projection);

        Vector middle1 = 0.5 * (delta1 + delta2);
        Vector middle2 = 0.5 * (delta2 + delta3);
        Vector middle3 = 0.5 * (delta3 + delta4);
        Vector middle4 = 0.5 * (delta4 + delta1);

        delta1 = GetDelta(nodes[nodeIndices[0]], nodes[nodeIndices[1]], projection);
        delta2 = GetDelta(nodes[nodeIndices[1]], nodes[nodeIndices[2]], projection);
        delta3 = GetDelta(nodes[nodeIndices[2]], nodes[nodeIndices[3]], projection);
        delta4 = GetDelta(nodes[nodeIndices[3]], nodes[nodeIndices[0]], projection);

        double xds1 = delta1.y() * middle1.x() - delta1.x() * middle1.y();
        double xds2 = delta2.y() * middle2.x() - delta2.x() * middle2.y();
        double xds3 = delta3.y() * middle3.x() - delta3.x() * middle3.y();
        double xds4 = delta4.y() * middle4.x() - delta4.x() * middle4.y();

        area = 0.5 * (xds1 + xds2 + xds3 + xds4);

        centreOfMass.x = xds1 * middle1.x();
        centreOfMass.y = xds1 * middle1.y();
        centreOfMass += xds2 * middle2;
        centreOfMass += xds3 * middle3;
        centreOfMass += xds4 * middle4;
    }
    else
    {

        for (UInt n = 0; n < numberOfPointsOpenedPolygon; ++n)
        {
            const auto nextNode = NextCircularForwardIndex(n, numberOfPointsOpenedPolygon);

            Vector delta = GetDelta(reference, nodes[nodeIndices[n]], projection);
            Vector deltaNext = GetDelta(reference, nodes[nodeIndices[nextNode]], projection);
            Vector middle = 0.5 * (delta + deltaNext);
            delta = GetDelta(nodes[nodeIndices[n]], nodes[nodeIndices[nextNode]], projection);

            // Rotate by 3pi/2
            Vector normal(delta.y(), -delta.x());
            double xds = dot(normal, middle);
            area += 0.5 * xds;

            centreOfMass += xds * middle;
        }
    }

    TraversalDirection direction = area > 0.0 ? TraversalDirection::AntiClockwise : TraversalDirection::Clockwise;

    area = std::abs(area) < minArea ? minArea : area;
    centreOfMass *= 1.0 / (3.0 * area);

    // TODO SHould this also apply to spheciral accurate?
    if (projection == Projection::spherical)
    {
        centreOfMass.y /= (constants::geometric::earth_radius * constants::conversion::degToRad);
        centreOfMass.x /= (constants::geometric::earth_radius * constants::conversion::degToRad * std::cos((centreOfMass.y + reference.y) * constants::conversion::degToRad));
    }

    centreOfMass += reference;

    return {std::abs(area), centreOfMass, direction};
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
