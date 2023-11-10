#include <cmath>

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

std::vector<meshkernel::Point> meshkernel::Polygon::Refine(const size_t startIndex, const size_t endIndex, const double refinementDistance) const
{

    if (startIndex > endIndex || endIndex >= m_nodes.size())
    {
        throw ConstraintError("The indices are not valid: {}, {}.", startIndex, endIndex);
    }

    std::vector<Point> refinedPolygon;
    refinedPolygon.reserve(m_nodes.size());

    // Add nodes before the section to be refined
    for (size_t i = 0; i <= startIndex; ++i)
    {
        refinedPolygon.emplace_back(m_nodes[i]);
    }

    // Refine each line segment.
    for (size_t i = startIndex; i < endIndex; ++i)
    {
        // Line segment starting point.
        Point p = m_nodes[i];
        const double segmentLength = ComputeDistance(m_nodes[i], m_nodes[i + 1], m_projection);
        // Refined segment step size.
        const Point delta = (m_nodes[i + 1] - m_nodes[i]) * refinementDistance / segmentLength;
        double lengthAlongInterval = refinementDistance;

        // Exit when the lengthAlongInterval is greater or equal than segmentLength
        // To prevent very small refined segment lengths, also exit when lengthAlongInterval is a small fraction less (defined by refinementTolerance)
        while (lengthAlongInterval < segmentLength && !IsEqual(lengthAlongInterval, segmentLength, constants::geometric::refinementTolerance))
        {
            p += delta;
            lengthAlongInterval += refinementDistance;
            refinedPolygon.emplace_back(p);
        }

        // Add last node to the refined polygon point sequence.
        refinedPolygon.emplace_back(m_nodes[i + 1]);
    }

    // Add nodes after the section to be refined
    for (size_t i = endIndex + 1; i < m_nodes.size(); ++i)
    {
        refinedPolygon.emplace_back(m_nodes[i]);
    }

    return refinedPolygon;
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
