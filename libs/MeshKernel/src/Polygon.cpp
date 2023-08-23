#include "MeshKernel/Polygon.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/LandBoundary.hpp"
#include "MeshKernel/Mesh.hpp"
#include "MeshKernel/Operations.hpp"

meshkernel::Polygon::Polygon(const std::vector<Point>& points,
                             Projection projection) : m_points(points), m_projection(projection)
{
    // TODO check number of points >= 3 (probably 4 to include the closure)
    // Polygon should not contain invalid points

    if (m_projection == Projection::spherical)
    {
        // TODO SHould this be called for spherical accurate too?
        TranslateSphericalCoordinates(m_points);
    }

    m_boundingBox.Reset(m_points);
}

meshkernel::Polygon::Polygon(std::vector<Point>&& points,
                             Projection projection) : m_points(points), m_projection(projection)
{
    // TODO check number of points >= 3 (probably 4 to include the closure)
    // Polygon should not contain invalid points

    if (m_projection == Projection::spherical)
    {
        // TODO SHould this be called for spherical accurate too?
        TranslateSphericalCoordinates(m_points);
    }

    m_boundingBox.Reset(m_points);
}

meshkernel::Polygon& meshkernel::Polygon::operator=(const Polygon& copy)
{

    if (this != &copy)
    {
        m_points = copy.m_points;
        m_projection = copy.m_projection;
        m_boundingBox = copy.m_boundingBox;
    }

    return *this;
}

meshkernel::Polygon& meshkernel::Polygon::operator=(Polygon&& copy)
{

    if (this != &copy)
    {
        m_points = std::move(copy.m_points);
        // TODO add undefined enum for Projection
        m_projection = copy.m_projection;
        // TODO should the BB be invalidated.
        m_boundingBox = copy.m_boundingBox;
    }

    return *this;
}

void meshkernel::Polygon::Reset(const std::vector<Point>& points,
                                Projection projection)
{
    m_projection = projection;
    m_points = points;

    if (m_projection == Projection::spherical)
    {
        // TODO SHould this be called for spherical accurate too?
        TranslateSphericalCoordinates(m_points);
    }

    m_boundingBox.Reset(m_points);
}

bool meshkernel::Polygon::ContainsCartesian(const Point& point) const
{

    int windingNumber = 0;

    for (size_t n = 0; n < m_points.size() - 1; n++)
    {
        // TODO always Cartesian
        // So Dx and Dy can be simplified (no branching)
        // Then for 2 or more points, return multiple cross product values
        const auto crossProductValue = crossProduct(m_points[n], m_points[n + 1], m_points[n], point, Projection::cartesian);

        if (IsEqual(crossProductValue, 0.0))
        {
            // point on the line
            return true;
        }

        if (m_points[n].y <= point.y) // an upward crossing
        {
            if (m_points[n + 1].y > point.y && crossProductValue > 0.0)

            {
                ++windingNumber; // have  a valid up intersect
            }
        }
        else
        {
            if (m_points[n + 1].y <= point.y && crossProductValue < 0.0) // a downward crossing
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
    bool isInPolygon = false;

    // get 3D polygon coordinates
    std::vector<Cartesian3DPoint> cartesian3DPoints;
    cartesian3DPoints.reserve(Size());

    for (UInt i = 0; i < m_points.size(); ++i)
    {
        cartesian3DPoints.emplace_back(SphericalToCartesian3D(m_points[i]));
    }

    // enlarge around polygon
    const double enlargementFactor = 1.000001;

    // TODO set to centre?
    Point polygonCenter;
    const Cartesian3DPoint polygonCenterCartesian3D{SphericalToCartesian3D(polygonCenter)};
    for (UInt i = 0; i < m_points.size(); ++i)
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
    for (UInt i = 0; i < m_points.size() - 1; ++i)
    {
        const auto nextNode = NextCircularForwardIndex(i, static_cast<UInt>(m_points.size()));
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

    if (inside == 1)
    {
        isInPolygon = true;
    }

    return inside == 1;
    return isInPolygon;
}

bool meshkernel::Polygon::Contains(const Point& pnt) const
{

    if (!pnt.IsValid())
    {
        throw ConstraintError("Point is not valid");
    }

    if (m_points.empty())
    {
        return true;
    }

    if (m_points.size() < Mesh::m_numNodesInTriangle)
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
    else
    {
        // projection = Projection::sphericalAccurate
        return ContainsSphericalAccurate(pnt);
    }
}

// TODO does this need start and end points, probably not all the polygon is to be snapped
void meshkernel::Polygon::SnapToLandBoundary(const size_t startIndex, const size_t endIndex, const LandBoundary& landBoundary)
{

    if (startIndex > endIndex || endIndex >= m_points.size())
    {
        throw ConstraintError(VariadicErrorMessage("Polygon::SnapToLandBoundary: The indices are not valid: {}, {}.", startIndex, endIndex));
    }

    // snap polygon section to land boundary
    for (size_t i = startIndex; i <= endIndex; ++i)
    {
        if (m_points[i].IsValid())
        {
            m_points[i] = landBoundary.FindNearestPoint(m_points[i], m_projection);
        }
    }

    if (m_projection == Projection::spherical)
    {
        // TODO Should this be called for spherical accurate too?
        // TODO WHere did I find this? Is it correct?
        TranslateSphericalCoordinates(m_points);
    }

    // Now update the bounding box
    m_boundingBox.Reset(m_points);
}

std::vector<double> meshkernel::Polygon::EdgeLengths() const
{
    std::vector<double> edgeLengths;
    edgeLengths.reserve(m_points.size());

    for (size_t p = 0; p < m_points.size(); ++p)
    {
        size_t firstNode = p;
        size_t secondNode = p + 1;

        if (secondNode == m_points.size())
        {
            secondNode = 0;
        }

        edgeLengths.emplace_back(ComputeDistance(m_points[firstNode], m_points[secondNode], m_projection));
    }

    return edgeLengths;
}

std::vector<meshkernel::Point> meshkernel::Polygon::Refine(const size_t startIndex, const size_t endIndex, const double refinementDistance) const
{

    // UInt polygonIndex;

    if (startIndex > endIndex || endIndex >= m_points.size())
    {
        throw ConstraintError(VariadicErrorMessage("The indices are not valid: {}, {}.", startIndex, endIndex));
    }

    //--------------------------------

    // const auto& [outerStart, outerEnd] = m_outer_polygons_indices[polygonIndex];

    // const auto edgeLengths = EdgeLengths();
    // std::vector<double> nodeLengthCoordinate(edgeLengths.size());
    // nodeLengthCoordinate[0] = 0.0;

    // for (UInt i = 1; i < edgeLengths.size(); ++i)
    // {
    //     nodeLengthCoordinate[i] = nodeLengthCoordinate[i - 1] + edgeLengths[i - 1];
    // }

    // Approximate number of nodes in the refined sections.
    // const UInt numNodesRefinedPart = static_cast<UInt>(std::ceil((nodeLengthCoordinate[endIndex] - nodeLengthCoordinate[startIndex]) / refinementDistance)) + endIndex - startIndex;
    // UInt numNodesNotRefinedPart = startIndex - outerStart + outerEnd - endIndex;
    // Approximate the number of nodes in the refined polygon.

    std::vector<Point> refinedPolygon;
    refinedPolygon.reserve(m_points.size());

    // Add nodes before the section to be refined
    for (size_t i = 0; i <= startIndex; ++i)
    {
        refinedPolygon.emplace_back(m_points[i]);
    }

    // Refine each line segment.
    for (size_t i = startIndex; i < endIndex; ++i)
    {
        // Line segment starting point.
        Point p = m_points[i];
        const double segmentLength = ComputeDistance(m_points[i], m_points[i + 1], m_projection);
        // Refined segment step size.
        const Point delta = (m_points[i + 1] - m_points[i]) * refinementDistance / segmentLength;
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
        refinedPolygon.emplace_back(m_points[i + 1]);
    }

    // Add nodes after the section to be refined
    for (size_t i = endIndex + 1; i < m_points.size(); ++i)
    {
        refinedPolygon.emplace_back(m_points[i]);
    }

    return refinedPolygon;
}

std::tuple<double, meshkernel::Point, bool> meshkernel::Polygon::FaceAreaAndCenterOfMass() const
{

    // TODO why size - 1? If open polygon?
    if (m_points.size() - 1 < Mesh::m_numNodesInTriangle)
    {
        throw std::invalid_argument("FaceAreaAndCenterOfMass: The polygon has less than 3 unique nodes.");
    }

    double area = 0.0;
    double xCenterOfMass = 0.0;
    double yCenterOfMass = 0.0;
    const double minArea = 1e-8;
    const Point reference = ReferencePoint(m_points, m_projection);
    const auto numberOfPointsOpenedPolygon = static_cast<UInt>(m_points.size()) - 1;

    for (UInt n = 0; n < numberOfPointsOpenedPolygon; n++)
    {
        const auto nextNode = NextCircularForwardIndex(n, numberOfPointsOpenedPolygon);
        double dx0 = GetDx(reference, m_points[n], m_projection);
        double dy0 = GetDy(reference, m_points[n], m_projection);
        const double dx1 = GetDx(reference, m_points[nextNode], m_projection);
        const double dy1 = GetDy(reference, m_points[nextNode], m_projection);

        const double xc = 0.5 * (dx0 + dx1);
        const double yc = 0.5 * (dy0 + dy1);

        dx0 = GetDx(m_points[n], m_points[nextNode], m_projection);
        dy0 = GetDy(m_points[n], m_points[nextNode], m_projection);

        // Rotate by pi/2
        const double dsx = dy0;
        const double dsy = -dx0;
        const double xds = xc * dsx + yc * dsy;
        area = area + 0.5 * xds;

        xCenterOfMass = xCenterOfMass + xds * xc;
        yCenterOfMass = yCenterOfMass + xds * yc;
    }

    bool isCounterClockWise = area > 0.0;

    area = std::abs(area) < minArea ? minArea : area;

    const double fac = 1.0 / (3.0 * area);
    xCenterOfMass = fac * xCenterOfMass;
    yCenterOfMass = fac * yCenterOfMass;

    // TODO SHould this also apply to spheciral accurate?
    if (m_projection == Projection::spherical)
    {
        yCenterOfMass = yCenterOfMass / (constants::geometric::earth_radius * constants::conversion::degToRad);
        xCenterOfMass = xCenterOfMass / (constants::geometric::earth_radius * constants::conversion::degToRad * std::cos((yCenterOfMass + reference.y) * constants::conversion::degToRad));
    }

    Point centreOfMass(xCenterOfMass, yCenterOfMass);
    centreOfMass += reference;

    return {area, centreOfMass, isCounterClockWise};
}

double meshkernel::Polygon::ClosedPerimeterLength() const
{
    if (m_points.size() <= 1)
    {
        return 0.0;
    }

    // TODO m_points.front, m_points.back
    // TOOD are all polygons closed?
    double perimeter = IsClosed() ? 0.0 : ComputeDistance(m_points[0], m_points[m_points.size() - 1], m_projection);

    for (size_t i = 1; i < m_points.size(); ++i)
    {
        perimeter += ComputeDistance(m_points[i - 1], m_points[i], m_projection);
    }

    return perimeter;
}

// meshkernel::Polygon meshkernel::Polygon::Displace(double displacement) const
// {
//     PointArray
// }
