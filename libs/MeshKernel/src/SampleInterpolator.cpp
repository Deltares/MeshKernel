//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2024.
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

#include "MeshKernel/SampleInterpolator.hpp"
#include "MeshKernel/Operations.hpp"

void meshkernel::SampleTriangulationInterpolator::SetDataSpan(const int propertyId, const std::span<const double>& sampleData)
{
    if (m_triangulation.NumberOfNodes() != sampleData.size())
    {
        throw ConstraintError("The sample data array does not have the same number of elements as the number of nodes in the triangulation: {} /= {}",
                              m_triangulation.NumberOfNodes(), sampleData.size());
    }

    m_sampleData[propertyId].assign(sampleData.begin(), sampleData.end());
}

void meshkernel::SampleTriangulationInterpolator::InterpolateSpan(const int propertyId, const std::span<const Point>& iterpolationNodes, std::span<double>& result) const
{
    if (!Contains(propertyId))
    {
        throw ConstraintError("Sample interpolator does not contain the id: {}.", propertyId);
    }

    if (iterpolationNodes.size() != result.size())
    {
        throw ConstraintError("The arrays for interpolation nodes and the results are different sizes: {} /= {}",
                              iterpolationNodes.size(), result.size());
    }

    const std::vector<double>& propertyValues = m_sampleData.at(propertyId);

    for (size_t i = 0; i < iterpolationNodes.size(); ++i)
    {
        result[i] = constants::missing::doubleValue;

        if (!iterpolationNodes[i].IsValid())
        {
            continue;
        }

        const UInt elementId = m_triangulation.FindNearestFace(iterpolationNodes[i]);

        if (elementId == constants::missing::uintValue)
        {
            continue;
        }

        if (m_triangulation.PointIsInElement(iterpolationNodes[i], elementId))
        {
            result[i] = InterpolateOnElement(elementId, iterpolationNodes[i], propertyValues);
        }
    }
}

double meshkernel::SampleTriangulationInterpolator::InterpolateOnElement(const UInt elementId, const Point& interpolationPoint, const std::vector<double>& sampleValues) const
{
    double result = constants::missing::doubleValue;

    auto [id1, id2, id3] = m_triangulation.GetNodeIds(elementId);

    if (sampleValues[id1] == constants::missing::doubleValue ||
        sampleValues[id2] == constants::missing::doubleValue ||
        sampleValues[id3] == constants::missing::doubleValue) [[unlikely]]
    {
        return result;
    }

    auto [p1, p2, p3] = m_triangulation.GetNodes(elementId);

    const double a11 = GetDx(p1, p2, m_triangulation.GetProjection());
    const double a21 = GetDy(p1, p2, m_triangulation.GetProjection());

    const double a12 = GetDx(p1, p3, m_triangulation.GetProjection());
    const double a22 = GetDy(p1, p3, m_triangulation.GetProjection());

    const double b1 = GetDx(p1, interpolationPoint, m_triangulation.GetProjection());
    const double b2 = GetDy(p1, interpolationPoint, m_triangulation.GetProjection());

    const double det = a11 * a22 - a12 * a21;

    if (std::abs(det) < 1e-12) [[unlikely]]
    {
        return result;
    }

    const double rlam = (a22 * b1 - a12 * b2) / det;
    const double rmhu = (a11 * b2 - a21 * b1) / det;

    result = sampleValues[id1] + rlam * (sampleValues[id2] - sampleValues[id1]) + rmhu * (sampleValues[id3] - sampleValues[id1]);

    return result;
}

double meshkernel::SampleTriangulationInterpolator::InterpolateValue(const int propertyId, const Point& evaluationPoint) const
{

    if (!Contains(propertyId))
    {
        throw ConstraintError("Sample interpolator does not contain the id: {}.", propertyId);
    }

    double result = constants::missing::doubleValue;

    if (!evaluationPoint.IsValid())
    {
        return result;
    }

    const UInt elementId = m_triangulation.FindNearestFace(evaluationPoint);

    if (elementId == constants::missing::uintValue)
    {
        return result;
    }

    if (m_triangulation.PointIsInElement(evaluationPoint, elementId))
    {
        const std::vector<double>& propertyValues = m_sampleData.at(propertyId);
        result = InterpolateOnElement(elementId, evaluationPoint, propertyValues);
    }

    return result;
}

void meshkernel::SampleAveragingInterpolator::SetDataSpan(const int propertyId, const std::span<const double>& sampleData)
{
    if (m_samplePoints.size() != sampleData.size())
    {
        throw ConstraintError("The sample data array does not have the same number of elements as the number of nodes in the triangulation: {} /= {}",
                              m_samplePoints.size(), sampleData.size());
    }

    m_sampleData[propertyId].assign(sampleData.begin(), sampleData.end());
}

void meshkernel::SampleAveragingInterpolator::InterpolateSpan(const int propertyId, const std::span<const Point>& iterpolationNodes, std::span<double>& result) const
{
    const std::vector<double>& propertyValues = m_sampleData.at(propertyId);

    for (size_t i = 0; i < iterpolationNodes.size(); ++i)
    {
        result[i] = propertyValues[i];
    }
}

double meshkernel::SampleAveragingInterpolator::GetSearchRadiusSquared(const std::vector<Point>& searchPolygon,
                                                                       const Point& interpolationPoint,
                                                                       const Projection projection) const
{
    double result = std::numeric_limits<double>::lowest();

    for (const auto& value : searchPolygon)
    {
        auto const squaredDistance = ComputeSquaredDistance(interpolationPoint, value, projection);
        std::cout << "search radius: " << value.x << ", " << value.y << "  " << squaredDistance << std::endl;
        result = std::max(result, squaredDistance);
    }

    return result;
}

void meshkernel::SampleAveragingInterpolator::GenerateSearchPolygon(const double relativeSearchRadius,
                                                                    const Point& interpolationPoint,
                                                                    std::vector<Point>& polygon,
                                                                    const Projection projection) const
{
    std::ranges::transform(std::begin(polygon),
                           std::end(polygon),
                           std::begin(polygon),
                           [&](const Point& p)
                           { return p * relativeSearchRadius + interpolationPoint * (1.0 - relativeSearchRadius); });

    if (projection == Projection::spherical)
    {
        const auto boundingBox = BoundingBox(polygon);
        const auto lowerLeft = boundingBox.lowerLeft();
        const auto upperRight = boundingBox.upperRight();

        if (upperRight.x - lowerLeft.x <= 180.0)
        {
            return;
        }

        auto const x_mean = 0.5 * (upperRight.x + lowerLeft.x);

        for (auto& value : polygon)
        {
            if (value.x < x_mean)
            {
                value.x = value.x + 360.0;
            }
        }
    }
}

double meshkernel::SampleAveragingInterpolator::GetSampleValueFromRTree(const int propertyId, const UInt index) const
{
    UInt sampleIndex = m_nodeRTree->GetQueryResult(index);
    return m_sampleData.at(propertyId)[sampleIndex];
}

double meshkernel::SampleAveragingInterpolator::ComputeInterpolationResultFromNeighbors(const int propertyId,
                                                                                        const Point& interpolationPoint,
                                                                                        const std::vector<Point>& searchPolygon,
                                                                                        const Projection projection,
                                                                                        std::vector<Sample>& sampleCache) const
{
    sampleCache.clear();

    const std::vector<double>& propertyData(m_sampleData.at(propertyId));

    Point polygonCentre{constants::missing::doubleValue, constants::missing::doubleValue};

    if (projection == Projection::sphericalAccurate)
    {
        polygonCentre = ComputeAverageCoordinate(searchPolygon, projection);
    }

    for (UInt i = 0; i < m_nodeRTree->GetQueryResultSize(); ++i)
    {
        auto const sampleIndex = m_nodeRTree->GetQueryResult(i);
        auto const sampleValue = propertyData[sampleIndex];

        if (sampleValue == constants::missing::doubleValue)
        {
            continue;
        }

        const Point& samplePoint(m_samplePoints[sampleIndex]);

        if (IsPointInPolygonNodes(samplePoint, searchPolygon, projection, polygonCentre))
        {
            sampleCache.emplace_back(samplePoint.x, samplePoint.y, sampleValue);
        }
    }

    return m_strategy->Calculate(interpolationPoint, sampleCache);
}

double meshkernel::SampleAveragingInterpolator::ComputeOnPolygon(const int propertyId,
                                                                 std::vector<Point>& polygon,
                                                                 const Point& interpolationPoint,
                                                                 const double relativeSearchRadius,
                                                                 const bool useClosestSampleIfNoneAvailable,
                                                                 const Projection projection,
                                                                 std::vector<Sample>& sampleCache) const
{

    if (!interpolationPoint.IsValid())
    {
        throw ConstraintError("Invalid interpolation point");
    }

    GenerateSearchPolygon(relativeSearchRadius, interpolationPoint, polygon, projection);

    const double searchRadiusSquared = GetSearchRadiusSquared(polygon, interpolationPoint, projection);

    if (searchRadiusSquared <= 0.0)
    {
        std::cout << " interpolationPoint " << interpolationPoint.x << ", " << interpolationPoint.y << "  " << polygon.size() << std::endl;
        throw ConstraintError("Search radius: {} <= 0", searchRadiusSquared);
    }

    m_nodeRTree->SearchPoints(interpolationPoint, searchRadiusSquared);

    if (!m_nodeRTree->HasQueryResults() && useClosestSampleIfNoneAvailable)
    {
        m_nodeRTree->SearchNearestPoint(interpolationPoint);
        return m_nodeRTree->HasQueryResults() ? GetSampleValueFromRTree(propertyId, 0) : constants::missing::doubleValue;
    }
    else if (m_nodeRTree->HasQueryResults())
    {
        return ComputeInterpolationResultFromNeighbors(propertyId, interpolationPoint, polygon, projection, sampleCache);
    }

    return constants::missing::doubleValue;
}

void meshkernel::SampleAveragingInterpolator::InterpolateSpan(const int propertyId, const Mesh2D& mesh, const Location location, std::span<double>& result) const
{
    std::ranges::fill(result, constants::missing::doubleValue);

    std::cout << "SampleAveragingInterpolator::InterpolateSpan " << std::endl;

    if (location == Location::Nodes)
    {
        std::cout << " location nodes: " << mesh.GetNumNodes() << std::endl;
        std::vector<Point> dualFacePolygon;
        dualFacePolygon.reserve(MaximumNumberOfEdgesPerNode);
        std::vector<Sample> sampleCache;
        sampleCache.reserve(100);

        for (UInt n = 0; n < mesh.GetNumNodes(); ++n)
        {
            mesh.MakeDualFace(n, m_interpolationParameters.m_relativeSearchRadius, dualFacePolygon);

            double resultValue = constants::missing::doubleValue;

            if (dualFacePolygon.size() > 0)
            {
                resultValue = ComputeOnPolygon(propertyId,
                                               dualFacePolygon,
                                               mesh.Node(n),
                                               m_interpolationParameters.m_relativeSearchRadius,
                                               m_interpolationParameters.m_useClosestIfNoneFound,
                                               mesh.m_projection,
                                               sampleCache);
            }

            result[n] = resultValue;
            std::cout << "averaging interpoaltion: " << mesh.Node(n).x << "   " << mesh.Node(n).y << "   " << resultValue << std::endl;
        }
    }
}

double meshkernel::SampleAveragingInterpolator::InterpolateValue(const int propertyId [[maybe_unused]], const Point& evaluationPoint [[maybe_unused]]) const
{
    return constants::missing::doubleValue;
}
