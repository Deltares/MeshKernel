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

#include "MeshKernel/SampleAveragingInterpolator.hpp"
#include "MeshKernel/Operations.hpp"

std::vector<meshkernel::Point> meshkernel::SampleAveragingInterpolator::CombineCoordinates(const std::span<const double> xNodes,
                                                                                           const std::span<const double> yNodes)
{
    std::vector<Point> result(xNodes.size());

    for (size_t i = 0; i < xNodes.size(); ++i)
    {
        result[i].x = xNodes[i];
        result[i].y = yNodes[i];
    }

    return result;
}

meshkernel::SampleAveragingInterpolator::SampleAveragingInterpolator(const std::span<const double> xNodes,
                                                                     const std::span<const double> yNodes,
                                                                     const Projection projection,
                                                                     const InterpolationParameters& interpolationParameters)
    : m_samplePoints(CombineCoordinates(xNodes, yNodes)),
      m_projection(projection),
      m_interpolationParameters(interpolationParameters),
      m_strategy(averaging::AveragingStrategyFactory::GetAveragingStrategy(interpolationParameters.m_method,
                                                                           interpolationParameters.m_minimumNumberOfSamples,
                                                                           projection)),
      m_nodeRTree(RTreeFactory::Create(projection))
{
    m_nodeRTree->BuildTree(m_samplePoints);
}

meshkernel::SampleAveragingInterpolator::SampleAveragingInterpolator(const std::span<const Point> nodes,
                                                                     const Projection projection,
                                                                     const InterpolationParameters& interpolationParameters)
    : m_samplePoints(nodes.begin(), nodes.end()),
      m_projection(projection),
      m_interpolationParameters(interpolationParameters),
      m_strategy(averaging::AveragingStrategyFactory::GetAveragingStrategy(interpolationParameters.m_method,
                                                                           interpolationParameters.m_minimumNumberOfSamples,
                                                                           projection)),
      m_nodeRTree(RTreeFactory::Create(projection))
{
    m_nodeRTree->BuildTree(m_samplePoints);
}

void meshkernel::SampleAveragingInterpolator::Interpolate(const int propertyId, const std::span<const Point> interpolationNodes, std::span<double> result) const
{
    const std::vector<double>& propertyValues = GetSampleData(propertyId);
    std::vector<Sample> sampleCache;
    sampleCache.reserve(100);

    double searchRadiusSquared = 1.0e5;

    for (size_t i = 0; i < interpolationNodes.size(); ++i)
    {
        m_nodeRTree->SearchPoints(interpolationNodes[i], searchRadiusSquared);

        double resultValue = constants::missing::doubleValue;

        if (!m_nodeRTree->HasQueryResults())
        {
            resultValue = constants::missing::doubleValue;
        }
        else
        {
            sampleCache.clear();

            for (UInt j = 0; j < m_nodeRTree->GetQueryResultSize(); ++j)
            {
                auto const sampleIndex = m_nodeRTree->GetQueryResult(j);
                auto const sampleValue = propertyValues[sampleIndex];

                if (sampleValue == constants::missing::doubleValue)
                {
                    continue;
                }

                const Point& samplePoint(m_samplePoints[sampleIndex]);

                sampleCache.emplace_back(samplePoint.x, samplePoint.y, sampleValue);
            }

            resultValue = m_strategy->Calculate(interpolationNodes[i], sampleCache);
        }

        result[i] = resultValue;
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
    return GetSampleData(propertyId)[sampleIndex];
}

double meshkernel::SampleAveragingInterpolator::ComputeInterpolationResultFromNeighbors(const int propertyId,
                                                                                        const Point& interpolationPoint,
                                                                                        const std::vector<Point>& searchPolygon,
                                                                                        const Projection projection,
                                                                                        std::vector<Sample>& sampleCache) const
{
    sampleCache.clear();

    const std::vector<double>& propertyData(GetSampleData(propertyId));

    Point polygonCentre{constants::missing::doubleValue, constants::missing::doubleValue};

    if (projection == Projection::sphericalAccurate)
    {
        polygonCentre = ComputeAverageCoordinate(searchPolygon, projection);
    }

    BoundingBox boundingBox(searchPolygon);

    for (UInt i = 0; i < m_nodeRTree->GetQueryResultSize(); ++i)
    {
        auto const sampleIndex = m_nodeRTree->GetQueryResult(i);
        auto const sampleValue = propertyData[sampleIndex];

        if (sampleValue == constants::missing::doubleValue)
        {
            continue;
        }

        const Point& samplePoint(m_samplePoints[sampleIndex]);

        if (IsPointInPolygonNodes(samplePoint, searchPolygon, projection, boundingBox, polygonCentre))
        {
            sampleCache.emplace_back(samplePoint.x, samplePoint.y, sampleValue);
        }
    }

    return m_strategy->Calculate(interpolationPoint, sampleCache);
}

double meshkernel::SampleAveragingInterpolator::ComputeOnPolygon(const int propertyId,
                                                                 std::vector<Point>& polygon,
                                                                 const Point& interpolationPoint,
                                                                 const Projection projection,
                                                                 std::vector<Sample>& sampleCache) const
{

    if (!interpolationPoint.IsValid())
    {
        throw ConstraintError("Invalid interpolation point");
    }

    GenerateSearchPolygon(m_interpolationParameters.m_relativeSearchRadius, interpolationPoint, polygon, projection);

    const double searchRadiusSquared = GetSearchRadiusSquared(polygon, interpolationPoint, projection);

    if (searchRadiusSquared <= 0.0)
    {
        throw ConstraintError("Search radius: {} <= 0", searchRadiusSquared);
    }

    m_nodeRTree->SearchPoints(interpolationPoint, searchRadiusSquared);

    if (!m_nodeRTree->HasQueryResults() && m_interpolationParameters.m_useClosestIfNoneFound)
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

void meshkernel::SampleAveragingInterpolator::InterpolateAtNodes(const int propertyId, const Mesh2D& mesh, std::span<double>& result) const
{
    std::vector<Point> dualFacePolygon;
    dualFacePolygon.reserve(MaximumNumberOfEdgesPerNode);
    std::vector<Sample> sampleCache;
    sampleCache.reserve(100);

    for (UInt n = 0; n < mesh.GetNumNodes(); ++n)
    {
        mesh.MakeDualFace(n, m_interpolationParameters.m_relativeSearchRadius, dualFacePolygon);

        double resultValue = constants::missing::doubleValue;

        if (!dualFacePolygon.empty())
        {
            resultValue = ComputeOnPolygon(propertyId,
                                           dualFacePolygon,
                                           mesh.Node(n),
                                           mesh.m_projection,
                                           sampleCache);
        }

        result[n] = resultValue;
    }
}

void meshkernel::SampleAveragingInterpolator::InterpolateAtEdgeCentres(const Mesh2D& mesh,
                                                                       const std::span<double>& nodeResult,
                                                                       std::span<double>& result) const
{
    std::ranges::fill(result, constants::missing::doubleValue);

    for (UInt e = 0; e < mesh.GetNumEdges(); ++e)
    {
        const auto& [first, second] = mesh.GetEdge(e);

        if (first == constants::missing::uintValue || second == constants::missing::uintValue)
        {
            continue;
        }

        const double& firstValue = nodeResult[first];
        const double& secondValue = nodeResult[second];

        if (firstValue != constants::missing::doubleValue && secondValue != constants::missing::doubleValue)
        {
            result[e] = 0.5 * (firstValue + secondValue);
        }
    }
}

void meshkernel::SampleAveragingInterpolator::InterpolateAtFaces(const int propertyId, const Mesh2D& mesh, std::span<double>& result) const
{
    std::vector<Point> polygonNodesCache(MaximumNumberOfEdgesPerNode + 1);
    std::ranges::fill(result, constants::missing::doubleValue);
    std::vector<Sample> sampleCache;
    sampleCache.reserve(100);

    for (UInt f = 0; f < mesh.GetNumFaces(); ++f)
    {
        polygonNodesCache.clear();

        for (UInt n = 0; n < mesh.GetNumFaceEdges(f); ++n)
        {
            polygonNodesCache.emplace_back(mesh.m_facesMassCenters[f] + (mesh.Node(mesh.m_facesNodes[f][n]) - mesh.m_facesMassCenters[f]) * m_interpolationParameters.m_relativeSearchRadius);
        }

        // Close the polygon
        polygonNodesCache.emplace_back(polygonNodesCache[0]);

        sampleCache.clear();
        result[f] = ComputeOnPolygon(propertyId,
                                     polygonNodesCache,
                                     mesh.m_facesMassCenters[f],
                                     mesh.m_projection,
                                     sampleCache);
    }
}

void meshkernel::SampleAveragingInterpolator::Interpolate(const int propertyId, const Mesh2D& mesh, const Location location, std::span<double> result) const
{
    std::ranges::fill(result, constants::missing::doubleValue);

    // Only allocate intermediateNodeResult if needed, i.e. when location is edges
    std::vector<double> intermediateNodeResult;
    std::span<double> nodeResult;

    using enum Location;

    if (location == Nodes)
    {
        nodeResult = std::span<double>(result.data(), result.size());
    }
    else if (location == Edges)
    {
        intermediateNodeResult.resize(mesh.GetNumNodes(), 0.0);
        nodeResult = std::span<double>(intermediateNodeResult.data(), intermediateNodeResult.size());
    }

    if (location == Nodes || location == Edges)
    {
        InterpolateAtNodes(propertyId, mesh, nodeResult);
    }

    if (location == Edges)
    {
        InterpolateAtEdgeCentres(mesh, nodeResult, result);
    }

    if (location == Faces)
    {
        InterpolateAtFaces(propertyId, mesh, result);
    }
}

double meshkernel::SampleAveragingInterpolator::InterpolateValue(const int propertyId [[maybe_unused]], const Point& evaluationPoint [[maybe_unused]]) const
{
    return constants::missing::doubleValue;
}
