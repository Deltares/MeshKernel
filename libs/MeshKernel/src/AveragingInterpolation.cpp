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

#include "MeshKernel/Utilities/RTreeFactory.hpp"
#include <MeshKernel/AveragingInterpolation.hpp>
#include <MeshKernel/AveragingStrategies/AveragingStrategyFactory.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Operations.hpp>

using meshkernel::AveragingInterpolation;

AveragingInterpolation::AveragingInterpolation(Mesh2D& mesh,
                                               std::vector<Sample>& samples,
                                               Method method,
                                               Location locationType,
                                               double relativeSearchRadius,
                                               bool useClosestSampleIfNoneAvailable,
                                               bool transformSamples,
                                               UInt minNumSamples)
    : m_mesh(mesh),
      m_samples(samples),
      m_interpolationLocation(locationType),
      m_relativeSearchRadius(relativeSearchRadius),
      m_useClosestSampleIfNoneAvailable(useClosestSampleIfNoneAvailable),
      m_transformSamples(transformSamples),
      m_samplesRtree(RTreeFactory::Create(mesh.m_projection)),
      m_strategy(averaging::AveragingStrategyFactory::GetAveragingStrategy(method, minNumSamples, m_mesh.m_projection))
{
    m_interpolationSampleCache.reserve(DefaultMaximumCacheSize);
}

void AveragingInterpolation::Compute()
{
    if (m_samples.empty())
    {
        throw AlgorithmError("AveragingInterpolation::Compute: No samples available.");
    }

    if (m_samplesRtree->Empty())
    {
        m_samplesRtree->BuildTree(m_samples);
    }

    if (m_interpolationLocation == Location::Nodes || m_interpolationLocation == Location::Edges)
    {
        m_nodeResults.resize(m_mesh.GetNumNodes(), constants::missing::doubleValue);
        std::ranges::fill(m_nodeResults, constants::missing::doubleValue);

        // make sure edge centers are computed
        m_mesh.ComputeEdgesCenters();

        std::vector<Point> dualFacePolygon;
        for (UInt n = 0; n < m_mesh.GetNumNodes(); ++n)
        {
            m_mesh.MakeDualFace(n, m_relativeSearchRadius, dualFacePolygon);
            const auto result = ComputeOnPolygon(dualFacePolygon, m_mesh.m_nodes[n]);
            m_nodeResults[n] = result;
        }
    }

    // for edges, an average of the nodal interpolated value is made
    if (m_interpolationLocation == Location::Edges)
    {
        m_edgeResults.resize(m_mesh.GetNumEdges(), constants::missing::doubleValue);
        std::ranges::fill(m_edgeResults, constants::missing::doubleValue);

        for (UInt e = 0; e < m_mesh.GetNumEdges(); ++e)
        {
            const auto& [first, second] = m_mesh.m_edges[e];

            const auto& firstValue = m_nodeResults[first];
            const auto& secondValue = m_nodeResults[second];

            if (!IsEqual(firstValue, constants::missing::doubleValue) && !IsEqual(secondValue, constants::missing::doubleValue))
            {
                m_edgeResults[e] = 0.5 * (firstValue + secondValue);
            }
        }
    }

    if (m_interpolationLocation == Location::Faces)
    {
        std::vector<bool> visitedSamples(m_samples.size(), false); ///< The visited samples
        std::vector<Point> polygonNodesCache(Mesh::m_maximumNumberOfNodesPerFace + 1);
        m_faceResults.resize(m_mesh.GetNumFaces(), constants::missing::doubleValue);
        std::ranges::fill(m_faceResults, constants::missing::doubleValue);

        for (UInt f = 0; f < m_mesh.GetNumFaces(); ++f)
        {
            polygonNodesCache.clear();

            for (UInt n = 0; n < m_mesh.GetNumFaceEdges(f); ++n)
            {
                polygonNodesCache.emplace_back(m_mesh.m_facesMassCenters[f] + (m_mesh.m_nodes[m_mesh.m_facesNodes[f][n]] - m_mesh.m_facesMassCenters[f]) * m_relativeSearchRadius);
            }
            polygonNodesCache.emplace_back(polygonNodesCache[0]);

            m_faceResults[f] = ComputeOnPolygon(polygonNodesCache, m_mesh.m_facesMassCenters[f]);

            if (m_transformSamples && m_faceResults[f] > 0)
            {
                // for certain algorithms we want to decrease the values of the samples (e.g. refinement)
                // it is difficult to do it otherwise without sharing or caching the query result
                for (UInt i = 0; i < m_samplesRtree->GetQueryResultSize(); ++i)
                {
                    if (const auto sample = m_samplesRtree->GetQueryResult(i); !visitedSamples[sample])
                    {
                        visitedSamples[sample] = true;
                        m_samples[sample].value -= 1;
                    }
                }
            }
        }
    }
}

std::vector<meshkernel::Point> AveragingInterpolation::GetSearchPolygon(std::vector<Point> const& polygon, Point const& interpolationPoint) const
{
    std::vector<Point> searchPolygon(polygon.size());
    std::ranges::transform(std::begin(polygon),
                           std::end(polygon),
                           begin(searchPolygon),
                           [&](Point const& p)
                           { return p * m_relativeSearchRadius + interpolationPoint * (1.0 - m_relativeSearchRadius); });

    if (m_mesh.m_projection == Projection::spherical)
    {
        const auto boundingBox = BoundingBox(searchPolygon);
        const auto lowerLeft = boundingBox.lowerLeft();
        const auto upperRight = boundingBox.upperRight();

        if (upperRight.x - lowerLeft.x <= 180.0)
            return searchPolygon;

        auto const x_mean = 0.5 * (upperRight.x + lowerLeft.x);

        for (auto& value : searchPolygon)
        {
            if (value.x < x_mean)
            {
                value.x = value.x + 360.0;
            }
        }
    }

    return searchPolygon;
}

double AveragingInterpolation::GetSearchRadiusSquared(std::vector<Point> const& searchPolygon,
                                                      Point const& interpolationPoint) const
{
    double result = std::numeric_limits<double>::lowest();

    for (const auto& value : searchPolygon)
    {
        auto const squaredDistance = ComputeSquaredDistance(interpolationPoint, value, m_mesh.m_projection);
        result = std::max(result, squaredDistance);
    }

    return result;
}

double AveragingInterpolation::GetSampleValueFromRTree(UInt const index)
{
    auto sample_index = m_samplesRtree->GetQueryResult(index);
    return m_samples[sample_index].value;
}

double AveragingInterpolation::ComputeInterpolationResultFromNeighbors(const Point& interpolationPoint,
                                                                       std::vector<Point> const& searchPolygon)
{
    m_interpolationSampleCache.clear();

    for (UInt i = 0; i < m_samplesRtree->GetQueryResultSize(); ++i)
    {
        auto const sampleIndex = m_samplesRtree->GetQueryResult(i);
        auto const sampleValue = m_samples[sampleIndex].value;

        if (sampleValue <= constants::missing::doubleValue)
        {
            continue;
        }

        Point samplePoint{m_samples[sampleIndex].x, m_samples[sampleIndex].y};

        if (IsPointInPolygonNodes(samplePoint, searchPolygon, m_mesh.m_projection))
        {
            m_interpolationSampleCache.emplace_back(samplePoint.x, samplePoint.y, sampleValue);
        }
    }

    return m_strategy->Calculate(interpolationPoint, m_interpolationSampleCache);
}

double AveragingInterpolation::ComputeOnPolygon(const std::vector<Point>& polygon,
                                                const Point& interpolationPoint)
{

    if (!interpolationPoint.IsValid())
    {
        throw std::invalid_argument("AveragingInterpolation::ComputeOnPolygon invalid interpolation point");
    }

    auto const searchPolygon = GetSearchPolygon(polygon, interpolationPoint);

    double const searchRadiusSquared = GetSearchRadiusSquared(searchPolygon, interpolationPoint);

    if (searchRadiusSquared <= 0.0)
    {
        throw std::invalid_argument("AveragingInterpolation::ComputeOnPolygon search radius <= 0");
    }

    m_samplesRtree->SearchPoints(interpolationPoint, searchRadiusSquared);
    if (!m_samplesRtree->HasQueryResults() && m_useClosestSampleIfNoneAvailable)
    {
        m_samplesRtree->SearchNearestPoint(interpolationPoint);
        return m_samplesRtree->HasQueryResults() ? GetSampleValueFromRTree(0) : constants::missing::doubleValue;
    }
    if (m_samplesRtree->HasQueryResults())
    {

        return ComputeInterpolationResultFromNeighbors(interpolationPoint, searchPolygon);
    }

    return constants::missing::doubleValue;
}
