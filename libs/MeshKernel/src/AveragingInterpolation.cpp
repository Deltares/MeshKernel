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

#include <MeshKernel/AveragingInterpolation.hpp>
#include <MeshKernel/AveragingStrategies/AveragingStrategyFactory.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/RTree.hpp>

using meshkernel::AveragingInterpolation;

AveragingInterpolation::AveragingInterpolation(std::shared_ptr<Mesh2D> mesh,
                                               std::vector<Sample>& samples,
                                               Method method,
                                               Mesh::Location locationType,
                                               double relativeSearchRadius,
                                               bool useClosestSampleIfNoneAvailable,
                                               bool transformSamples,
                                               size_t minNumSamples)
    : m_mesh(mesh),
      m_samples(samples),
      m_method(method),
      m_interpolationLocation(locationType),
      m_relativeSearchRadius(relativeSearchRadius),
      m_useClosestSampleIfNoneAvailable(useClosestSampleIfNoneAvailable),
      m_transformSamples(transformSamples),
      m_minNumSamples(minNumSamples)
{
}

void AveragingInterpolation::Compute()
{
    if (m_samples.empty())
    {
        throw AlgorithmError("AveragingInterpolation::Compute: No samples available.");
    }

    m_visitedSamples.resize(m_samples.size());

    // build sample r-tree for searches
    m_samplesRtree.BuildTree(m_samples);

    if (m_interpolationLocation == Mesh::Location::Nodes || m_interpolationLocation == Mesh::Location::Edges)
    {
        std::vector<Point> dualFacePolygon;
        m_results.resize(m_mesh->GetNumNodes(), constants::missing::doubleValue);
        std::fill(m_results.begin(), m_results.end(), constants::missing::doubleValue);

        // make sure edge centers are computed
        m_mesh->ComputeEdgesCenters();
        for (size_t n = 0; n < m_mesh->GetNumNodes(); ++n)
        {
            m_mesh->MakeDualFace(n, m_relativeSearchRadius, dualFacePolygon);

            const auto result = ComputeOnPolygon(dualFacePolygon, m_mesh->m_nodes[n]);

            // flag the visited samples
            for (size_t i = 0; i < m_samplesRtree.GetQueryResultSize(); ++i)
            {
                const auto sample = m_samplesRtree.GetQueryResult(i);
                m_visitedSamples[sample] = true;
            }

            m_results[n] = result;
        }
    }

    // for edges, an average of the nodal interpolated value is made
    if (m_interpolationLocation == Mesh::Location::Edges)
    {
        if (m_results.size() != m_mesh->GetNumNodes())
        {
            throw std::runtime_error("Number of nodes not matching!");
        }

        const auto nodeResults = m_results;
        m_results.resize(m_mesh->GetNumEdges(), constants::missing::doubleValue);
        std::fill(m_results.begin(), m_results.end(), constants::missing::doubleValue);

        for (size_t e = 0; e < m_mesh->GetNumEdges(); ++e)
        {
            const auto first = m_mesh->m_edges[e].first;
            const auto second = m_mesh->m_edges[e].second;

            const auto firstValue = nodeResults[first];
            const auto secondValue = nodeResults[second];

            if (!IsEqual(firstValue, constants::missing::doubleValue) && !IsEqual(secondValue, constants::missing::doubleValue))
            {
                m_results[e] = (firstValue + secondValue) * 0.5;
            }
        }
    }

    if (m_interpolationLocation == Mesh::Location::Faces)
    {
        m_results.resize(m_mesh->GetNumFaces(), constants::missing::doubleValue);
        std::fill(m_results.begin(), m_results.end(), constants::missing::doubleValue);

        std::vector<Point> polygonNodesCache(Mesh::m_maximumNumberOfNodesPerFace + 1);
        std::fill(m_visitedSamples.begin(), m_visitedSamples.end(), false);
        for (size_t f = 0; f < m_mesh->GetNumFaces(); ++f)
        {
            polygonNodesCache.clear();
            const auto numFaceNodes = m_mesh->GetNumFaceEdges(f);

            for (size_t n = 0; n < numFaceNodes; ++n)
            {
                polygonNodesCache.emplace_back(m_mesh->m_facesMassCenters[f] + (m_mesh->m_nodes[m_mesh->m_facesNodes[f][n]] - m_mesh->m_facesMassCenters[f]) * m_relativeSearchRadius);
            }
            polygonNodesCache.emplace_back(polygonNodesCache[0]);

            m_results[f] = ComputeOnPolygon(polygonNodesCache, m_mesh->m_facesMassCenters[f]);

            if (m_transformSamples && m_results[f] > 0)
            {
                // for certain algorithms we want to decrease the values of the samples (e.g. refinement)
                // it is difficult to do it otherwise without sharing or caching the query result
                for (size_t i = 0; i < m_samplesRtree.GetQueryResultSize(); ++i)
                {
                    if (const auto sample = m_samplesRtree.GetQueryResult(i); !m_visitedSamples[sample])
                    {
                        m_visitedSamples[sample] = true;
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
    std::transform(std::begin(polygon),
                   std::end(polygon),
                   begin(searchPolygon),
                   [&](Point const& p)
                   { return p * m_relativeSearchRadius + interpolationPoint * (1.0 - m_relativeSearchRadius); });

    if (m_mesh->m_projection == Projection::spherical)
    {
        auto [lowerLeft, upperRight] = GetBoundingBox(searchPolygon);

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
        auto const squaredDistance = ComputeSquaredDistance(interpolationPoint, value, m_mesh->m_projection);
        result = std::max(result, squaredDistance);
    }

    return result;
}

double AveragingInterpolation::GetSampleValueFromRTree(size_t const index)
{
    auto const sample_index = m_samplesRtree.GetQueryResult(index);
    return m_samples[sample_index].value;
}

double AveragingInterpolation::ComputeInterpolationResultFromNeighbors(std::unique_ptr<averaging::AveragingStrategy> strategy,
                                                                       std::vector<Point> const& searchPolygon)
{

    for (size_t i = 0; i < m_samplesRtree.GetQueryResultSize(); ++i)
    {
        auto const sampleIndex = m_samplesRtree.GetQueryResult(i);
        auto const sampleValue = m_samples[sampleIndex].value;

        if (sampleValue <= constants::missing::doubleValue)
        {
            continue;
        }

        Point samplePoint{m_samples[sampleIndex].x, m_samples[sampleIndex].y};
        if (IsPointInPolygonNodes(samplePoint, searchPolygon, m_mesh->m_projection))
        {
            strategy->Add(samplePoint, sampleValue);
        }
    }

    return strategy->Calculate();
}

double AveragingInterpolation::ComputeOnPolygon(const std::vector<Point>& polygon,
                                                Point const interpolationPoint)
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

    m_samplesRtree.SearchPoints(interpolationPoint, searchRadiusSquared);
    if (!m_samplesRtree.HasQueryResults() && m_useClosestSampleIfNoneAvailable)
    {
        m_samplesRtree.SearchNearestPoint(interpolationPoint);
        return m_samplesRtree.HasQueryResults() ? GetSampleValueFromRTree(0) : constants::missing::doubleValue;
    }
    if (m_samplesRtree.HasQueryResults())
    {

        auto strategy = averaging::AveragingStrategyFactory::GetAveragingStrategy(m_method, m_minNumSamples, interpolationPoint, m_mesh->m_projection);
        return ComputeInterpolationResultFromNeighbors(std::move(strategy), searchPolygon);
    }

    return constants::missing::doubleValue;
}
