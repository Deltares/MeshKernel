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

#include <stdexcept>
#include <tuple>

#include <MeshKernel/AveragingInterpolation.hpp>
#include <MeshKernel/AveragingStrategies/AveragingStrategyFactory.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/RTree.hpp>

namespace meshkernel
{

    AveragingInterpolation::AveragingInterpolation(std::shared_ptr<Mesh2D> mesh,
                                                   std::vector<Sample>& samples,
                                                   Method method,
                                                   MeshLocations locationType,
                                                   double relativeSearchRadius,
                                                   bool useClosestSampleIfNoneAvailable,
                                                   bool transformSamples) : m_mesh(mesh),
                                                                            m_samples(samples),
                                                                            m_method(method),
                                                                            m_interpolationLocation(locationType),
                                                                            m_relativeSearchRadius(relativeSearchRadius),
                                                                            m_useClosestSampleIfNoneAvailable(useClosestSampleIfNoneAvailable),
                                                                            m_transformSamples(transformSamples)
    {
    }

    void AveragingInterpolation::Compute()
    {
        if (m_samples.empty())
        {
            throw AlgorithmError("TriangulationInterpolation::Compute: No samples available.");
        }

        m_visitedSamples.resize(m_samples.size());

        // build sample r-tree for searches
        m_samplesRtree.BuildTree(m_samples);

        auto interpolatedResults = ComputeOnLocations();

        //for edges, an average of the nodal interpolated value is made
        if (m_interpolationLocation == MeshLocations::Edges)
        {
            m_results.resize(m_mesh->GetNumEdges(), doubleMissingValue);
            for (auto e = 0; e < m_mesh->GetNumEdges(); ++e)
            {
                const auto first = m_mesh->m_edges[e].first;
                const auto second = m_mesh->m_edges[e].second;

                const auto firstValue = interpolatedResults[first];
                const auto secondValue = interpolatedResults[second];

                if (!IsEqual(firstValue, doubleMissingValue) && !IsEqual(secondValue, doubleMissingValue))
                {
                    m_results[e] = (firstValue + secondValue) * 0.5;
                }
                else if (IsEqual(firstValue, doubleMissingValue))
                {
                    m_results[e] = secondValue;
                }
                else if (IsEqual(secondValue, doubleMissingValue))
                {
                    m_results[e] = firstValue;
                }
            }
            return;
        }

        //for the other cases, the interpolated values are already at the correct location
        m_results = std::move(interpolatedResults);
    }
    std::vector<double> AveragingInterpolation::ComputeOnFaces()
    {
        std::vector<double> interpolatedResults(m_mesh->GetNumFaces(), doubleMissingValue);
        std::vector<Point> polygonNodesCache(maximumNumberOfNodesPerFace + 1);
        std::fill(m_visitedSamples.begin(), m_visitedSamples.end(), false);
        for (auto f = 0; f < m_mesh->GetNumFaces(); ++f)
        {
            polygonNodesCache.clear();
            const auto numFaceNodes = m_mesh->GetNumFaceEdges(f);

            for (auto n = 0; n < numFaceNodes; ++n)
            {
                polygonNodesCache.emplace_back(m_mesh->m_facesMassCenters[f] + (m_mesh->m_nodes[m_mesh->m_facesNodes[f][n]] - m_mesh->m_facesMassCenters[f]) * m_relativeSearchRadius);
            }
            polygonNodesCache.emplace_back(polygonNodesCache[0]);

            interpolatedResults[f] = ComputeOnPolygon(polygonNodesCache, m_mesh->m_facesMassCenters[f]);

            if (m_transformSamples && interpolatedResults[f] > 0)
            {
                // for certain algorithms we want to decrease the values of the samples (e.g. refinement)
                // it is difficult to do it otherwise without sharing or caching the query result
                for (auto i = 0; i < m_samplesRtree.GetQueryResultSize(); i++)
                {
                    if (const auto sample = m_samplesRtree.GetQueryResult(i); !m_visitedSamples[sample])
                    {
                        m_visitedSamples[sample] = true;
                        m_samples[sample].value -= 1;
                    }
                }
            }
        }
        return interpolatedResults;
    }

    std::vector<double> AveragingInterpolation::ComputeOnNodesOrEdges()
    {
        std::vector<Point> dualFacePolygon;
        std::vector<double> interpolatedResults(m_mesh->GetNumNodes(), doubleMissingValue);
        // make sure edge centers are computed
        m_mesh->ComputeEdgesCenters();
        for (auto n = 0; n < m_mesh->GetNumNodes(); ++n)
        {
            m_mesh->MakeDualFace(n, m_relativeSearchRadius, dualFacePolygon);

            const auto result = ComputeOnPolygon(dualFacePolygon, m_mesh->m_nodes[n]);

            // flag the visited samples
            for (auto i = 0; i < m_samplesRtree.GetQueryResultSize(); i++)
            {
                const auto sample = m_samplesRtree.GetQueryResult(i);
                m_visitedSamples[sample] = true;
            }

            interpolatedResults[n] = result;
        }
        return interpolatedResults;
    }

    std::vector<double> AveragingInterpolation::ComputeOnLocations()
    {
        switch (m_interpolationLocation)
        {
        case MeshLocations::Faces:
        {
            return ComputeOnFaces();
        }
        case MeshLocations::Nodes:
        case MeshLocations::Edges:
        {
            return ComputeOnNodesOrEdges();
        }
        default:
            return std::vector<double>{};
        }
    }

    void ValidateInterpolationPoint(Point const& interpolationPoint)
    {
        if (!interpolationPoint.IsValid())
        {
            throw std::invalid_argument("AveragingInterpolation::ComputeOnPolygon invalid interpolation point");
        }
    }

    void AdjustForSpherical(std::vector<Point>& searchPolygon)
    {
        auto [lowerLeft, upperRight] = GetBoundingBox(searchPolygon);

        if (upperRight.x - lowerLeft.x <= 180.0)
            return;

        auto const x_mean = 0.5 * (upperRight.x + lowerLeft.x);

        lowerLeft.x = std::numeric_limits<double>::max();
        upperRight.x = std::numeric_limits<double>::lowest();

        for (auto& value : searchPolygon)
        {
            if (value.x < x_mean)
            {
                value.x = value.x + 360.0;
                lowerLeft.x = std::min(lowerLeft.x, value.x);
                upperRight.x = std::max(upperRight.x, value.x);
            }
        }
    }

    std::vector<Point> AveragingInterpolation::GetSearchPolygon(std::vector<Point> const& polygon, Point const& interpolationPoint)
    {
        std::vector<Point> searchPolygon(polygon.size());
        std::transform(std::begin(polygon),
                       std::end(polygon),
                       begin(searchPolygon),
                       [&](Point const& p) { return p * m_relativeSearchRadius + interpolationPoint * (1.0 - m_relativeSearchRadius); });

        if (m_mesh->m_projection == Projection::spherical)
        {
            AdjustForSpherical(searchPolygon);
        }

        return searchPolygon;
    }

    double getSearchRadiusSquared(std::vector<Point> const& searchPolygon,
                                  Point const& interpolationPoint,
                                  Projection meshProjection)
    {
        double result = std::numeric_limits<double>::lowest();

        for (const auto& value : searchPolygon)
        {
            auto const squaredDistance = ComputeSquaredDistance(interpolationPoint, value, meshProjection);
            result = std::max(result, squaredDistance);
        }

        return result;
    }

    double AveragingInterpolation::GetSampleValueFromRTree(size_t const index)
    {
        auto const sample_index = m_samplesRtree.GetQueryResult(index);
        return m_samples[sample_index].value;
    }

    double AveragingInterpolation::ComputeInterpolationResultFromNeighbors(std::unique_ptr<averaging::AveragingStrategy> strategy, std::vector<Point> const& searchPolygon)
    {

        for (auto i = 0; i < m_samplesRtree.GetQueryResultSize(); i++)
        {
            auto const sampleIndex = m_samplesRtree.GetQueryResult(i);
            auto const sampleValue = m_samples[sampleIndex].value;

            if (sampleValue <= doubleMissingValue)
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

        std::vector<Point> const searchPolygon = GetSearchPolygon(polygon, interpolationPoint);

        double const searchRadiusSquared = getSearchRadiusSquared(searchPolygon, interpolationPoint, m_mesh->m_projection);

        if (searchRadiusSquared <= 0.0)
        {
            throw std::invalid_argument("AveragingInterpolation::ComputeOnPolygon search radius <= 0");
        }

        m_samplesRtree.NearestNeighborsOnSquaredDistance(interpolationPoint, searchRadiusSquared);
        if (!m_samplesRtree.HasQueryResults() && m_useClosestSampleIfNoneAvailable)
        {
            m_samplesRtree.NearestNeighbors(interpolationPoint);
            return m_samplesRtree.HasQueryResults() ? GetSampleValueFromRTree(0) : doubleMissingValue;
        }
        if (m_samplesRtree.HasQueryResults())
        {

            auto strategy = averaging::AveragingStrategyFactory::GetAveragingStrategy(m_method, interpolationPoint, m_mesh->m_projection);
            return ComputeInterpolationResultFromNeighbors(std::move(strategy), searchPolygon);
        }

        return doubleMissingValue;
    }
} // namespace meshkernel
