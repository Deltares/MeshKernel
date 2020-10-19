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

#include "Operations.cpp"
#include "Mesh.hpp"
#include "SpatialTrees.hpp"
#include "Averaging.hpp"

meshkernel::Averaging::Averaging(std::shared_ptr<Mesh> mesh,
                                 std::vector<Sample>& samples,
                                 AveragingMethod averagingMethod,
                                 InterpolationLocation locationType,
                                 double relativeSearchRadius,
                                 bool useClosestSampleIfNoneAvailable,
                                 bool transformSamples) : m_mesh(mesh),
                                                          m_samples(samples),
                                                          m_averagingMethod(averagingMethod),
                                                          m_interpolationLocation(locationType),
                                                          m_relativeSearchRadius(relativeSearchRadius),
                                                          m_useClosestSampleIfNoneAvailable(useClosestSampleIfNoneAvailable),
                                                          m_transformSamples(transformSamples)
{
    m_visitedSamples.resize(samples.size());
    // build sample rtree for searches
    m_samplesRtree.BuildTree(m_samples);
}

bool meshkernel::Averaging::Compute()
{

    std::fill(m_visitedSamples.begin(), m_visitedSamples.end(), false);

    if (m_interpolationLocation == Faces)
    {
        m_results.resize(m_mesh->GetNumFaces(), doubleMissingValue);
        std::vector<Point> polygonNodesCache(maximumNumberOfNodesPerFace + 1);

        for (int f = 0; f < m_mesh->GetNumFaces(); ++f)
        {
            polygonNodesCache.clear();
            const auto numFaceNodes = m_mesh->GetNumFaceEdges(f);

            for (int n = 0; n < numFaceNodes; ++n)
            {
                polygonNodesCache.push_back(m_mesh->m_facesMassCenters[f] + (m_mesh->m_nodes[m_mesh->m_facesNodes[f][n]] - m_mesh->m_facesMassCenters[f]) * m_relativeSearchRadius);
            }
            polygonNodesCache.push_back(polygonNodesCache[0]);

            double result = 0.0;
            auto successful = ComputeOnPolygon(polygonNodesCache, m_mesh->m_facesMassCenters[f], result);

            if (!successful)
            {
                return false;
            }

            m_results[f] = result;

            if (m_transformSamples && result > 0)
            {
                // for certain algorithms we want to decrease the values of the samples (e.g. refinement)
                // it is difficult to do it otherwise without sharing or caching the query result
                for (int i = 0; i < m_samplesRtree.GetQueryResultSize(); i++)
                {
                    const auto sample = m_samplesRtree.GetQuerySampleIndex(i);
                    if (!m_visitedSamples[sample])
                    {
                        m_visitedSamples[sample] = true;
                        m_samples[sample].value -= 1;
                    }
                }
            }
        }

        return true;
    }

    std::vector<double> interpolatedResults(m_mesh->GetNumNodes(), doubleMissingValue);
    if (m_interpolationLocation == Nodes || m_interpolationLocation == Edges)
    {
        std::vector<Point> dualFacePolygon;
        m_mesh->ComputeEdgesCenters();
        for (int n = 0; n < m_mesh->GetNumNodes(); ++n)
        {
            m_mesh->MakeDualFace(n, m_relativeSearchRadius, dualFacePolygon);

            double result = 0.0;
            const bool successful = ComputeOnPolygon(dualFacePolygon, m_mesh->m_nodes[n], result);

            // flag the visited samples
            for (int i = 0; i < m_samplesRtree.GetQueryResultSize(); i++)
            {
                const auto sample = m_samplesRtree.GetQuerySampleIndex(i);
                m_visitedSamples[sample] = true;
            }

            if (!successful)
            {
                return false;
            }

            interpolatedResults[n] = result;
        }
    }

    if (m_interpolationLocation == Edges)
    {
        m_results.resize(m_mesh->GetNumEdges(), doubleMissingValue);
        for (int e = 0; e < m_mesh->GetNumEdges(); ++e)
        {
            const auto first = m_mesh->m_edges[e].first;
            const auto second = m_mesh->m_edges[e].second;

            const auto firstValue = interpolatedResults[first];
            const auto secondValue = interpolatedResults[second];

            const bool isFirstValueValid = std::abs(firstValue - doubleMissingValue) > std::numeric_limits<double>::epsilon();
            const bool isSecondValueValid = std::abs(secondValue - doubleMissingValue) > std::numeric_limits<double>::epsilon();

            if (isFirstValueValid && isSecondValueValid)
            {
                m_results[e] = (firstValue + secondValue) * 0.5;
            }
            else if (!isFirstValueValid)
            {
                m_results[e] = secondValue;
            }
            else if (!isSecondValueValid)
            {
                m_results[e] = firstValue;
            }
        }
    }
    if (m_interpolationLocation == Nodes)
    {
        m_results = std::move(interpolatedResults);
    }

    return true;
}

bool meshkernel::Averaging::ComputeOnPolygon(const std::vector<Point>& polygon,
                                             Point interpolationPoint,
                                             double& result)
{

    if (!interpolationPoint.IsValid())
    {
        return false;
    }

    std::vector<Point> searchPolygon(polygon.size());

    // averaging settings
    double minx = std::numeric_limits<double>::max();
    double maxx = std::numeric_limits<double>::lowest();
    double miny = std::numeric_limits<double>::max();
    double maxy = std::numeric_limits<double>::lowest();

    for (int i = 0; i < polygon.size(); i++)
    {
        searchPolygon[i] = polygon[i] * 1.01 + interpolationPoint * (1.0 - 1.01);
        minx = std::min(minx, searchPolygon[i].x);
        maxx = std::max(maxx, searchPolygon[i].x);
        miny = std::min(miny, searchPolygon[i].y);
        maxy = std::max(maxy, searchPolygon[i].y);
    }

    if (m_mesh->m_projection == Projections::spherical && maxx - minx > 180.0)
    {

        double xmean = 0.5 * (maxx + minx);
        minx = std::numeric_limits<double>::max();
        maxx = std::numeric_limits<double>::lowest();
        for (int i = 0; i < polygon.size(); i++)
        {
            if (searchPolygon[i].x < xmean)
            {
                searchPolygon[i].x = searchPolygon[i].x + 360.0;
                minx = std::min(minx, searchPolygon[i].x);
                maxx = std::max(maxx, searchPolygon[i].x);
            }
        }
    }

    result = doubleMissingValue;
    double searchRadiusSquared = std::numeric_limits<double>::lowest();
    for (int i = 0; i < polygon.size() - 1; i++)
    {
        auto const squaredDistance = ComputeSquaredDistance(interpolationPoint, searchPolygon[i], m_mesh->m_projection);
        searchRadiusSquared = std::max(searchRadiusSquared, squaredDistance);
    }
    if (searchRadiusSquared <= 0.0)
    {
        return true;
    }

    m_samplesRtree.NearestNeighboursOnSquaredDistance(interpolationPoint, searchRadiusSquared);
    if (m_samplesRtree.GetQueryResultSize() == 0)
    {
        if (m_useClosestSampleIfNoneAvailable)
        {
            // use the closest sample if none available
            m_samplesRtree.NearestNeighbour(interpolationPoint);
            if (m_samplesRtree.GetQueryResultSize() > 0)
            {
                const auto sampleIndex = m_samplesRtree.GetQuerySampleIndex(0);
                if (sampleIndex >= 0)
                {
                    result = m_samples[sampleIndex].value;
                }
            }
        }
        return true;
    }

    int numValidSamplesInPolygon = 0;
    double wall = 0.0;
    bool firstValidSampleFound = false;
    double closestSquaredDistance = std::numeric_limits<double>::max();

    for (int i = 0; i < m_samplesRtree.GetQueryResultSize(); i++)
    {
        //do stuff based on the averaging method
        const auto sampleIndex = m_samplesRtree.GetQuerySampleIndex(i);
        const auto sampleValue = m_samples[sampleIndex].value;
        if (sampleValue <= doubleMissingValue)
        {
            continue;
        }

        Point samplePoint{m_samples[sampleIndex].x, m_samples[sampleIndex].y};
        // assume here polygon has a size equal to numPolygonNodes + 1
        bool isInPolygon = IsPointInPolygonNodes(samplePoint, searchPolygon, 0, int(searchPolygon.size() - 1));
        if (isInPolygon)
        {
            if (m_averagingMethod == SimpleAveraging)
            {
                if (!firstValidSampleFound)
                {
                    firstValidSampleFound = true;
                    result = 0.0;
                }
                result += sampleValue;
                numValidSamplesInPolygon++;
            }
            if (m_averagingMethod == ClosestPoint)
            {
                const auto squaredDistance = ComputeSquaredDistance(interpolationPoint, samplePoint, m_mesh->m_projection);
                if (squaredDistance < closestSquaredDistance)
                {
                    closestSquaredDistance = squaredDistance;
                    result = sampleValue;
                }
            }
            if (m_averagingMethod == Max)
            {
                if (!firstValidSampleFound)
                {
                    firstValidSampleFound = true;
                    result = std::numeric_limits<double>::lowest();
                }
                result = std::max(result, sampleValue);
            }
            if (m_averagingMethod == Min)
            {
                if (!firstValidSampleFound)
                {
                    firstValidSampleFound = true;
                    result = std::numeric_limits<double>::max();
                }
                result = std::min(result, sampleValue);
            }
            if (m_averagingMethod == InverseWeightDistance)
            {
                if (!firstValidSampleFound)
                {
                    firstValidSampleFound = true;
                    result = 0.0;
                }
                const double distance = std::max(0.01, Distance(interpolationPoint, samplePoint, m_mesh->m_projection));
                const double weight = 1.0 / distance;
                wall += weight;
                numValidSamplesInPolygon++;
                result += weight * sampleValue;
            }
            if (m_averagingMethod == MinAbsValue)
            {
                if (!firstValidSampleFound)
                {
                    firstValidSampleFound = true;
                    result = sampleValue;
                }
                result = std::min(std::abs(result), std::abs(sampleValue));
            }
        }
    }

    if (m_averagingMethod == SimpleAveraging && numValidSamplesInPolygon > 0)
    {
        if (result > doubleMissingValue)
        {
            result /= numValidSamplesInPolygon;
        }
    }

    if (m_averagingMethod == InverseWeightDistance && numValidSamplesInPolygon > 0)
    {
        result /= wall;
    }

    return true;
}