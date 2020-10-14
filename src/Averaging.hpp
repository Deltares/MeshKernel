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

namespace meshkernel
{
    // Forward declarations
    class Polygons;

    class Averaging
    {

    public:
        /// @brief
        /// @param mesh
        /// @param samples
        /// @param averagingMethod
        /// @param locationType
        /// @param relativeSearchRadius
        explicit Averaging(std::shared_ptr<Mesh> mesh, std::vector<Sample> samples, AveragingMethod averagingMethod, LocationTypes locationType, double relativeSearchRadius) : m_mesh(mesh),
                                                                                                                                                                                m_samples(samples),
                                                                                                                                                                                m_averagingMethod(averagingMethod),
                                                                                                                                                                                m_locationTypes(locationType),
                                                                                                                                                                                m_relativeSearchRadius(relativeSearchRadius)
        {
        }

        bool Compute()
        {

            m_samplesRtree.BuildTree(m_samples);
            if (m_locationTypes == CircumCenters)
            {
                m_results.resize(m_mesh->GetNumFaces(), doubleMissingValue);
                std::vector<Point> polygonNodesCache(maximumNumberOfNodesPerFace);

                for (int f = 0; f < m_mesh->GetNumFaces(); ++f)
                {
                    polygonNodesCache.clear();
                    const auto numFaceNodes = m_mesh->GetNumFaceEdges(f);
                    for (int n = 0; n < numFaceNodes; ++n)
                    {
                        polygonNodesCache.push_back(m_mesh->m_facesCircumcenters[f] + (m_mesh->m_nodes[m_mesh->m_facesNodes[f][n]] - m_mesh->m_facesCircumcenters[f]) * m_relativeSearchRadius);
                    }

                    double result = 0.0;
                    bool successful = AveragingOnPolygon(m_samples,
                                                         numFaceNodes,
                                                         polygonNodesCache,
                                                         m_mesh->m_facesCircumcenters[f],
                                                         m_mesh->m_projection,
                                                         m_samplesRtree,
                                                         m_averagingMethod,
                                                         result);

                    if (!successful)
                    {
                        return false;
                    }

                    m_results[f] = result;
                }
            }

            if (m_locationTypes == Nodes || m_locationTypes == Edges)
            {
                std::vector<Point> dualFacePolygon;
                std::vector<double> interpolatedResults(m_mesh->GetNumFaces(), doubleMissingValue);
                for (int n = 0; n < m_mesh->GetNumNodes(); ++n)
                {
                    m_mesh->MakeDualFace(n, m_relativeSearchRadius, dualFacePolygon);
                    double result = 0.0;
                    bool successful = AveragingOnPolygon(m_samples,
                                                         dualFacePolygon.size(),
                                                         dualFacePolygon,
                                                         m_mesh->m_nodes[n],
                                                         m_mesh->m_projection,
                                                         m_samplesRtree,
                                                         m_averagingMethod,
                                                         result);
                    if (!successful)
                    {
                        return false;
                    }

                    interpolatedResults[n] = result;
                }

                if (m_locationTypes == Edges)
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
                        if (isFirstValueValid)
                        {
                            m_results[e] = firstValue;
                        }
                        if (isSecondValueValid)
                        {
                            m_results[e] = secondValue;
                        }
                    }
                }
                if (m_locationTypes == Nodes)
                {
                    m_results = std::move(interpolatedResults);
                }
            }

            // copy values back
            return true;
        };

        bool SetResults(double** results);

    private:
        std::shared_ptr<Mesh> m_mesh;
        std::vector<Sample> m_samples;
        AveragingMethod m_averagingMethod;
        LocationTypes m_locationTypes;
        double m_relativeSearchRadius;

        SpatialTrees::RTree m_samplesRtree;
        std::vector<double> m_results;
    };
} // namespace meshkernel