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

#include <vector>

#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Network1D.hpp>
#include <MeshKernel/Operations.hpp>

meshkernel::Network1D::Network1D(std::vector<std::vector<Point>> const& polyLines,
                                 Projection projection) : m_polyLines(polyLines), m_projection(projection)
{
    m_chainages.resize(m_polyLines.size());

    // start and end polyline chainages should always be accounted for
    for (auto i = 0; i < m_polyLines.size(); ++i)
    {
        auto const nodalChainages = ComputePolyLineNodalChainages(m_polyLines[i], projection);
        m_chainages[i].push_back(nodalChainages.front());
        m_chainages[i].push_back(nodalChainages.back());
    }
}

void meshkernel::Network1D::ComputeFixedChainages(std::vector<std::vector<double>> const& fixedChainagesByPolyline,
                                                  double minFaceSize,
                                                  double fixedChainagesOffset)
{

    if (m_polyLines.size() != fixedChainagesByPolyline.size())
    {
        throw std::invalid_argument("Network1D::ComputeFixedChainages: The polyline vector and the fixed chainages vector size must be the same");
    }

    for (auto p = 0; p < m_polyLines.size(); ++p)
    {
        if (fixedChainagesByPolyline[p].empty())
        {
            continue;
        }
        // For the current polyline, the start and end chainages should always be present
        if (m_chainages[p].size() < 2)
        {
            continue;
        }

        double const startChainage = m_chainages[p][0];
        double const endChainage = m_chainages[p][1];

        double previousChainage = startChainage;
        bool previousChainageIsAFixedPoint = IsEqual(previousChainage, fixedChainagesByPolyline[p].front()) ? true : false;
        for (auto const& fixedChainage : fixedChainagesByPolyline[p])
        {
            if (const auto chainageBeforeFixedPoint = fixedChainage - fixedChainagesOffset; chainageBeforeFixedPoint - previousChainage >= minFaceSize && chainageBeforeFixedPoint > startChainage)
            {
                m_chainages[p].emplace_back(chainageBeforeFixedPoint);
                previousChainage = chainageBeforeFixedPoint;
                previousChainageIsAFixedPoint = true;
            }
            else if (previousChainageIsAFixedPoint)
            {
                //center the gridpoint between two fixed points
                m_chainages[p].back() = (chainageBeforeFixedPoint + previousChainage) * 0.5;
                previousChainage = m_chainages[p].back();
            }

            if (const auto chainageAfterFixedPoint = fixedChainage + fixedChainagesOffset; chainageAfterFixedPoint - previousChainage >= minFaceSize && chainageAfterFixedPoint < endChainage)
            {
                m_chainages[p].emplace_back(chainageAfterFixedPoint);
                previousChainage = chainageAfterFixedPoint;
                previousChainageIsAFixedPoint = true;
            }
        }
    }
}

void meshkernel::Network1D::ComputeOffsettedChainages(double offset)
{

    for (auto p = 0; p < m_polyLines.size(); ++p)
    {
        // Sort whatever is there
        std::sort(m_chainages[p].begin(), m_chainages[p].end());
        std::vector<double> chainagesAtInterval;
        for (auto i = 1; i < m_chainages[p].size(); ++i)
        {
            double const segmentLength = m_chainages[p][i] - m_chainages[p][i - 1];
            if (segmentLength < offset)
            {
                continue;
            }
            auto const numberOfNewSegments = static_cast<size_t>(std::ceil(segmentLength / offset));
            for (auto j = 1; j < numberOfNewSegments; j++)
            {
                chainagesAtInterval.push_back(m_chainages[p][i - 1] + j * (segmentLength / static_cast<double>(numberOfNewSegments)));
            }
        }

        // add the newly computed chainages at interval are added
        std::copy(chainagesAtInterval.begin(), chainagesAtInterval.end(), std::back_inserter(m_chainages[p]));
    }
}

std::vector<std::vector<meshkernel::Point>> meshkernel::Network1D::ComputeDiscretizationsFromChainages()
{
    std::vector<std::vector<Point>> result;
    for (auto p = 0; p < m_polyLines.size(); ++p)
    {
        result.emplace_back(ComputePolyLineDiscretization(m_polyLines[p], m_chainages[p], m_projection));
    }
    return result;
}