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

#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh1D.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Polygons.hpp>
#include <vector>

meshkernel::Mesh1D::Mesh1D(const std::vector<Edge>& edges,
                           const std::vector<Point>& nodes,
                           Projection projection) : Mesh(edges, nodes, projection){};

meshkernel::Mesh1D::Mesh1D(std::vector<std::vector<Point>> const& polyLines,
                           std::vector<std::vector<double>> const& fixedChainages,
                           double offset,
                           double minFaceSize,
                           double offsetFixedChainages,
                           Projection projection)
{
    if (offset <= minFaceSize)
    {
        throw std::invalid_argument("Mesh1D::Mesh1D: Offset cannot be smaller than merging length");
    }
    if (polyLines.size() != fixedChainages.size())
    {
        throw std::invalid_argument("Mesh1D::Mesh1D: The polyline vector and the fixed chainages vector size should be the same");
    }

    std::vector<Edge> edges;
    std::vector<Point> nodes;
    size_t numNodes = 0;
    for (auto i = 0; i < polyLines.size(); ++i)
    {
        auto const polyLineNodalChainages = ComputePolyLineChainages(polyLines[i], projection);
        std::vector<double> chainages;

        // Start and end are always discretized
        chainages.reserve(fixedChainages[i].size() * 2);
        chainages.push_back(polyLineNodalChainages.front());
        chainages.push_back(polyLineNodalChainages.back());

        // Add chainages at fixed nodes
        ComputeFixedChainages(fixedChainages[i],
                              minFaceSize,
                              offsetFixedChainages,
                              chainages);

        // Add chainages at fixed offsets
        ComputeChainagesAtInterval(offset, chainages);

        // Compute the discretization
        auto const computedDiscretization = RefinePolyLine(polyLines[i], chainages, projection);
        if (computedDiscretization.empty())
        {
            continue;
        }
        // add the new computed nodes
        std::copy(computedDiscretization.begin(), computedDiscretization.end(), back_inserter(nodes));

        // add the new computed edges
        for (auto i = numNodes; i < nodes.size() - 1; ++i)
        {
            edges.emplace_back(i, i + 1);
        }
        // Poly lines are separated. If the end of one polyline coincides with the start of another, the two nodes will be merged later on.
        numNodes = numNodes + nodes.size();
    }

    // Sets the edges, nodes and projections
    m_edges = edges;
    m_nodes = nodes;
    m_projection = projection;

    // Perform node administration to fill the internal arrays
    AdministrateNodesEdges();

    // If there are computational nodes at a distance smaller than  the threshold, these are eliminated
    const Polygons polygon{};
    MergeNodesInPolygon(polygon, minFaceSize);
}

void meshkernel::Mesh1D::ComputeFixedChainages(std::vector<double> const& fixedChainages,
                                               double minFaceSize,
                                               double offsetFromFixedChainages,
                                               std::vector<double>& chainages)
{
    if (fixedChainages.empty())
    {
        return;
    }

    if (chainages.size() < 2)
    {
        throw std::invalid_argument("ComputeFixedChainages: Start and end chainages not present");
    }

    double const startChainage = chainages[0];
    double const endChainage = chainages[1];

    double previousChainage = startChainage;
    bool previousChainageIsAFixedPoint = IsEqual(previousChainage, fixedChainages.front()) ? true : false;
    for (auto s = 0; s < fixedChainages.size(); ++s)
    {
        const auto chainageBeforeFixedPoint = fixedChainages[s] - offsetFromFixedChainages;

        if (chainageBeforeFixedPoint - previousChainage >= minFaceSize && chainageBeforeFixedPoint > startChainage)
        {
            chainages.emplace_back(chainageBeforeFixedPoint);
            previousChainage = chainageBeforeFixedPoint;
            previousChainageIsAFixedPoint = true;
        }
        else if (previousChainageIsAFixedPoint)
        {
            //center the gridpoint between two fixed points
            chainages.back() = (chainageBeforeFixedPoint + previousChainage) * 0.5;
            previousChainage = chainages.back();
        }

        const auto chainageAfterFixedPoint = fixedChainages[s] + offsetFromFixedChainages;
        if (chainageAfterFixedPoint - previousChainage >= minFaceSize && chainageAfterFixedPoint < endChainage)
        {
            chainages.emplace_back(chainageAfterFixedPoint);
            previousChainage = chainageAfterFixedPoint;
            previousChainageIsAFixedPoint = true;
        }
    }
}

void meshkernel::Mesh1D::ComputeChainagesAtInterval(double offset, std::vector<double>& chainages)
{
    // Sort whatever is there
    std::sort(chainages.begin(), chainages.end());
    std::vector<double> chainagesAtInterval;
    for (auto i = 1; i < chainages.size(); ++i)
    {
        double const segmentLength = chainages[i] - chainages[i - 1];
        if (segmentLength < offset)
        {
            continue;
        }
        auto const numberOfNewSegments = static_cast<size_t>(std::ceil(segmentLength / offset));
        for (auto j = 1; j < numberOfNewSegments; j++)
        {
            chainagesAtInterval.push_back(chainages[i - 1] + j * (segmentLength / numberOfNewSegments));
        }
    }

    // add the newly computed chainages at interval are added
    std::copy(chainagesAtInterval.begin(), chainagesAtInterval.end(), std::back_inserter(chainages));
}