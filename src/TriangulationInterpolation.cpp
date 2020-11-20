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
#include <MeshKernel/TriangulationInterpolation.hpp>
#include <MeshKernel/TriangulationWrapper.hpp>

#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Operations.cpp>
#include <MeshKernel/SpatialTrees.hpp>
#include <MeshKernel/Exceptions.hpp>

meshkernel::TriangulationInterpolation::TriangulationInterpolation(const std::vector<Point>& m_locations,
                                                                   const std::vector<Sample>& samples,
                                                                   Projections projection) : m_locations(m_locations),
                                                                                             m_samples(samples),
                                                                                             m_projection(projection){};

void meshkernel::TriangulationInterpolation::Compute()
{
    // allocate and initialize result vector
    m_results.resize(m_locations.size(), doubleMissingValue);

    // no samples available, return
    if (m_samples.empty())
    {
        throw AlgorithmError("TriangulationInterpolation::Compute: No samples available.");
    }

    // triangulate samples
    const auto numPolygonNodes = static_cast<int>(m_samples.size());
    TriangulationWrapper triangulationWrapper;

    triangulationWrapper.Compute(m_samples,
                                 numPolygonNodes,
                                 TriangulationWrapper::TriangulationOptions::TriangulatePointsAndGenerateFaces,
                                 0.0,
                                 0);

    // no triangles formed, return
    if (triangulationWrapper.m_numFaces < 1)
    {
        throw AlgorithmError("TriangulationInterpolation::Compute: Triangulation of samples produced no triangles.");
    }

    // for each triangle compute the bounding circumcenter, bounding closed polygon, and the values at the nodes of each triangle
    std::vector<Point> trianglesCircumcenters(triangulationWrapper.m_numFaces, {doubleMissingValue, doubleMissingValue});
    std::vector<std::vector<Point>> triangles(triangulationWrapper.m_numFaces, std::vector<Point>(4, {doubleMissingValue, doubleMissingValue}));
    std::vector<std::vector<double>> values(triangulationWrapper.m_numFaces, std::vector<double>(4, doubleMissingValue));

    for (int f = 0; f < triangulationWrapper.m_numFaces; ++f)
    {
        // compute triangle polygons
        for (int n = 0; n < 3; ++n)
        {
            auto const node = triangulationWrapper.m_faceNodes[f][n];
            triangles[f][n] = {m_samples[node].x, m_samples[node].y};
            values[f][n] = m_samples[node].value;
        }
        triangles[f][3] = triangles[f][0];
        values[f][3] = values[f][0];

        trianglesCircumcenters[f] = ComputeAverageCoordinate(triangles[f], 3, m_projection);
    }

    SpatialTrees::RTree samplesRtree;
    samplesRtree.BuildTree(trianglesCircumcenters);

    // compute the sample bounding box
    Point lowerLeft;
    Point upperRight;
    GetBoundingBox(m_samples, lowerLeft, upperRight);

    // loop over locations
    for (int n = 0; n < m_locations.size(); ++n)
    {
        if (!IsValueInBoundingBox(m_locations[n], lowerLeft, upperRight) ||
            !IsEqual(m_results[n], doubleMissingValue))
        {
            continue;
        }

        // find the nearest triangle
        samplesRtree.NearestNeighbour(m_locations[n]);
        if (samplesRtree.GetQueryResultSize() <= 0)
        {
            continue;
        }
        auto triangle = samplesRtree.GetQuerySampleIndex(0);

        // search for the triangle where the location is included
        bool isInTriangle = false;
        int numFacesSearched = 0;
        while (!isInTriangle && numFacesSearched < 2 * triangulationWrapper.m_numFaces && triangle >= 0 && triangle < triangulationWrapper.m_numFaces)
        {

            isInTriangle = IsPointInPolygonNodes(m_locations[n], triangles[triangle], 0, 3, m_projection, trianglesCircumcenters[triangle]);

            // valid triangle found, no need to search further
            if (isInTriangle)
            {
                break;
            }

            // proceed to next triangle, which is adjacent to the edge that is cut by the line from the current triangle to the point location
            numFacesSearched++;
            for (int i = 0; i < 3; ++i)
            {
                const auto edge = triangulationWrapper.m_faceEdges[triangle][i];
                if (triangulationWrapper.m_edgesFaces[edge][1] == 0)
                {
                    continue;
                }

                // there is no valid other triangle
                const auto otherTriangle = triangulationWrapper.m_edgesFaces[edge][0] + triangulationWrapper.m_edgesFaces[edge][1] - triangle;
                const auto k1 = triangulationWrapper.m_edgeNodes[edge][0];
                const auto k2 = triangulationWrapper.m_edgeNodes[edge][1];
                Point intersection;
                double crossProduct;
                double firstRatio;
                double secondRatio;
                const auto areCrossing = AreLinesCrossing(trianglesCircumcenters[triangle],
                                                          m_locations[n],
                                                          {m_samples[k1].x, m_samples[k1].y},
                                                          {m_samples[k2].x, m_samples[k2].y},
                                                          false,
                                                          intersection,
                                                          crossProduct,
                                                          firstRatio,
                                                          secondRatio,
                                                          m_projection);

                if (areCrossing)
                {
                    triangle = otherTriangle;
                    break;
                }
            }
        }

        if (isInTriangle && triangle >= 0 && triangle < triangulationWrapper.m_numFaces)
        {
            // Perform linear interpolation
            m_results[n] = LinearInterpolationInTriangle(m_locations[n], triangles[triangle], values[triangle], m_projection);
        }
    }
}
