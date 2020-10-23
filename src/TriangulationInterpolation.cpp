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
#include "TriangulationInterpolation.hpp"
#include "TriangulationWrapper.cpp"

#include "Entities.hpp"
#include "Mesh.hpp"
#include "Operations.cpp"
#include "SpatialTrees.hpp"

meshkernel::TriangulationInterpolation::TriangulationInterpolation(const std::shared_ptr<Mesh>& mesh,
                                                                   const std::vector<Sample>& samples,
                                                                   InterpolationLocation locationType) : m_mesh(mesh),
                                                                                                         m_samples(samples),
                                                                                                         m_interpolationLocation(locationType){};

void meshkernel::TriangulationInterpolation::Compute()
{
    // first triangulate
    double averageTriangleArea = 0;
    int numPolygonNodes = int(m_samples.size()); // open polygon
    int numberOfTriangles = 0;

    TriangulationWrapper triangulationWrapper;

    triangulationWrapper.Compute(m_samples,
                                 numPolygonNodes,
                                 3,
                                 averageTriangleArea,
                                 numberOfTriangles);

    std::vector<Point> locations;
    int numLocations = 0;
    if (m_interpolationLocation == Nodes)
    {
        locations = m_mesh->m_nodes;
        numLocations = m_mesh->GetNumNodes();
    }
    if (m_interpolationLocation == Edges)
    {
        m_mesh->ComputeEdgesCenters();
        locations = m_mesh->m_edgesCenters;
        numLocations = m_mesh->GetNumEdges();
    }
    if (m_interpolationLocation == Faces)
    {
        locations = m_mesh->m_facesMassCenters;
        numLocations = m_mesh->GetNumFaces();
    }

    m_results.resize(numLocations, doubleMissingValue);

    // no triangles formed
    if (triangulationWrapper.m_numFaces < 1)
    {
        return;
    }

    // for each triangle compute the bounding circumcenter, bounding closed polygon, and the values
    std::vector<Point> trianglesCircumcenters(triangulationWrapper.m_numFaces, {doubleMissingValue, doubleMissingValue});
    std::vector<std::vector<Point>> triangles(triangulationWrapper.m_numFaces, std::vector<Point>(4, {doubleMissingValue, doubleMissingValue}));
    std::vector<std::vector<double>> values(triangulationWrapper.m_numFaces, std::vector<double>(4, doubleMissingValue));

    for (int f = 0; f < triangulationWrapper.m_numFaces; ++f)
    {
        double xCircumcenter = 0.0;
        double yCircumcenter = 0.0;
        for (int n = 0; n < 3; ++n)
        {
            auto const node = triangulationWrapper.m_faceNodes[f][n];
            xCircumcenter += m_samples[node].x;
            yCircumcenter += m_samples[node].y;
            triangles[f][n] = {m_samples[node].x, m_samples[node].y};
            values[f][n] = m_samples[node].value;
        }
        triangles[f][3] = triangles[f][0];
        values[f][3] = values[f][0];

        trianglesCircumcenters[f] = {xCircumcenter * oneThird,
                                     yCircumcenter * oneThird};
    }

    SpatialTrees::RTree samplesRtree;
    samplesRtree.BuildTree(trianglesCircumcenters);

    Point lowerLeft;
    Point upperRight;
    GetBoundingBox(m_samples, lowerLeft, upperRight);

    // Loop over nodes
    for (int n = 0; n < numLocations; ++n)
    {
        if (!IsValueInBoundingBox(locations[n], lowerLeft, upperRight) ||
            !IsDifferenceLessThanEpsilon(m_results[n], doubleMissingValue))
        {
            continue;
        }

        samplesRtree.NearestNeighbour(locations[n]);
        if (samplesRtree.GetQueryResultSize() <= 0)
        {
            continue;
        }
        auto triangle = samplesRtree.GetQuerySampleIndex(0);

        bool isInTriangle = false;
        int numFacesSearched = 0;
        while (!isInTriangle && numFacesSearched < 2 * triangulationWrapper.m_numFaces)
        {
            isInTriangle = IsPointInPolygonNodes(locations[n], triangles[triangle], 0, 3);
            if (isInTriangle)
            {
                break;
            }

            // proceed to next triangle, which is adjacent to the edge that is cut by the line from the current triangle to the query point
            numFacesSearched++;
            for (int i = 0; i < 3; ++i)
            {
                const auto edge = triangulationWrapper.m_faceEdges[triangle][i];
                if (triangulationWrapper.m_edgesFaces[edge][1] == 0)
                {
                    continue;
                }

                const auto otherTriangle = triangulationWrapper.m_edgesFaces[edge][0] + triangulationWrapper.m_edgesFaces[edge][1] - triangle;
                const auto k1 = triangulationWrapper.m_edgeNodes[edge][0];
                const auto k2 = triangulationWrapper.m_edgeNodes[edge][1];
                bool adimensional = false;
                Point intersection;
                double crossProduct;
                double firstRatio;
                double secondRatio;
                bool areCrossing = AreLinesCrossing(trianglesCircumcenters[triangle],
                                                    locations[n],
                                                    {m_samples[k1].x, m_samples[k1].y},
                                                    {m_samples[k2].x, m_samples[k2].y},
                                                    adimensional,
                                                    intersection,
                                                    crossProduct,
                                                    firstRatio,
                                                    secondRatio,
                                                    m_mesh->m_projection);

                if (areCrossing)
                {
                    triangle = otherTriangle;
                    break;
                }
            }
        }
        // Perform linear interpolation
        m_results[n] = LinearInterpolationInTriangle(locations[n], triangles[triangle], values[triangle], m_mesh->m_projection);
    }
}
