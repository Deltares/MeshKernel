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
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/RTree.hpp>
#include <MeshKernel/TriangulationInterpolation.hpp>
#include <MeshKernel/TriangulationWrapper.hpp>

using meshkernel::TriangulationInterpolation;

TriangulationInterpolation::TriangulationInterpolation(const std::vector<Point>& m_locations,
                                                       const std::vector<Sample>& samples,
                                                       Projection projection) : m_locations(m_locations),
                                                                                m_samples(samples),
                                                                                m_projection(projection) {}

void TriangulationInterpolation::Compute()
{
    // allocate and initialize result vector
    m_results.resize(m_locations.size(), constants::missing::doubleValue);

    // no samples available, return
    if (m_samples.empty())
    {
        throw AlgorithmError("TriangulationInterpolation::Compute: No samples available.");
    }

    // triangulate samples
    TriangulationWrapper triangulationWrapper;
    triangulationWrapper.Compute(m_samples,
                                 TriangulationWrapper::TriangulationOptions::TriangulatePointsAndGenerateFaces,
                                 0.0,
                                 0);

    triangulationWrapper.BuildTriangulation();

    // no triangles formed, return
    if (triangulationWrapper.GetNumFaces() < 1)
    {
        throw AlgorithmError("TriangulationInterpolation::Compute: Triangulation of samples produced no triangles.");
    }

    // for each triangle compute the bounding circumcenter, bounding closed polygon, and the values at the nodes of each triangle
    std::vector<Point> trianglesCircumcenters(triangulationWrapper.GetNumFaces(), {constants::missing::doubleValue, constants::missing::doubleValue});
    std::vector<std::vector<Point>> triangles(triangulationWrapper.GetNumFaces(), std::vector<Point>(4));
    std::vector<std::vector<double>> values(triangulationWrapper.GetNumFaces(), std::vector<double>(4, constants::missing::doubleValue));

    for (auto f = 0; f < triangulationWrapper.GetNumFaces(); ++f)
    {
        // compute triangle polygons
        for (size_t n = 0; n < Mesh::m_numNodesInTriangle; ++n)
        {
            auto const node = triangulationWrapper.GetFaceNode(f, n);
            triangles[f][n] = {m_samples[node].x, m_samples[node].y};
            values[f][n] = m_samples[node].value;
        }
        triangles[f][3] = triangles[f][0];
        values[f][3] = values[f][0];

        trianglesCircumcenters[f] = ComputeAverageCoordinate(triangles[f], m_projection);
    }

    RTree samplesRtree;
    samplesRtree.BuildTree(trianglesCircumcenters);

    // compute the sample bounding box
    const auto [lowerLeft, upperRight] = GetBoundingBox(m_samples);

    // loop over locations
    for (size_t n = 0; n < m_locations.size(); ++n)
    {
        if (!IsValueInBoundingBox(m_locations[n], lowerLeft, upperRight) ||
            !IsEqual(m_results[n], constants::missing::doubleValue))
        {
            continue;
        }

        // find the nearest triangle
        samplesRtree.SearchNearestPoint(m_locations[n]);
        if (!samplesRtree.HasQueryResults())
        {
            continue;
        }
        auto triangle = samplesRtree.GetQueryResult(0);

        // search for the triangle where the location is included
        bool isInTriangle = false;
        int numFacesSearched = 0;
        while (!isInTriangle && numFacesSearched < 2 * triangulationWrapper.GetNumFaces() && triangle != constants::missing::sizetValue && static_cast<int>(triangle) < triangulationWrapper.GetNumFaces())
        {

            isInTriangle = IsPointInPolygonNodes(m_locations[n], triangles[triangle], m_projection, trianglesCircumcenters[triangle]);

            // valid triangle found, no need to search further
            if (isInTriangle)
            {
                break;
            }

            // proceed to next triangle, which is adjacent to the edge that is cut by the line from the current triangle to the point location
            numFacesSearched++;
            for (size_t i = 0; i < Mesh::m_numNodesInTriangle; ++i)
            {
                const auto edge = triangulationWrapper.GetFaceEdge(triangle, i);
                if (triangulationWrapper.GetEdgeFace(edge, 1) == 0)
                {
                    continue;
                }

                // there is no valid other triangle
                const auto otherTriangle = triangle == triangulationWrapper.GetEdgeFace(edge, 0) ? triangulationWrapper.GetEdgeFace(edge, 1) : triangulationWrapper.GetEdgeFace(edge, 0);
                const auto k1 = triangulationWrapper.GetEdgeNode(edge, 0);
                const auto k2 = triangulationWrapper.GetEdgeNode(edge, 1);
                Point intersection;
                double crossProduct;
                double firstRatio;
                double secondRatio;
                const auto areCrossing = AreSegmentsCrossing(trianglesCircumcenters[triangle],
                                                             m_locations[n],
                                                             {m_samples[k1].x, m_samples[k1].y},
                                                             {m_samples[k2].x, m_samples[k2].y},
                                                             false,
                                                             m_projection,
                                                             intersection,
                                                             crossProduct,
                                                             firstRatio,
                                                             secondRatio);

                if (areCrossing)
                {
                    triangle = otherTriangle;
                    break;
                }
            }
        }

        if (isInTriangle && triangle != constants::missing::sizetValue && static_cast<int>(triangle) < triangulationWrapper.GetNumFaces())
        {
            // Perform linear interpolation
            m_results[n] = LinearInterpolationInTriangle(m_locations[n], triangles[triangle], values[triangle], m_projection);
        }
    }
}
