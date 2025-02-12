//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2025.
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

#include "MeshKernel/MeshFaceCenters.hpp"
#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Operations.hpp"

std::vector<meshkernel::Point> meshkernel::MeshFaceCenters::ComputeCircumcenters(const Mesh& mesh)
{
    std::vector<Point> faceCentres(mesh.GetNumFaces());
    ComputeCircumcenters(mesh, faceCentres);

    return faceCentres;
}

void meshkernel::MeshFaceCenters::ComputeCircumcenters(const Mesh& mesh, std::span<Point> faceCentres)
{
    if (faceCentres.size() != mesh.GetNumFaces())
    {
        throw ConstraintError("array for faceCentres values is not the correct size");
    }

    const auto numFaces = static_cast<int>(mesh.GetNumFaces());

    std::vector<UInt> numEdgeFacesCache;
    numEdgeFacesCache.reserve(constants::geometric::maximumNumberOfEdgesPerFace);
    std::vector<Point> polygonNodesCache;
    polygonNodesCache.reserve(constants::geometric::maximumNumberOfEdgesPerFace);

#pragma omp parallel for private(numEdgeFacesCache, polygonNodesCache)
    for (int f = 0; f < numFaces; f++)
    {

        UInt numberOfInteriorEdges = 0;
        const auto numberOfFaceNodes = mesh.GetNumFaceEdges(f);

        for (UInt n = 0; n < numberOfFaceNodes; ++n)
        {
            if (!mesh.IsEdgeOnBoundary(mesh.m_facesEdges[f][n]))
            {
                numberOfInteriorEdges += 1;
            }
        }

        if (numberOfInteriorEdges == 0)
        {
            faceCentres[f] = mesh.m_facesMassCenters[f];
        }
        else
        {
            // need to account for spherical coordinates. Build a polygon around a face
            mesh.ComputeFaceClosedPolygon(f, polygonNodesCache);
            numEdgeFacesCache.clear();

            for (UInt n = 0; n < numberOfFaceNodes; ++n)
            {
                numEdgeFacesCache.emplace_back(mesh.m_edgesNumFaces[mesh.m_facesEdges[f][n]]);
            }

            faceCentres[f] = ComputeFaceCircumenter(polygonNodesCache, numEdgeFacesCache, mesh.m_projection);
        }
    }
}

meshkernel::Point meshkernel::MeshFaceCenters::ComputeFaceCircumenter(std::vector<Point>& polygon,
                                                                      const std::vector<UInt>& edgesNumFaces,
                                                                      const Projection& projection)
{
    static constexpr double weightCircumCenter = 1.0; ///< Weight circum center

    std::array<Point, constants::geometric::maximumNumberOfNodesPerFace> middlePoints;
    std::array<Point, constants::geometric::maximumNumberOfNodesPerFace> normals;
    UInt pointCount = 0;

    const auto numNodes = static_cast<UInt>(polygon.size()) - 1;

    Point centerOfMass{0.0, 0.0};
    for (UInt n = 0; n < numNodes; ++n)
    {
        centerOfMass.x += polygon[n].x;
        centerOfMass.y += polygon[n].y;
    }

    centerOfMass /= static_cast<double>(numNodes);

    auto result = centerOfMass;
    if (numNodes == constants::geometric::numNodesInTriangle)
    {
        result = CircumcenterOfTriangle(polygon[0], polygon[1], polygon[2], projection);
    }
    else if (!edgesNumFaces.empty())
    {
        UInt numValidEdges = CountNumberOfValidEdges(edgesNumFaces, numNodes);

        if (numValidEdges > 1)
        {
            ComputeMidPointsAndNormals(polygon, edgesNumFaces, numNodes, middlePoints, normals, pointCount, projection);
            result = ComputeCircumCentre(centerOfMass, pointCount, middlePoints, normals, projection);
        }
    }

    for (UInt n = 0; n < numNodes; ++n)
    {
        polygon[n] = weightCircumCenter * polygon[n] + (1.0 - weightCircumCenter) * centerOfMass;
    }

    // The circumcenter is included in the face, then return the calculated circumcentre
    if (IsPointInPolygonNodes(result, polygon, projection))
    {
        return result;
    }

    // If the circumcenter is not included in the face,
    // the circumcenter will be placed at the intersection between an edge and the segment connecting the mass center with the circumcentre.
    for (UInt n = 0; n < numNodes; ++n)
    {
        const auto nextNode = NextCircularForwardIndex(n, numNodes);

        const auto [areLineCrossing,
                    intersection,
                    crossProduct,
                    firstRatio,
                    secondRatio] = AreSegmentsCrossing(centerOfMass, result, polygon[n], polygon[nextNode], false, projection);

        if (areLineCrossing)
        {
            result = intersection;
            break;
        }
    }

    return result;
}
