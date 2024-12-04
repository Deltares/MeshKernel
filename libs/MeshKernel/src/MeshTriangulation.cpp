//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2024.
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

#include "MeshKernel/MeshTriangulation.hpp"
#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Utilities/Utilities.hpp"

#include <iostream>
#include <limits>

extern "C"
{
    /// @brief Function of the Triangle library
    ///
    /// \see https://www.cs.cmu.edu/~quake/triangle.html
    void Triangulation(int jatri,
                       double const* const xs,
                       double const* const ys,
                       int ns,
                       int* const indx,
                       int* const numtri,
                       int* const edgeidx,
                       int* const numedge,
                       int* const triedge,
                       double* const xs3,
                       double* const ys3,
                       int* const ns3,
                       double trisize);
}

void meshkernel::MeshTriangulation::Compute(const std::span<const double>& xNodes,
                                            const std::span<const double>& yNodes)
{

    if (xNodes.size() != yNodes.size())
    {
        throw ConstraintError("x- and y-points are not the same size: {} /= {}", xNodes.size(), yNodes.size());
    }

    if (xNodes.size() < 3)
    {
        throw ConstraintError("There are not enough points in mesh for a triangulation: {}", xNodes.size());
    }

    m_numEdges = 0;
    m_numFaces = 0;

    UInt estimatedNumberOfTriangles = 0;

    int numInputNodes = static_cast<int>(xNodes.size());

    if (estimatedNumberOfTriangles == 0)
    {
        estimatedNumberOfTriangles = static_cast<UInt>(xNodes.size()) * 6 + 10;
    }

    double averageTriangleArea = 0.0;

    int intTriangulationOption = 3;
    int numFaces = -1;
    int numEdges = 0;

    std::vector<int> faceNodes;
    std::vector<int> edgeNodes;
    std::vector<int> faceEdges;

    // If the number of estimated triangles is not sufficient, triangulation must be repeated
    while (numFaces < 0)
    {
        numFaces = static_cast<int>(estimatedNumberOfTriangles);

        faceNodes.resize(estimatedNumberOfTriangles * 3);
        std::ranges::fill(faceNodes, 0);

        edgeNodes.resize(estimatedNumberOfTriangles * 2);
        std::ranges::fill(edgeNodes, 0);

        faceEdges.resize(estimatedNumberOfTriangles * 3);
        std::ranges::fill(faceEdges, 0);

        Triangulation(intTriangulationOption,
                      xNodes.data(),
                      yNodes.data(),
                      numInputNodes,
                      faceNodes.data(), // INDX
                      &numFaces,
                      edgeNodes.data(), // EDGEINDX
                      &numEdges,
                      faceEdges.data(), // TRIEDGE
                      nullptr, nullptr, nullptr,
                      averageTriangleArea);

        if (estimatedNumberOfTriangles > 0)
        {
            estimatedNumberOfTriangles = static_cast<UInt>(-numFaces);
        }
    }

    m_numFaces = static_cast<UInt>(numFaces);
    m_numEdges = static_cast<UInt>(numEdges);

    m_faceNodes.resize(3 * m_numFaces);
    m_edgeNodes.resize(2 * m_numEdges);
    m_faceEdges.resize(3 * m_numFaces);

    std::transform(faceNodes.begin(), faceNodes.begin() + 3 * m_numFaces, m_faceNodes.begin(), [](const int val)
                   { return static_cast<UInt>(val - 1); });

    std::transform(faceEdges.begin(), faceEdges.begin() + 3 * m_numFaces, m_faceEdges.begin(), [](const int val)
                   { return static_cast<UInt>(val - 1); });

    std::transform(edgeNodes.begin(), edgeNodes.begin() + 2 * m_numEdges, m_edgeNodes.begin(), [](const int val)
                   { return static_cast<UInt>(val - 1); });

    m_elementCentres.resize(m_numFaces);

    for (UInt i = 0; i < m_numFaces; ++i)
    {
        m_elementCentres[i] = ComputeAverageCoordinate(GetNodes(i), m_projection);
    }

    if (m_numFaces > 0)
    {
        m_elementCentreRTree = RTreeFactory::Create(m_projection);
        m_elementCentreRTree->BuildTree(m_elementCentres);

        m_nodeRTree = RTreeFactory::Create(m_projection);
        m_nodeRTree->BuildTree(m_nodes);
    }

    m_edgesFaces.resize(m_numEdges, {constants::missing::uintValue, constants::missing::uintValue});
    size_t count = 0;

    for (UInt i = 0; i < m_numFaces; ++i)
    {

        for (UInt j = 0; j < constants::geometric::numNodesInTriangle; ++j)
        {
            UInt edge = m_faceEdges[count];
            ++count;

            if (m_edgesFaces[edge][0] == constants::missing::uintValue)
            {
                m_edgesFaces[edge][0] = i;
            }
            else
            {
                m_edgesFaces[edge][1] = i;
            }
        }
    }

    m_nodesEdges.resize(m_nodes.size());

    for (UInt e = 0; e < m_numEdges; ++e)
    {
        Edge edge = GetEdge(e);

        m_nodesEdges[edge.first].push_back(e);
        m_nodesEdges[edge.second].push_back(e);
    }
}

meshkernel::UInt meshkernel::MeshTriangulation::FindNearestFace(const Point& pnt) const
{
    m_elementCentreRTree->SearchNearestPoint(pnt);
    // return m_elementCentreRTree->HasQueryResults() ? m_elementCentreRTree->GetQueryResult(0) : constants::missing::uintValue;

    bool printIt = 728900.0 <= pnt.x && pnt.x <= 729100.0 && -5.60855e6 <= pnt.y && pnt.y <= -5.60845e6;

    if (m_elementCentreRTree->HasQueryResults())
    {
        UInt faceId = m_elementCentreRTree->GetQueryResult(0);

        if (printIt)
        {
            std::cout << "face id = " << faceId << std::endl;
        }

        if (PointIsInElement(pnt, faceId))
        {

            if (printIt)
            {
                auto nodes = GetNodes(faceId);

                std::cout << "face nodes: {"
                          << nodes[0].x << ", " << nodes[0].y << "} {"
                          << nodes[1].x << ", " << nodes[1].y << "} {"
                          << nodes[2].x << ", " << nodes[2].y << "} " << std::endl;

                auto nodeIds = GetNodeIds(faceId);

                std::cout << "face node ids: " << nodeIds[0] << ", " << nodeIds[1] << ", " << nodeIds[2] << std::endl;
            }

            return faceId;
        }
        else
        {
            const auto edgeIds = GetEdgeIds(faceId);

            BoundedArray<3 * MaximumNumberOfEdgesPerNode> elementsChecked;
            elementsChecked.push_back(faceId);

            // TODO collect the 3 neighbour elements ids, so that later
            // if we have to check the elements that are connected to a
            // node then we do not have to check those elements we have
            // checked already

            for (UInt i = 0; i < edgeIds.size(); ++i)
            {
                const auto [neighbour1, neighbour2] = GetFaceIds(edgeIds[i]);

                if (printIt)
                {
                    std::cout << "neighbours = " << neighbour1 << "  " << neighbour2 << std::endl;
                }

                if (neighbour1 != faceId && neighbour1 != constants::missing::uintValue)
                {
                    if (PointIsInElement(pnt, neighbour1))
                    {
                        return neighbour1;
                    }

                    elementsChecked.push_back(neighbour1);
                }

                if (neighbour2 != faceId && neighbour2 != constants::missing::uintValue)
                {
                    if (PointIsInElement(pnt, neighbour2))
                    {
                        return neighbour2;
                    }

                    elementsChecked.push_back(neighbour2);
                }
            }

            // Point not in direct neighbours of the element

            m_nodeRTree->SearchNearestPoint(pnt);

            if (m_nodeRTree->HasQueryResults())
            {
                UInt nodeId = m_nodeRTree->GetQueryResult(0);

                if (nodeId == constants::missing::uintValue)
                {
                    return constants::missing::uintValue;
                }

                for (UInt n = 0; n < m_nodesEdges[nodeId].size(); ++n)
                {
                    const auto faces = GetFaceIds(m_nodesEdges[nodeId][n]);

                    if (printIt)
                    {
                        std::cout << "faces = " << faces[0] << "  " << faces[1] << std::endl;
                    }

                    if (faces[0] != constants::missing::uintValue && !elementsChecked.contains(faces[0]))
                    {

                        if (PointIsInElement(pnt, faces[0]))
                        {
                            return faces[0];
                        }

                        elementsChecked.push_back(faces[0]);
                    }

                    if (faces[1] != constants::missing::uintValue && !elementsChecked.contains(faces[1]))
                    {

                        if (PointIsInElement(pnt, faces[1]))
                        {
                            return faces[1];
                        }

                        elementsChecked.push_back(faces[1]);
                    }
                }
            }
        }
    }

    // else
    {
        return constants::missing::uintValue;
    }

    // if (m_elementCentreRTree->HasQueryResults())
    // {

    //     for (UInt i = 0; i < m_elementCentreRTree->GetQueryResultSize(); ++i )
    //     {
    //         UInt id = m_elementCentreRTree->GetQueryResult(i);

    //         if (id < NumberOfFaces() && PointIsInElement (pnt, id))
    //         {
    //             return id;
    //         }

    //     }
    // }

    // // else
    // {
    //     return constants::missing::uintValue;
    // }
}

std::array<meshkernel::Point, 3> meshkernel::MeshTriangulation::GetNodes(const UInt faceId) const
{
    if (faceId == constants::missing::uintValue)
    {
        throw ConstraintError("Invalid face id");
    }

    if (faceId >= m_numFaces)
    {
        throw ConstraintError("face id out of range: {} >= {}", faceId, m_numFaces);
    }

    UInt faceIndexStart = 3 * faceId;

    return {GetNode(m_faceNodes[faceIndexStart]), GetNode(m_faceNodes[faceIndexStart + 1]), GetNode(m_faceNodes[faceIndexStart + 2])};
}

std::array<meshkernel::UInt, 3> meshkernel::MeshTriangulation::GetNodeIds(const UInt faceId) const
{
    if (faceId == constants::missing::uintValue)
    {
        throw ConstraintError("Invalid face id");
    }

    if (faceId >= m_numFaces)
    {
        throw ConstraintError("face id out of range: {} >= {}", faceId, m_numFaces);
    }

    UInt faceIndexStart = 3 * faceId;

    return {m_faceNodes[faceIndexStart], m_faceNodes[faceIndexStart + 1], m_faceNodes[faceIndexStart + 2]};
}

std::array<meshkernel::UInt, 3> meshkernel::MeshTriangulation::GetEdgeIds(const UInt faceId) const
{
    if (faceId == constants::missing::uintValue)
    {
        throw ConstraintError("Invalid face id");
    }

    if (faceId >= m_numFaces)
    {
        throw ConstraintError("face id out of range: {} >= {}", faceId, m_numFaces);
    }

    UInt faceIndexStart = 3 * faceId;

    return {m_faceEdges[faceIndexStart], m_faceEdges[faceIndexStart + 1], m_faceEdges[faceIndexStart + 2]};
}

bool meshkernel::MeshTriangulation::PointIsInElement(const Point& pnt, const UInt faceId) const
{
    if (faceId == constants::missing::uintValue)
    {
        throw ConstraintError("Invalid face id");
    }

    if (faceId >= m_numFaces)
    {
        throw ConstraintError("face id out of range: {} >= {}", faceId, m_numFaces);
    }

    if (!pnt.IsValid())
    {
        throw ConstraintError("Point is not valid");
    }

    return IsPointInTriangle(pnt, GetNodes(faceId), m_projection);
}

void meshkernel::MeshTriangulation::Print(std::ostream& out) const
{

    std::vector<Point> nodes(NumberOfNodes());
    std::vector<Edge> edges(NumberOfEdges());

    for (UInt i = 0; i < NumberOfNodes(); ++i)
    {
        nodes[i] = GetNode(i);
    }

    for (UInt i = 0; i < NumberOfEdges(); ++i)
    {
        edges[i] = GetEdge(i);
    }

    meshkernel::Print(nodes, edges, out);
}
