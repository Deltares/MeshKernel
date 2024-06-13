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

#include <algorithm>

#include "MeshKernel/CasulliDeRefinement.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/UndoActions/AddNodeAction.hpp"
#include "MeshKernel/UndoActions/CompoundUndoAction.hpp"
#include "MeshKernel/UndoActions/FullUnstructuredGridUndo.hpp"

std::unique_ptr<meshkernel::UndoAction> meshkernel::CasulliDeRefinement::Compute(Mesh2D& mesh) const
{
    Polygons emptyPolygon;
    return Compute(mesh, emptyPolygon);
}

std::unique_ptr<meshkernel::UndoAction> meshkernel::CasulliDeRefinement::Compute(Mesh2D& mesh, const Polygons& polygon) const
{
    std::unique_ptr<FullUnstructuredGridUndo> refinementAction = FullUnstructuredGridUndo::Create(mesh);
    const std::vector<Point> meshNodes(mesh.Nodes());
    const std::vector<Edge> meshEdges(mesh.Edges());

    if (DoDeRefinement(mesh, polygon))
    {
        mesh.DeleteInvalidNodesAndEdges();
        mesh.Administrate();
    }
    else
    {
        // The de-refinement failed, restore the original mesh
        refinementAction->Restore();
        // Reset the undo action to a null action
        refinementAction.reset(nullptr);
        throw AlgorithmError("Unable to compute the Casulli de-refinement");
    }

    return refinementAction;
}

void meshkernel::CasulliDeRefinement::FindDirectlyConnectedCells(const Mesh2D& mesh,
                                                                 const UInt elementId,
                                                                 std::vector<UInt>& connected) const
{
    connected.clear();

    for (UInt i = 0; i < mesh.m_numFacesNodes[elementId]; ++i)
    {
        UInt edgeId = mesh.m_facesEdges[elementId][i];

        if (mesh.m_edgesNumFaces[edgeId] < 2)
        {
            continue;
        }

        UInt neighbouringElementId = elementId == mesh.m_edgesFaces[edgeId][0] ? mesh.m_edgesFaces[edgeId][1] : mesh.m_edgesFaces[edgeId][0];
        connected.push_back(neighbouringElementId);
    }
}

void meshkernel::CasulliDeRefinement::FindIndirectlyConnectedCells(const Mesh2D& mesh,
                                                                   const UInt elementId,
                                                                   const std::vector<UInt>& directlyConnected,
                                                                   std::vector<UInt>& indirectlyConnected) const
{
    indirectlyConnected.clear();

    for (UInt i = 0; i < mesh.m_numFacesNodes[elementId]; ++i)
    {
        UInt nodeId = mesh.m_facesNodes[elementId][i];

        for (UInt j = 0; j < mesh.m_nodesNumEdges[nodeId]; ++j)
        {
            UInt edgeId = mesh.m_nodesEdges[nodeId][j];

            for (UInt k = 0; k < mesh.m_edgesNumFaces[edgeId]; ++k)
            {
                UInt otherElementId = mesh.m_edgesFaces[edgeId][k];

                if (elementId == otherElementId ||
                    std::ranges::find(directlyConnected, otherElementId) != directlyConnected.end() ||
                    std::ranges::find(indirectlyConnected, otherElementId) != indirectlyConnected.end())
                {
                    continue;
                }

                // Add new cell
                indirectlyConnected.push_back(otherElementId);
            }
        }
    }
}

void meshkernel::CasulliDeRefinement::FindAdjacentCells(const Mesh2D& mesh,
                                                        const std::vector<UInt>& directlyConnected,
                                                        const std::vector<UInt>& indirectlyConnected,
                                                        std::vector<std::array<int, 2>>& kne) const
{
    std::ranges::fill(kne, std::array<int, 2>{constants::missing::intValue, constants::missing::intValue});

    for (UInt i = 0; i < directlyConnected.size(); ++i)
    {
        UInt elementId = directlyConnected[i];

        for (UInt j = 0; j < mesh.m_numFacesNodes[elementId]; ++j)
        {
            UInt edgeId = mesh.m_facesEdges[elementId][j];

            if (mesh.m_edgesNumFaces[edgeId] < 2)
            {
                continue;
            }

            UInt neighbouringElementId = elementId == mesh.m_edgesFaces[edgeId][0] ? mesh.m_edgesFaces[edgeId][1] : mesh.m_edgesFaces[edgeId][0];

            for (UInt k = 0; k < directlyConnected.size(); ++k)
            {
                if (directlyConnected[k] == neighbouringElementId)
                {
                    if (kne[i][0] == constants::missing::intValue)
                    {
                        kne[i][0] = -static_cast<int>(neighbouringElementId);
                    }
                    else
                    {
                        kne[i][1] = -static_cast<int>(neighbouringElementId);
                    }

                    neighbouringElementId = constants::missing::uintValue;
                }
            }

            if (neighbouringElementId == constants::missing::uintValue)
            {
                continue;
            }

            for (UInt k = 0; k < indirectlyConnected.size(); ++k)
            {
                if (indirectlyConnected[k] == neighbouringElementId)
                {
                    if (kne[i][0] == constants::missing::intValue)
                    {
                        kne[i][0] = static_cast<int>(neighbouringElementId);
                    }
                    else
                    {
                        kne[i][1] = static_cast<int>(neighbouringElementId);
                    }
                }
            }
        }
    }
}

void meshkernel::CasulliDeRefinement::FindSurroundingCells(const Mesh2D& mesh,
                                                           const UInt elementId,
                                                           std::vector<UInt>& directlyConnected,
                                                           std::vector<UInt>& indirectlyConnected,
                                                           std::vector<std::array<int, 2>>& kne) const
{
    FindDirectlyConnectedCells(mesh, elementId, directlyConnected);
    FindIndirectlyConnectedCells(mesh, elementId, directlyConnected, indirectlyConnected);
    FindAdjacentCells(mesh, directlyConnected, indirectlyConnected, kne);
}

bool meshkernel::CasulliDeRefinement::ElementIsSeed(const Mesh2D& mesh,
                                                    const std::vector<int>& nodeTypes,
                                                    const UInt element) const
{
    bool isSeed = true;

    for (UInt i = 0; i < mesh.m_numFacesNodes[element]; ++i)
    {
        if (nodeTypes[mesh.m_facesNodes[element][i]] == 0)
        {
            isSeed = false;
            break;
        }
    }

    return isSeed;
}

meshkernel::UInt meshkernel::CasulliDeRefinement::FindElementSeedIndex(const Mesh2D& mesh,
                                                                       const std::vector<int>& nodeTypes) const
{
    UInt seedIndex = constants::missing::uintValue;

    for (UInt e = 0; e < mesh.Edges().size(); ++e)
    {
        if (mesh.m_edgesNumFaces[e] != 1)
        {
            continue;
        }

        if (const Edge& edge = mesh.GetEdge(e); nodeTypes[edge.first] != 2 || nodeTypes[edge.second] != 2)
        {
            continue;
        }

        UInt elementId = mesh.m_edgesFaces[e][0];

        if (mesh.m_numFacesNodes[elementId] != constants::geometric::numNodesInQuadrilateral)
        {
            continue;
        }

        if (ElementIsSeed(mesh, nodeTypes, elementId))
        {
            seedIndex = elementId;
            break;
        }
    }

    // No seed index found, select the first quadrilateral inside the selecting polygon
    if (seedIndex == constants::missing::uintValue)
    {
        for (UInt face = 0; face < mesh.GetNumFaces(); ++face)
        {

            if (mesh.m_numFacesNodes[face] != constants::geometric::numNodesInQuadrilateral)
            {
                continue;
            }

            if (!ElementIsSeed(mesh, nodeTypes, face))
            {
                continue;
            }

            seedIndex = face;
            break;
        }
    }

    // still no element found, so take the first
    if (seedIndex == constants::missing::uintValue)
    {
        seedIndex = 0;
    }

    return seedIndex;
}

void meshkernel::CasulliDeRefinement::AddElementToList(const Mesh& mesh, const std::vector<UInt>& frontList, std::vector<UInt>& frontListCopy, const UInt elementId) const
{
    if (elementId != constants::missing::uintValue)
    {

        if (mesh.m_numFacesNodes[elementId] != constants::geometric::numNodesInQuadrilateral)
        {
            return;
        }

        // If the id is not already on the frontList then add it to the copy.
        if (std::ranges::find(frontList, elementId) == frontList.end())
        {
            frontListCopy.push_back(elementId);
        }
    }
}

std::vector<meshkernel::CasulliDeRefinement::ElementType> meshkernel::CasulliDeRefinement::InitialiseElementType(const Mesh2D& mesh,
                                                                                                                 const std::vector<int>& nodeTypes) const
{
    UInt seedElement = FindElementSeedIndex(mesh, nodeTypes);
    UInt iterationCount = 0;

    std::vector<UInt> directlyConnected;
    std::vector<UInt> indirectlyConnected;
    std::vector<std::array<int, 2>> kne(maximumSize);
    std::vector<UInt> frontIndex;
    std::vector<UInt> frontIndexCopy;

    directlyConnected.reserve(maximumSize);
    indirectlyConnected.reserve(maximumSize);
    frontIndex.reserve(maximumSize);
    frontIndexCopy.reserve(maximumSize);

    std::vector<ElementType> cellMask(mesh.GetNumFaces(), ElementType::Unassigned);

    cellMask[seedElement] = ElementType::WasNodeFirst;
    frontIndex.push_back(seedElement);

    while (frontIndex.size() > 0 && iterationCount < maxIterationCount)
    {
        ++iterationCount;
        frontIndexCopy.clear();

        for (UInt i = 0; i < frontIndex.size(); ++i)
        {
            UInt elementId = frontIndex[i];

            // get the connected cells
            FindSurroundingCells(mesh, elementId, directlyConnected, indirectlyConnected, kne);

            if (cellMask[elementId] == ElementType::WasNodeFirst)
            {

                for (UInt j = 0; j < directlyConnected.size(); ++j)
                {
                    UInt connectedElementId = directlyConnected[j];

                    if (mesh.m_numFacesNodes[connectedElementId] != constants::geometric::numNodesInQuadrilateral)
                    {
                        continue;
                    }

                    if ((cellMask[connectedElementId] != ElementType::WasNodeFirst && cellMask[connectedElementId] != ElementType::WasNodeAfter) &&
                        (cellMask[connectedElementId] != ElementType::WasEdgeFirst && cellMask[connectedElementId] != ElementType::WasEdgeAfter))
                    {
                        cellMask[connectedElementId] = ElementType::WasEdgeFirst;
                        AddElementToList(mesh, frontIndex, frontIndexCopy, connectedElementId);
                    }
                }

                for (UInt j = 0; j < indirectlyConnected.size(); ++j)
                {
                    UInt connectedElementId = indirectlyConnected[j];

                    if (mesh.m_numFacesNodes[connectedElementId] != constants::geometric::numNodesInQuadrilateral)
                    {
                        continue;
                    }

                    if (cellMask[connectedElementId] != ElementType::WasCell)
                    {
                        cellMask[connectedElementId] = ElementType::WasCell;
                    }
                }

                cellMask[elementId] = ElementType::WasNodeAfter;
            }
            else if (cellMask[elementId] == ElementType::WasEdgeFirst)
            {

                for (UInt j = 0; j < directlyConnected.size(); ++j)
                {
                    UInt connectedElementId = directlyConnected[j];

                    if (mesh.m_numFacesNodes[connectedElementId] != constants::geometric::numNodesInQuadrilateral)
                    {
                        continue;
                    }

                    if ((cellMask[connectedElementId] != ElementType::WasCell) &&
                        (cellMask[connectedElementId] != ElementType::WasNodeFirst && cellMask[connectedElementId] != ElementType::WasNodeAfter) &&
                        (cellMask[connectedElementId] != ElementType::WasEdgeFirst && cellMask[connectedElementId] != ElementType::WasEdgeAfter))
                    {
                        cellMask[connectedElementId] = ElementType::WasNodeFirst;
                        AddElementToList(mesh, frontIndex, frontIndexCopy, connectedElementId);
                    }
                }

                for (UInt j = 0; j < indirectlyConnected.size(); ++j)
                {
                    UInt connectedElementId = indirectlyConnected[j];

                    if (mesh.m_numFacesNodes[connectedElementId] != constants::geometric::numNodesInQuadrilateral)
                    {
                        continue;
                    }

                    if ((cellMask[connectedElementId] != ElementType::WasEdgeFirst && cellMask[connectedElementId] != ElementType::WasEdgeAfter) &&
                        (cellMask[connectedElementId] != ElementType::WasNodeFirst && cellMask[connectedElementId] != ElementType::WasNodeAfter) &&
                        cellMask[connectedElementId] != ElementType::WasCell)
                    {
                        cellMask[connectedElementId] = ElementType::WasEdgeFirst;
                        AddElementToList(mesh, frontIndex, frontIndexCopy, connectedElementId);
                    }
                }

                cellMask[elementId] = ElementType::WasEdgeAfter;
            }
        }

        frontIndex = frontIndexCopy;
    }

    return cellMask;
}

bool meshkernel::CasulliDeRefinement::DoDeRefinement(Mesh2D& mesh, const Polygons& polygon) const
{
    std::vector<UInt> directlyConnected;
    std::vector<UInt> indirectlyConnected;
    std::vector<std::array<int, 2>> kne(maximumSize);
    std::vector<int> nodeTypes(ComputeNodeTypes(mesh, polygon));
    directlyConnected.reserve(maximumSize);
    indirectlyConnected.reserve(maximumSize);

    std::vector<ElementType> cellMask(InitialiseElementType(mesh, nodeTypes));
    mesh.ComputeCircumcentersMassCentersAndFaceAreas(true);

    for (UInt k = 0; k < cellMask.size(); ++k)
    {
        if (cellMask[k] == ElementType::WasNodeAfter && mesh.m_numFacesNodes[k] > 0)
        {
            FindSurroundingCells(mesh, k, directlyConnected, indirectlyConnected, kne);

            if (!DeleteElement(mesh, nodeTypes, polygon, k, directlyConnected, indirectlyConnected, kne))
            {
                return false;
            }
        }
    }

    return true;
}

bool meshkernel::CasulliDeRefinement::ElementCannotBeDeleted(const Mesh2D& mesh,
                                                             const std::vector<int>& nodeTypes,
                                                             const Polygons& polygon,
                                                             const UInt elementId) const
{
    bool noGo = false;

    // Firstly, check and see if
    //    the cell to be deleted is not a corner cell, but has a link that
    //       is internal and whose both nodes are marked as non-internal nodes
    // return if so

    for (UInt i = 0; i < mesh.m_numFacesNodes[elementId]; ++i)
    {
        UInt edgeId = mesh.m_facesEdges[elementId][i];

        if (mesh.GetEdge(edgeId).first == constants::missing::uintValue ||
            mesh.GetEdge(edgeId).second == constants::missing::uintValue)
        {
            noGo = true;
            break;
        }

        if (mesh.m_edgesNumFaces[edgeId] == 2 &&
            nodeTypes[mesh.GetEdge(edgeId).first] != 1 &&
            nodeTypes[mesh.GetEdge(edgeId).second] != 1)
        {
            noGo = true;
            break;
        }
    }

    //--------------------------------

    for (UInt i = 0; i < mesh.m_numFacesNodes[elementId]; ++i)
    {
        UInt nodeId = mesh.m_facesNodes[elementId][i];

        if (nodeTypes[nodeId] == 3 && mesh.m_nodesNumEdges[nodeId] <= 2)
        {
            noGo = false;
            break;
        }
    }

    //--------------------------------

    // check if all nodes are in the selecting polygon
    for (UInt i = 0; i < mesh.m_numFacesNodes[elementId]; ++i)
    {
        UInt nodeId = mesh.m_facesNodes[elementId][i];

        if (!polygon.IsPointInAnyPolygon(mesh.Node(nodeId)))
        {
            noGo = true;
            break;
        }
    }

    return noGo;
}

meshkernel::Point meshkernel::CasulliDeRefinement::ComputeNewNodeCoordinates(const Mesh2D& mesh,
                                                                             const std::vector<int>& nodeTypes,
                                                                             const UInt elementId) const
{

    Point newNode(0.0, 0.0);
    double factor = 0.0;

    for (UInt i = 0; i < mesh.m_numFacesNodes[elementId]; ++i)
    {
        double fac = 1.0;
        UInt nodeId = mesh.m_facesNodes[elementId][i];

        if (nodeTypes[nodeId] == 2 || nodeTypes[nodeId] == 4)
        {
            fac = 1.0e45;
        }
        else if (nodeTypes[nodeId] == 3)
        {
            factor = 1.0;
            newNode = mesh.Node(nodeId);
            break;
        }

        newNode += fac * mesh.Node(nodeId);
        factor += fac;
    }

    newNode /= factor;

    return newNode;
}

bool meshkernel::CasulliDeRefinement::UpdateDirectlyConnectedTriangleElements(Mesh2D& mesh,
                                                                              const UInt index,
                                                                              const UInt connectedElementId,
                                                                              const std::vector<std::array<int, 2>>& kne) const
{

    for (UInt j = 0; j < mesh.m_numFacesNodes[connectedElementId]; ++j)
    {
        UInt edgeId = mesh.m_facesEdges[connectedElementId][j];

        if (mesh.m_edgesNumFaces[edgeId] < 2 && !CleanUpEdge(mesh, edgeId))
        {
            return false;
        }
    }

    // Find adjacent direct neighbours
    UInt otherEdgeId = constants::missing::uintValue;

    for (UInt i = 0; i < 2; ++i)
    {
        UInt leftElementId = kne[index][i] == constants::missing::intValue || kne[index][i] < 0 ? constants::missing::uintValue : static_cast<UInt>(kne[index][i]);

        if (leftElementId == constants::missing::uintValue)
        {
            continue;
        }

        UInt oppositeSide = 1 - i;
        UInt rightElementId = kne[index][oppositeSide] == constants::missing::intValue || kne[index][oppositeSide] < 0 ? constants::missing::uintValue : static_cast<UInt>(kne[index][oppositeSide]);

        if (leftElementId == constants::missing::uintValue || rightElementId == constants::missing::uintValue)
        {
            continue;
        }

        UInt j;
        UInt edgeId = constants::missing::uintValue;

        // find the common link
        for (j = 0; j < mesh.m_numFacesNodes[leftElementId]; ++j)
        {
            edgeId = mesh.m_facesEdges[leftElementId][j];

            if (mesh.m_edgesNumFaces[edgeId] < 2)
            {
                continue;
            }

            if (mesh.m_edgesFaces[edgeId][0] == connectedElementId && mesh.m_edgesFaces[edgeId][1] == leftElementId)
            {
                if (rightElementId != constants::missing::uintValue)
                {
                    mesh.m_edgesFaces[edgeId][0] = rightElementId;
                }
                else
                {
                    mesh.m_edgesFaces[edgeId][0] = mesh.m_edgesFaces[edgeId][1];
                    mesh.m_edgesFaces[edgeId][1] = constants::missing::uintValue;
                    mesh.m_edgesNumFaces[edgeId] = 1;
                }

                break;
            }
            else if (mesh.m_edgesFaces[edgeId][1] == connectedElementId && mesh.m_edgesFaces[edgeId][0] == leftElementId)
            {
                if (rightElementId != constants::missing::uintValue)
                {
                    mesh.m_edgesFaces[edgeId][1] = rightElementId;
                }
                else
                {
                    mesh.m_edgesFaces[edgeId][1] = constants::missing::uintValue;
                    mesh.m_edgesNumFaces[edgeId] = 1;
                }

                break;
            }
        }

        if (otherEdgeId != constants::missing::uintValue)
        {
            mesh.m_facesEdges[leftElementId][j] = otherEdgeId;

            if (!CleanUpEdge(mesh, edgeId))
            {
                return false;
            }
        }

        otherEdgeId = edgeId;
    }

    // deactivate cell
    mesh.m_numFacesNodes[connectedElementId] = 0;
    return true;
}

void meshkernel::CasulliDeRefinement::UpdateDirectlyConnectedNonTriangleElements(Mesh2D& mesh,
                                                                                 const UInt index,
                                                                                 const UInt elementId,
                                                                                 const UInt connectedElementId) const
{
    // polygons of degree higher than three: remove node and link

    for (UInt j = 0; j < mesh.m_numFacesNodes[connectedElementId]; ++j)
    {
        UInt edgeId = mesh.m_facesEdges[connectedElementId][j];

        if (mesh.m_edgesNumFaces[edgeId] < 2)
        {
            continue;
        }

        if (mesh.m_edgesFaces[edgeId][0] == elementId || mesh.m_edgesFaces[edgeId][1] == elementId)
        {
            UInt ndum = mesh.m_numFacesNodes[connectedElementId] - 1;

            std::shift_left(mesh.m_facesEdges[connectedElementId].begin() + j, mesh.m_facesEdges[connectedElementId].begin() + ndum, 1);

            // remove one node per removed link
            // take the first node that has not been removed before, but not the node that is kept,
            // which is the first of the center cell

            UInt i = 0;

            while (i < mesh.m_numFacesNodes[connectedElementId] &&
                   mesh.m_facesNodes[connectedElementId][i] != mesh.GetEdge(edgeId).first &&
                   mesh.m_facesNodes[connectedElementId][i] != mesh.GetEdge(edgeId).second &&
                   mesh.m_facesNodes[connectedElementId][i] != mesh.m_facesNodes[elementId][0])
            {
                ++i;
            }

            if (index < mesh.m_numFacesNodes[connectedElementId])
            {
                std::shift_left(mesh.m_facesNodes[connectedElementId].begin() + i, mesh.m_facesNodes[connectedElementId].begin() + ndum, 1);
            }
            else
            {
                throw AlgorithmError("No node found in Casulli de-refinement");
            }

            mesh.m_numFacesNodes[connectedElementId] = ndum;
        }
    }
}

bool meshkernel::CasulliDeRefinement::UpdateDirectlyConnectedElements(Mesh2D& mesh,
                                                                      const UInt elementId,
                                                                      const std::vector<UInt>& directlyConnected,
                                                                      const std::vector<std::array<int, 2>>& kne) const
{
    for (UInt k = 0; k < directlyConnected.size(); ++k)
    {
        UInt connectedElementId = directlyConnected[k];

        if (mesh.m_numFacesNodes[connectedElementId] < constants::geometric::numNodesInQuadrilateral)
        {
            if (!UpdateDirectlyConnectedTriangleElements(mesh, k, connectedElementId, kne))
            {
                return false;
            }
        }
        else
        {
            UpdateDirectlyConnectedNonTriangleElements(mesh, k, elementId, connectedElementId);
        }
    }

    return true;
}

int meshkernel::CasulliDeRefinement::GetNodeCode(const Mesh2D& mesh,
                                                 const std::vector<int>& nodeTypes,
                                                 const UInt elementId) const
{
    int nodeCode = std::numeric_limits<int>::lowest();

    for (UInt i = 0; i < mesh.m_numFacesNodes[elementId]; ++i)
    {
        nodeCode = std::max(nodeCode, nodeTypes[mesh.m_facesNodes[elementId][i]]);
    }

    return nodeCode;
}

void meshkernel::CasulliDeRefinement::RedirectNodesOfConnectedElements(Mesh2D& mesh, const UInt elementId, const UInt nodeId, const std::vector<UInt>& indirectlyConnected) const
{
    for (UInt i = 0; i < indirectlyConnected.size(); ++i)
    {
        UInt connectedElementId = indirectlyConnected[i];

        if (mesh.m_numFacesNodes[connectedElementId] < 3)
        {
            mesh.m_numFacesNodes[connectedElementId] = 0;
        }

        for (UInt j = 0; j < mesh.m_numFacesNodes[connectedElementId]; ++j)
        {
            // Another node id of the element
            UInt faceNodeId = mesh.m_facesNodes[connectedElementId][j];

            for (UInt k = 1; k < mesh.m_numFacesNodes[elementId]; ++k)
            {
                if (faceNodeId == mesh.m_facesNodes[elementId][k])
                {
                    mesh.m_facesNodes[elementId][j] = nodeId;
                }
            }
        }
    }
}

bool meshkernel::CasulliDeRefinement::RemoveUnwantedBoundaryNodes(Mesh2D& mesh,
                                                                  const std::vector<int>& nodeTypes,
                                                                  const Polygons& polygon,
                                                                  const std::vector<UInt>& indirectlyConnected) const
{
    for (UInt kk = 0; kk < indirectlyConnected.size(); ++kk)
    {
        UInt connectedElementId = indirectlyConnected[kk];
        bool continueOuterLoop = false;

        if (mesh.m_numFacesNodes[connectedElementId] < 3)
        {
            mesh.m_numFacesNodes[connectedElementId] = 0;
        }

        for (UInt i = 0; i < mesh.m_numFacesNodes[connectedElementId]; ++i)
        {
            UInt edgeId = mesh.m_facesEdges[connectedElementId][i];

            if (mesh.m_edgesNumFaces[edgeId] == 1)
            {
                UInt im1 = NextCircularBackwardIndex(i, static_cast<UInt>(mesh.m_facesNodes[connectedElementId].size()));
                UInt previousEdgeId = mesh.m_facesEdges[connectedElementId][im1];

                if (mesh.m_edgesNumFaces[previousEdgeId] == 1)
                {
                    UInt nodeId = mesh.FindCommonNode(edgeId, previousEdgeId);

                    // [sic] weird
                    if (nodeId == constants::missing::uintValue)
                    {
                        continue;
                    }

                    // this node may be outside polygon: ignore
                    if (nodeTypes[nodeId] != 3 && polygon.IsPointInAnyPolygon(mesh.Node(nodeId)))
                    {
                        if (mesh.m_numFacesNodes[connectedElementId] > 3)
                        {
                            std::shift_left(mesh.m_facesNodes[connectedElementId].begin() + i, mesh.m_facesNodes[connectedElementId].end(), 1);
                            std::shift_left(mesh.m_facesEdges[connectedElementId].begin() + i, mesh.m_facesEdges[connectedElementId].end(), 1);

                            --mesh.m_numFacesNodes[connectedElementId];

                            // redirect node of the link that is kept
                            if (mesh.GetEdge(previousEdgeId).first == nodeId)
                            {
                                mesh.GetEdge(previousEdgeId).first = mesh.GetEdge(edgeId).first + mesh.GetEdge(edgeId).second - nodeId;
                            }
                            else
                            {
                                mesh.GetEdge(previousEdgeId).second = mesh.GetEdge(edgeId).first + mesh.GetEdge(edgeId).second - nodeId;
                            }

                            // delete other link

                            if (CleanUpEdge(mesh, edgeId))
                            {
                                return false;
                            }

                            mesh.SetNode(nodeId, {constants::missing::doubleValue, constants::missing::doubleValue});
                            continueOuterLoop = true;
                            break;
                        }
                        else
                        {
                            mesh.m_numFacesNodes[connectedElementId] = 0;

                            if (!CleanUpEdge(mesh, edgeId) || !CleanUpEdge(mesh, previousEdgeId))
                            {
                                return false;
                            }

                            // previous-previous edgeId (not a real word)
                            UInt antePreviousEdgeId = std::accumulate(mesh.m_facesEdges[connectedElementId].begin(), mesh.m_facesEdges[connectedElementId].begin() + 3, 0) - edgeId - previousEdgeId;

                            if (mesh.m_edgesNumFaces[antePreviousEdgeId] > 1)
                            {
                                if (mesh.m_edgesFaces[antePreviousEdgeId][0] == connectedElementId)
                                {
                                    mesh.m_edgesNumFaces[antePreviousEdgeId] = 1;
                                    mesh.m_edgesFaces[antePreviousEdgeId][0] = mesh.m_edgesFaces[antePreviousEdgeId][1];
                                }
                                else if (mesh.m_edgesFaces[antePreviousEdgeId][1] == connectedElementId)
                                {
                                    mesh.m_edgesNumFaces[antePreviousEdgeId] = 1;
                                }
                            }

                            continueOuterLoop = true;
                            break;
                        }
                    }
                }
            }
        }

        if (continueOuterLoop)
        {
            continue;
        }
    }

    return true;
}

bool meshkernel::CasulliDeRefinement::DeleteElement(Mesh2D& mesh,
                                                    std::vector<int>& nodeTypes,
                                                    const Polygons& polygon,
                                                    const UInt elementId,
                                                    const std::vector<UInt>& directlyConnected,
                                                    const std::vector<UInt>& indirectlyConnected,
                                                    const std::vector<std::array<int, 2>>& kne) const
{
    if (directlyConnected.empty() || indirectlyConnected.empty())
    {
        return true;
    }

    if (ElementCannotBeDeleted(mesh, nodeTypes, polygon, elementId))
    {
        return true;
    }

    Point newNode(ComputeNewNodeCoordinates(mesh, nodeTypes, elementId));

    for (UInt i = 0; i < mesh.m_numFacesNodes[elementId]; ++i)
    {
        mesh.SetNode(mesh.m_facesNodes[elementId][i], newNode);
    }

    UInt N = mesh.m_numFacesNodes[elementId];

    if (!UpdateDirectlyConnectedElements(mesh, elementId, directlyConnected, kne))
    {
        return false;
    }

    //--------------------------------
    // Set the node code

    UInt nodeId = mesh.m_facesNodes[elementId][0];
    nodeTypes[nodeId] = GetNodeCode(mesh, nodeTypes, elementId);

    // merge nodes
    mesh.SetNode(nodeId, newNode);

    for (UInt kk = 1; kk < N; ++kk)
    {
        [[maybe_unused]] auto undoAction = mesh.MergeTwoNodes(mesh.m_facesNodes[elementId][kk], nodeId);
    }

    // redirect nodes of indirectly connected cells, deactivate polygons of degree smaller than three
    RedirectNodesOfConnectedElements(mesh, elementId, nodeId, indirectlyConnected);

    // remove unwanted boundary node: a non-corner node that is shared by two boundary links
    if (!RemoveUnwantedBoundaryNodes(mesh, nodeTypes, polygon, indirectlyConnected))
    {
        return false;
    }
    // redirect nodes of directly connected cells and deactivate polygons of degree smaller than three
    RedirectNodesOfConnectedElements(mesh, elementId, nodeId, directlyConnected);

    // deactivate links
    for (UInt i = 0; i < mesh.m_numFacesNodes[elementId]; ++i)
    {
        UInt L = mesh.m_facesEdges[elementId][i];

        if (!CleanUpEdge(mesh, L))
        {
            return false;
        }
    }

    // deactivate cell
    mesh.m_numFacesNodes[elementId] = 0;
    return true;
}

bool meshkernel::CasulliDeRefinement::CleanUpEdge(Mesh2D& mesh, const UInt edgeId) const
{

    for (UInt i = 0; i < 2; ++i)
    {
        UInt nodeId = EdgeNodeIndex(mesh.GetEdge(edgeId), i);

        if (nodeId == constants::missing::uintValue)
        {
            // already cleaned up
            continue;
        }

        UInt startOffset = constants::missing::uintValue;

        for (UInt j = 0; j < mesh.m_nodesNumEdges[nodeId]; ++j)
        {
            if (mesh.m_nodesEdges[nodeId][j] == edgeId)
            {
                startOffset = j;
                break;
            }
        }

        if (startOffset != constants::missing::uintValue)
        {
            std::shift_left(mesh.m_nodesEdges[nodeId].begin() + startOffset, mesh.m_nodesEdges[nodeId].end(), 1);
        }
        else
        {
            // Error: no link found
            return false;
        }

        --mesh.m_nodesNumEdges[nodeId];
    }

    mesh.GetEdge(edgeId).first = constants::missing::uintValue;
    mesh.GetEdge(edgeId).second = constants::missing::uintValue;
    return true;
}

std::vector<int> meshkernel::CasulliDeRefinement::ComputeNodeTypes(const Mesh2D& mesh, const Polygons& polygon) const
{
    std::vector<int> nodeTypes(mesh.GetNumNodes(), 0);

    for (UInt i = 0; i < mesh.GetNumNodes(); ++i)
    {
        if (polygon.IsPointInAnyPolygon(mesh.Node(i)))
        {
            nodeTypes[i] = mesh.m_nodesTypes[i];
        }
    }

    return nodeTypes;
}

std::vector<meshkernel::Point> meshkernel::CasulliDeRefinement::ElementsToDelete(const Mesh2D& mesh, const Polygons& polygon) const
{
    std::vector<int> nodeTypes(ComputeNodeTypes(mesh, polygon));
    std::vector<ElementType> cellMask(InitialiseElementType(mesh, nodeTypes));
    std::vector<Point> elementCentres;
    elementCentres.reserve(cellMask.size());

    for (UInt k = 0; k < cellMask.size(); ++k)
    {
        if (cellMask[k] == ElementType::WasNodeAfter && mesh.m_numFacesNodes[k] > 0)
        {
            bool toDelete = false;

            for (UInt j = 0; j < mesh.m_numFacesNodes[k]; ++j)
            {
                if (nodeTypes[mesh.m_facesNodes[k][j]] > 0)
                {
                    toDelete = true;
                    break;
                }
            }

            if (toDelete)
            {
                elementCentres.push_back(mesh.m_facesMassCenters[k]);
            }
        }
    }

    return elementCentres;
}
