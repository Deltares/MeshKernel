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
#include <iostream> // REMOVE
using namespace std;

#include "MeshKernel/CasulliDeRefinement.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/UndoActions/AddNodeAction.hpp"
#include "MeshKernel/UndoActions/CompoundUndoAction.hpp"

std::unique_ptr<meshkernel::UndoAction> meshkernel::CasulliDeRefinement::Compute(Mesh2D& mesh)
{
    Polygons emptyPolygon;
    return Compute(mesh, emptyPolygon);
}

std::unique_ptr<meshkernel::UndoAction> meshkernel::CasulliDeRefinement::Compute(Mesh2D& mesh, const Polygons& polygon)
{
    std::vector<EdgeNodes> newNodes(mesh.GetNumEdges(), {constants::missing::uintValue, constants::missing::uintValue, constants::missing::uintValue, constants::missing::uintValue});
    std::vector<NodeMask> nodeMask(InitialiseNodeMask(mesh, polygon));
    std::unique_ptr<CompoundUndoAction> refinementAction = CompoundUndoAction::Create();

    [[maybe_unused]] UInt elementSeedIndex = FindElementSeedIndex(mesh, polygon);

    DoDeRefinement(mesh, polygon);

    // [[maybe_unsued]] const UInt numNodes = mesh.GetNumNodes();
    // [[maybe_unsued]] const UInt numEdges = mesh.GetNumEdges();
    // [[maybe_unsued]] const UInt numFaces = mesh.GetNumFaces();

    // refinementAction->Add(ComputeNewNodes(mesh, newNodes, nodeMask));
    // refinementAction->Add(ConnectNewNodes(mesh, newNodes, numNodes, numEdges, numFaces, nodeMask));
    // refinementAction->Add(Administrate(mesh, numNodes, nodeMask));
    return refinementAction;
}

void meshkernel::CasulliDeRefinement::FindSurroundingCells(Mesh2D& mesh,
                                                           const Polygons& polygon [[maybe_unused]],
                                                           const UInt kCell,
                                                           UInt& nDirect,
                                                           UInt& nIndirect,
                                                           std::vector<UInt>& kDirect,
                                                           std::vector<UInt>& kIndirect,
                                                           std::vector<std::array<int, 2>>& kne)
{

    // Find directly connected cells

    std::fill(kne.begin(), kne.end(), std::array<int, 2>{constants::missing::intValue, constants::missing::intValue});

    kDirect.clear();
    kIndirect.clear();

    nDirect = 0;
    nIndirect = 0;

    for (UInt kk = 0; kk < mesh.m_numFacesNodes[kCell]; ++kk)
    {
        UInt L = mesh.m_facesEdges[kCell][kk];

        if (mesh.m_edgesNumFaces[L] < 2)
        {
            continue;
        }

        UInt kCell2 = mesh.m_edgesFaces[L][0] + mesh.m_edgesFaces[L][1] - kCell; // The other cell?
        bool alreadyVisited = false;

        for (UInt kkk = 0; kkk < kDirect.size(); ++kkk)
        {
            if (kDirect[kkk] == kCell2)
            {
                alreadyVisited = true;
                break;
            }
        }

        if (alreadyVisited)
        {
            continue;
        }

        kDirect.push_back(kCell2);
    }

    // find the cells indirectly connected cells

    for (UInt kk = 0; kk < mesh.m_numFacesNodes[kCell]; ++kk)
    {
        UInt k1 = mesh.m_facesNodes[kCell][kk];

        for (UInt kkk = 0; kkk < mesh.m_nodesNumEdges[k1]; ++kkk)
        {
            UInt L = mesh.m_nodesEdges[k1][kkk];
            bool isFound = false;

            for (UInt i = 0; i < mesh.m_edgesNumFaces[L]; ++i)
            {
                UInt kCell2 = mesh.m_edgesFaces[L][i];

                if (kCell == kCell2)
                {
                    continue;
                }

                isFound = false;

                for (UInt kkkk = 0; kkkk < kDirect.size(); ++kkkk)
                {
                    if (kCell2 == kDirect[kkkk])
                    {
                        isFound = true;
                        break;
                    }
                }

                if (isFound)
                {
                    continue;
                }

                isFound = false;

                for (UInt kkkk = 0; kkkk < kIndirect.size(); ++kkkk)
                {
                    if (kCell2 == kIndirect[kkkk])
                    {
                        isFound = true;
                        break;
                    }
                }

                if (isFound)
                {
                    continue;
                }

                // Add new cell
                kIndirect.push_back(kCell2);
            }
        }
    }

    // Find the adjacent cells

    for (UInt i = 0; i < kDirect.size(); ++i)
    {
        UInt kcell1 = kDirect[i];

        for (UInt j = 0; j < mesh.m_numFacesNodes[kcell1]; ++j)
        {
            UInt L = mesh.m_facesEdges[kcell1][j];

            if (mesh.m_edgesNumFaces[L] < 2)
            {
                continue;
            }

            // std::cout << " kcell2 " << mesh.m_edgesFaces[L][0] + mesh.m_edgesFaces[L][1] << "   " << kcell1 << "  "
            //           << int(mesh.m_edgesFaces[L][0] + mesh.m_edgesFaces[L][1]) - int(kcell1) << std::endl;
            UInt kCell2 = mesh.m_edgesFaces[L][0] + mesh.m_edgesFaces[L][1] - kcell1;

            for (UInt kk = 0; kk < kDirect.size(); ++kk)
            {
                if (kDirect[kk] == kCell2)
                {
                    if (kne[i][0] == constants::missing::intValue)
                    {
                        kne[i][0] = -static_cast<int>(kCell2);
                    }
                    else
                    {
                        kne[i][1] = -static_cast<int>(kCell2);
                    }

                    kCell2 = constants::missing::uintValue;
                }
            }

            if (kCell2 == constants::missing::uintValue)
            {
                continue;
            }

            for (UInt kk = 0; kk < kIndirect.size(); ++kk)
            {
                if (kIndirect[kk] == kCell2)
                {
                    if (kne[i][0] == constants::missing::intValue)
                    {
                        kne[i][0] = static_cast<int>(kCell2);
                    }
                    else
                    {
                        kne[i][1] = static_cast<int>(kCell2);
                    }
                }
            }
        }
    }
}

bool meshkernel::CasulliDeRefinement::ElementIsSeed(Mesh2D& mesh, const Polygons& polygon [[maybe_unused]], const UInt face)
{
    bool isFace = true;

    for (UInt i = 0; i < mesh.m_numFacesNodes[face]; ++i)
    {
        if (mesh.m_nodesTypes[mesh.m_facesNodes[face][i]] == 0)
        {
            isFace = false;
            break;
        }
    }

    return isFace;
}

meshkernel::UInt meshkernel::CasulliDeRefinement::FindElementSeedIndex(Mesh2D& mesh, const Polygons& polygon)
{
    UInt seedIndex = constants::missing::uintValue;

    for (UInt e = 0; e < mesh.Edges().size(); ++e)
    {
        if (mesh.m_edgesNumFaces[e] != 1)
        {
            continue;
        }

        const Edge& edge = mesh.GetEdge(e);

        if (mesh.m_nodesTypes[edge.first] != 2 || mesh.m_nodesTypes[edge.second] != 2)
        {
            continue;
        }

        UInt k1 = mesh.m_edgesFaces[e][0];

        if (mesh.m_numFacesNodes[k1] != constants::geometric::numNodesInQuadrilateral)
        // if (mesh.m_facesNodes[k1].size() != constants::geometric::numNodesInQuadrilateral)
        {
            continue;
        }

        if (!ElementIsSeed(mesh, polygon, k1))
        {
            continue;
        }

        seedIndex = k1;
        break;
    }

    // No seed index found, select the first quadrilateral inside the selecting polygon
    if (seedIndex == constants::missing::uintValue)
    {
        for (UInt face = 0; face < mesh.GetNumFaces(); ++face)
        {

            if (mesh.m_numFacesNodes[face] != constants::geometric::numNodesInQuadrilateral)
            // if (mesh.m_facesNodes[face].size() != constants::geometric::numNodesInQuadrilateral)
            {
                continue;
            }

            if (!ElementIsSeed(mesh, polygon, face))
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

void meshkernel::CasulliDeRefinement::UpdateFrontList(const Mesh& mesh, const std::vector<UInt>& frontList, std::vector<UInt>& frontListCopy, const UInt kNew)
{

    if (kNew != constants::missing::uintValue)
    {

        if (mesh.m_numFacesNodes[kNew] != 4)
        {
            return;
        }

        for (UInt i = 0; i < frontList.size(); ++i)
        {
            if (frontList[i] == kNew)
            {
                return;
            }
        }

        frontListCopy.push_back(kNew);
    }
}

void meshkernel::CasulliDeRefinement::DoDeRefinement(Mesh2D& mesh, const Polygons& polygon)
{
    [[maybe_unused]] UInt seedElement = FindElementSeedIndex(mesh, polygon);
    [[maybe_unused]] UInt iterationCount = 0;
    [[maybe_unused]] UInt nMax = 10;              // fix
    [[maybe_unused]] UInt numFront = 1;           // fix
    [[maybe_unused]] UInt maxIterationCount = 10; // fix

    [[maybe_unused]] UInt maxNumFront = 10; // fix

    [[maybe_unused]] UInt nDirect = 0;                          // fix
    [[maybe_unused]] UInt nIndirect = 0;                        // fix
    [[maybe_unused]] std::vector<UInt> kDirect;                 // fix
    [[maybe_unused]] std::vector<UInt> kIndirect;               // fix
    [[maybe_unused]] std::vector<std::array<int, 2>> kne(nMax); // fix
    [[maybe_unused]] std::vector<UInt> frontIndex;              // fix
    [[maybe_unused]] std::vector<UInt> frontIndexCopy;          // fix

    std::vector<ElementMask> cellMask(mesh.GetNumFaces(), ElementMask::Unassigned);

    cellMask[seedElement] = ElementMask::A;
    frontIndex.push_back(seedElement);

    while (frontIndex.size() > 0 && iterationCount < maxIterationCount)
    {
        ++iterationCount;
        frontIndexCopy.clear();

        for (UInt i = 0; i < frontIndex.size(); ++i)
        {
            UInt k = frontIndex[i];
            UInt kOther = constants::missing::uintValue;

            FindSurroundingCells(mesh, polygon, k, nDirect, nIndirect, kDirect, kIndirect, kne);

            if (cellMask[k] == ElementMask::A)
            {

                for (UInt j = 0; j < kDirect.size(); ++j)
                {
                    kOther = kDirect[j];

                    if (mesh.m_numFacesNodes[kOther] != 4)
                    {
                        continue;
                    }

                    if ((cellMask[kOther] != ElementMask::A && cellMask[kOther] != ElementMask::NotA) && (cellMask[kOther] != ElementMask::B && cellMask[kOther] != ElementMask::NotB))
                    {
                        cellMask[kOther] = ElementMask::B;
                        UpdateFrontList(mesh, frontIndex, frontIndexCopy, kOther);
                    }
                }

                for (UInt j = 0; j < kIndirect.size(); ++j)
                {
                    kOther = kIndirect[j];

                    if (mesh.m_numFacesNodes[kOther] != 4)
                    {
                        continue;
                    }

                    if (cellMask[kOther] != ElementMask::C)
                    {
                        cellMask[kOther] = ElementMask::C;
                    }
                }

                cellMask[k] = ElementMask::NotA;
            }
            else if (cellMask[k] == ElementMask::B)
            {

                for (UInt j = 0; j < kDirect.size(); ++j)
                {
                    kOther = kDirect[j];

                    if (mesh.m_numFacesNodes[kOther] != 4)
                    {
                        continue;
                    }

                    if ((cellMask[kOther] != ElementMask::C) &&
                        (cellMask[kOther] != ElementMask::A && cellMask[kOther] != ElementMask::NotA) &&
                        (cellMask[kOther] != ElementMask::B && cellMask[kOther] != ElementMask::NotB))
                    {
                        cellMask[kOther] = ElementMask::A;
                        UpdateFrontList(mesh, frontIndex, frontIndexCopy, kOther);
                    }
                }

                for (UInt j = 0; j < kIndirect.size(); ++j)
                {
                    kOther = kIndirect[j];

                    if (mesh.m_numFacesNodes[kOther] != 4)
                    {
                        continue;
                    }

                    if ((cellMask[kOther] != ElementMask::B && cellMask[kOther] != ElementMask::NotB) &&
                        (cellMask[kOther] != ElementMask::A && cellMask[kOther] != ElementMask::NotA) &&
                        cellMask[kOther] != ElementMask::C)
                    {
                        cellMask[kOther] = ElementMask::B;
                        UpdateFrontList(mesh, frontIndex, frontIndexCopy, kOther);
                    }
                }

                cellMask[k] = ElementMask::NotB;
            }
        }

        frontIndex = frontIndexCopy;
    }

    //--------------------------------

    for (UInt k = 0; k < cellMask.size(); ++k)
    {
        if (cellMask[k] == ElementMask::NotA && mesh.m_numFacesNodes[k] > 0)
        {
            FindSurroundingCells(mesh, polygon, k, nDirect, nIndirect, kDirect, kIndirect, kne);

            [[maybe_unused]] UInt k1 = mesh.m_facesNodes[k][0];
            bool deletedCell;

            DeleteCell(mesh, polygon, k, kDirect, kIndirect, kne, deletedCell);
        }
    }
}

void meshkernel::CasulliDeRefinement::DeleteCell(Mesh2D& mesh,
                                                 const Polygons& polygon,
                                                 const UInt k,
                                                 const std::vector<UInt>& kDirect,
                                                 const std::vector<UInt>& kIndirect,
                                                 const std::vector<std::array<int, 2>>& kne [[maybe_unused]],
                                                 bool& jaDeleted)
{
    jaDeleted = false;

    if (kDirect.size() == 0 || kIndirect.size() == 0)
    {
        return;
    }

    bool noGo = false;
    UInt kk = 0;

    // Firstly, check and see if
    //    the cell to be deleted is not a corner cell, but has a link that
    //       is internal and whose both nodes are marked as non-internal nodes
    // return if so
    while (!noGo && kk < mesh.m_numFacesNodes[k])
    {
        ++kk;
        UInt klin = mesh.m_facesEdges[k][kk];

        if (mesh.GetEdge(klin).first == constants::missing::uintValue ||
            mesh.GetEdge(klin).second == constants::missing::uintValue)
        {
            noGo = true;
            break;
        }

        if (mesh.m_edgesNumFaces[klin] == 2)
        {
            if (mesh.m_nodesTypes[mesh.GetEdge(klin).first] != 1 &&
                mesh.m_nodesTypes[mesh.GetEdge(klin).second] != 1)
            {
                noGo = true;
            }
        }
    }

    //--------------------------------

    //  check if the cell is a cornercell
    kk = 0;

    // TOOD Change to for loop and break
    while ((kk < mesh.m_numFacesNodes[k] - 1) && noGo)
    // while ((kk < mesh.m_facesNodes[k].size() - 1) && noGo)
    {
        ++kk;
        UInt knod = mesh.m_facesNodes[k][kk];

        if (mesh.m_nodesTypes[knod] == 3 && mesh.m_nodesNumEdges[knod] < 2)
        {
            noGo = false;
        }
    }

    //--------------------------------

    //  check if all nodes are in the selecting polygon

    for (UInt kk = 0; kk < mesh.m_numFacesNodes[k]; ++kk)
    {
        UInt k1 = mesh.m_facesNodes[k][kk];

        if (!polygon.IsPointInAnyPolygon(mesh.Node(k1)))
        {
            noGo = true;
            break;
        }
    }

    //--------------------------------

    if (noGo)
    {
        jaDeleted = false;
        return;
    }

    //--------------------------------

    std::vector<UInt> savedNodes;
    std::vector<Edge> savedEdges;
    std::vector<UInt> saveCells;

    StoreMesh(mesh, k, kDirect, kIndirect, savedNodes, savedEdges, saveCells);

    //--------------------------------

    UInt N = mesh.m_numFacesNodes[k];
    Point p(0.0, 0.0);

    double factor = 0.0;

    for (UInt kk = 0; kk < N; ++kk)
    {
        double fac = 1.0;
        UInt k1 = mesh.m_facesNodes[k][kk];

        if (mesh.m_nodesTypes[k1] == 2 || mesh.m_nodesTypes[k1] == 4)
        {
            fac = 1.0e45;
        }
        else if (mesh.m_nodesTypes[k1] == 3)
        {
            factor = 1.0;
            p = mesh.Node(k1);
            break;
        }

        p += fac * mesh.Node(k1);
        factor += fac;
    }

    p /= factor;

    // std::transform(mesh.m_facesNodes[k].begin(), mesh.m_facesNodes[k].end(),
    //                [&grid = mesh, newP = p](const UInt id)
    //                { mesh.Node(id) = newP; });

    for (UInt i = 0; i < mesh.m_numFacesNodes[k]; ++i)
    // for (UInt i = 0; i < mesh.m_facesNodes[k].size(); ++i)
    {
        std::cout << "reset node: from " << mesh.m_facesNodes[k][i] << "  " << mesh.Node(mesh.m_facesNodes[k][i]).x << ", " << mesh.Node(mesh.m_facesNodes[k][i]).y << "    to " << p.x << ", " << p.y << std::endl;
        auto undoAction = mesh.ResetNode(mesh.m_facesNodes[k][i], p);
    }

    for (UInt kk = 0; kk < kDirect.size(); ++kk)
    {
        UInt kcell1 = kDirect[kk];

        if (mesh.m_numFacesNodes[kcell1] < 4)
        {
            for (UInt j = 0; j < mesh.m_numFacesNodes[kcell1]; ++j)
            {
                UInt L = mesh.m_facesEdges[kcell1][j];

                if (mesh.m_edgesNumFaces[L] < 2)
                {
                    CleanUpLink(mesh, L);
                }
            }
        }
    }

    //--------------------------------

    // alter directly connected cells

    for (UInt kk = 0; kk < kDirect.size(); ++kk)
    {
        UInt kcell1 = kDirect[kk];

        if (mesh.m_numFacesNodes[kcell1] < 4)
        {

            for (UInt j = 0; j < mesh.m_numFacesNodes[kcell1]; ++j)
            {
                UInt L = mesh.m_facesEdges[kcell1][j];

                if (mesh.m_edgesNumFaces[L] < 2)
                {
                    CleanUpLink(mesh, L);
                }
            }

            // Find adjacent direct neighbours
            UInt L1 = 0;

            for (UInt i = 0; i < 2; ++i)
            {
                UInt kcL = kne[kk][i] == constants::missing::intValue ? constants::missing::uintValue : kne[kk][i];

                if (kcL == constants::missing::uintValue)
                {
                    continue;
                }
                // if (kcL == constants::missing::uintValue)
                // {
                //     continue;
                // }

                UInt iR = i == 1 ? 0 : 1; // 1 - i

                UInt kcR = kne[kk][iR] == constants::missing::intValue ? constants::missing::uintValue : kne[kk][iR];
                // UInt kcR = kne[kk][iR];

                if (kcL == constants::missing::uintValue || kcR == constants::missing::uintValue)
                {
                    continue;
                }

                UInt j;
                UInt L = constants::missing::uintValue;

                // find the common link
                for (j = 0; j < mesh.m_numFacesNodes[kcL]; ++j)
                {
                    L = mesh.m_facesEdges[kcL][j];

                    if (mesh.m_edgesNumFaces[L] < 2)
                    {
                        continue;
                    }

                    if (mesh.m_edgesFaces[L][0] == kcell1 && mesh.m_edgesFaces[L][1] == kcL)
                    {
                        if (kcR != 0)
                        {
                            mesh.m_edgesFaces[L][0] = kcR;
                        }
                        else
                        {
                            mesh.m_edgesFaces[L][0] = mesh.m_edgesFaces[L][1];
                            mesh.m_edgesFaces[L][1] = constants::missing::uintValue;
                            mesh.m_edgesNumFaces[L] = 1;
                        }

                        break;
                    }
                    else if (mesh.m_edgesFaces[L][1] == kcell1 && mesh.m_edgesFaces[L][0] == kcL)
                    {
                        if (kcR != 0)
                        {
                            mesh.m_edgesFaces[L][1] = kcR;
                        }
                        else
                        {
                            mesh.m_edgesFaces[L][1] = constants::missing::uintValue;
                            mesh.m_edgesNumFaces[L] = 1;
                        }
                    }

                    break;
                }
                if (L1 != constants::missing::uintValue)
                {
                    mesh.m_facesEdges[kcL][j] = L1;
                    CleanUpLink(mesh, L);
                }

                L1 = L;
            }

            // deactivate cell
            mesh.m_numFacesNodes[kcell1] = 0;
        }
        else
        {
            for (UInt j = 0; j < mesh.m_numFacesNodes[kcell1]; ++j)
            {
                UInt L = mesh.m_facesEdges[kcell1][j];

                if (mesh.m_edgesNumFaces[L] < 2)
                {
                    continue;
                }

                if (mesh.m_edgesFaces[L][0] == k || mesh.m_edgesFaces[L][1] == k)
                {
                    UInt ndum = mesh.m_numFacesNodes[kcell1] - 1;

                    // TODO use ShiftArray(i, mesh.m_XXX);
                    // ShiftArray(j, mesh.m_facesEdges[kcell1]);

                    for (UInt jj = j; jj < ndum; ++jj)
                    {
                        mesh.m_facesEdges[kcell1][jj] = mesh.m_facesEdges[kcell1][jj + 1];
                    }

                    // remove one node per removed link
                    // take the first node that has not been removed before, but not the node that is kept,
                    // which is the first of the center cell

                    UInt i = 0;

                    while (i < mesh.m_numFacesNodes[kcell1] &&
                           mesh.m_facesNodes[kcell1][i] != mesh.GetEdge(L).first &&
                           mesh.m_facesNodes[kcell1][i] != mesh.GetEdge(L).second &&
                           mesh.m_facesNodes[kcell1][i] != mesh.m_facesNodes[k][0])
                    {
                        ++i;
                    }

                    if (kk < mesh.m_numFacesNodes[kcell1])
                    {
                        // TODO use ShiftArray(i, mesh.m_XXX);
                        // ShiftArray(j, mesh.m_facesNodes[kcell1]);

                        for (UInt jj = i; jj < ndum; ++jj)
                        {
                            mesh.m_facesNodes[kcell1][jj] = mesh.m_facesNodes[kcell1][jj + 1];
                        }
                    }
                    else
                    {
                        // No node found
                    }

                    mesh.m_numFacesNodes[kcell1] = ndum;
                }
            }
        }
    }

    //--------------------------------
    // Set the node code

    UInt k1 = mesh.m_facesNodes[k][0];
    int maxVal = std::numeric_limits<int>::lowest();

    for (UInt i = 0; i < mesh.m_numFacesNodes[k]; ++i)
    // for (UInt i = 0; i < mesh.m_facesNodes[k].size(); ++i)
    {
        maxVal = std::max(maxVal, mesh.m_nodesTypes[mesh.m_facesNodes[k][i]]);
    }

    // TODO check maxVal != std::numeric_limits<int>::lowest()
    // TODO or same as mesh.m_facesNodes[k].size() != 0

    mesh.m_nodesTypes[k1] = maxVal;

    // merge nodes

    std::cout << "reset node: from " << k1 << "  " << mesh.Node(k1).x << ", " << mesh.Node(k1).y << "    to " << p.x << ", " << p.y << std::endl;
    auto undoAction = mesh.ResetNode(k1, p);
    // mesh.Node(k1) = p;

    for (UInt kk = 1; kk < N; ++kk)
    {
        [[maybe_unused]] auto undoAction = mesh.MergeTwoNodes(mesh.m_facesNodes[k][kk], k1);
    }

    // redirect nodes of indirectly connected cells, deactivate polygons of degree smaller than three and
    // remove unwanted boundary node

    for (UInt kk = 0; kk < kIndirect.size(); ++kk)
    {
        UInt kcell = kIndirect[kk];

        if (mesh.m_numFacesNodes[kcell] < 3)
        {
            mesh.m_numFacesNodes[kcell] = 0;
        }

        for (UInt i = 0; i < mesh.m_numFacesNodes[kcell]; ++i)
        {
            UInt k2 = mesh.m_facesNodes[kcell][i];
            // UInt L = mesh.m_facesEdges[kcell][i];

            for (UInt j = 1; j < mesh.m_numFacesNodes[k]; ++j)
            {
                if (k2 == mesh.m_facesNodes[k][j])
                {
                    mesh.m_facesNodes[kcell][i] = k1;
                }
            }
        }
    }

    // remove unwanted boundary node: a non-corner node that is shared by two boundary links

    for (UInt kk = 0; kk < kIndirect.size(); ++kk)
    {
        UInt kcell = kIndirect[kk];
        bool continueOuterLoop = false;

        if (mesh.m_numFacesNodes[kcell] < 3)
        {
            mesh.m_numFacesNodes[kcell] = 0;
        }

        for (UInt i = 0; i < mesh.m_numFacesNodes[kcell]; ++i)
        // for (UInt i = 0; i < mesh.m_facesNodes[kcell].size(); ++i)
        {
            UInt K2 = mesh.m_facesNodes[kcell][i];
            UInt L = mesh.m_facesEdges[kcell][i];

            if (mesh.m_edgesNumFaces[L] == 1)
            {
                UInt im1 = RotateIndex(i, mesh.m_facesNodes[kcell], false /*forward*/);
                UInt L1 = mesh.m_facesEdges[kcell][im1];

                if (mesh.m_edgesNumFaces[L1] == 1)
                {
                    K2 = FindCommonNode(mesh, L, L1);

                    // weird
                    if (K2 == constants::missing::uintValue)
                    {
                        continue;
                    }

                    // this node may be outside polygon: ignore
                    if (mesh.m_nodesTypes[K2] != 3 && polygon.IsPointInAnyPolygon(mesh.Node(K2)))
                    {
                        ShiftArray(i, mesh.m_facesNodes[kcell]);
                        ShiftArray(i, mesh.m_facesEdges[kcell]);
                        --mesh.m_numFacesNodes[kcell];
                        // redirect node of the link that is kept

                        if (mesh.GetEdge(L1).first == K2)
                        {
                            mesh.GetEdge(L1).first = mesh.GetEdge(L).first + mesh.GetEdge(L).second - K2;
                        }
                        else
                        {
                            mesh.GetEdge(L1).second = mesh.GetEdge(L).first + mesh.GetEdge(L).second - K2;
                        }

                        // delete other link

                        CleanUpLink(mesh, L);
                        std::cout << "reset node: from " << K2 << "  " << mesh.Node(K2).x << ", " << mesh.Node(K2).y << "    to " << constants::missing::doubleValue << ", " << constants::missing::doubleValue << std::endl;
                        auto undoAcction = mesh.ResetNode(K2, {constants::missing::doubleValue, constants::missing::doubleValue});
                        continueOuterLoop = true;
                        break;
                    }
                    else
                    {
                        mesh.m_numFacesNodes[kcell] = 0;
                        CleanUpLink(mesh, L);
                        CleanUpLink(mesh, L1);
                        UInt L2 = std::accumulate(mesh.m_facesEdges[kcell].begin(), mesh.m_facesEdges[kcell].begin() + 3, 0) - L - L1;

                        if (mesh.m_edgesNumFaces[L2] > 1)
                        {
                            if (mesh.m_edgesFaces[L2][0] == kcell)
                            {
                                mesh.m_edgesNumFaces[L2] = 1;
                                mesh.m_edgesFaces[L2][0] = mesh.m_edgesFaces[L2][1];
                            }
                            else if (mesh.m_edgesFaces[L2][1] == kcell)
                            {
                                mesh.m_edgesNumFaces[L2] = 1;
                            }
                        }

                        continueOuterLoop = true;
                        break;
                    }
                }
            }
        }

        if (continueOuterLoop)
        {
            continue;
        }
    }

    // redirect nodes of directly connected cells and deactivate polygons of degree smaller than three

    for (UInt kk = 0; kk < kDirect.size(); ++kk)
    {
        UInt kcell = kDirect[kk];

        if (mesh.m_numFacesNodes[kcell] < 3)
        {
            mesh.m_numFacesNodes[kcell] = 0;
        }

        for (UInt i = 0; i < mesh.m_numFacesNodes[kcell]; ++i)
        {
            UInt k2 = mesh.m_facesNodes[kcell][i];
            // UInt L = mesh.m_facesEdges[kcell][i];

            for (UInt j = 1; j < mesh.m_numFacesNodes[kcell]; ++j)
            {
                if (k2 == mesh.m_facesNodes[kcell][j])
                {
                    mesh.m_facesNodes[kcell][i] = k1;
                }
            }
        }
    }

    // deactivate links
    for (UInt kk = 0; kk < mesh.m_numFacesNodes[k]; ++kk)
    {
        UInt L = mesh.m_facesEdges[k][kk];
        CleanUpLink(mesh, L);
    }

    // deactivate cell
    mesh.m_numFacesNodes[k] = 0;

    // jadeleted = 1;

    // 1234 continue

    if (!noGo)
    {
        //    check and see if
        //        the altered indirectly connected cells are non convex
        //     restore and return if so

        //    check convexity
        noGo = false;
        UInt kk = 0;

        // TODO change to for loop
        while (kk < kIndirect.size())
        {
            if (IsConvexCell(mesh, kIndirect[kk]))
            {
                noGo = true;
                break;
            }

            ++kk;
        }

        // if (noGo)
        // {
        //     ja = false;
        // }

        // if (ja)
        // {
        //     // local_netrestore ();
        // }
    }
}

meshkernel::UInt meshkernel::CasulliDeRefinement::FindCommonNode(const Mesh2D& mesh, const UInt L1, const UInt L2)
{
    const Edge& e1 = mesh.GetEdge(L1);
    const Edge& e2 = mesh.GetEdge(L2);

    if (e1.first == e2.first || e1.first == e2.second)
    {
        return e1.first;
    }

    if (e1.second == e2.first || e1.second == e2.second)
    {
        return e1.second;
    }

    return constants::missing::uintValue;
}

void meshkernel::CasulliDeRefinement::CleanUpLink(Mesh2D& mesh, const UInt L)
{

    std::cout << " cleanup link: " << L << "  ";
    if (mesh.GetEdge(L).first != constants::missing::uintValue)
    {
        std::cout << mesh.Node(mesh.GetEdge(L).first).x << ",  " << mesh.Node(mesh.GetEdge(L).first).y << " ---  ";
    }
    else
    {
        std::cout << " null point ---  ";
    }

    if (mesh.GetEdge(L).second != constants::missing::uintValue)
    {
        std::cout << mesh.Node(mesh.GetEdge(L).second).x << ",  " << mesh.Node(mesh.GetEdge(L).second).y << "  ";
    }
    else
    {
        std::cout << " null point ---  ";
    }
    std::cout << std::endl;

    for (UInt i = 0; i < 2; ++i)
    {
        UInt k = EdgeNodeIndex(mesh.GetEdge(L), i); // i == 0 ? mesh.Edges(L).first : mesh.Edges(L).second;

        if (k == constants::missing::uintValue)
        {
            // already cleaned up
            continue;
        }

        // UInt N = mesh.m_nodesNumEdges[k];
        UInt j = constants::missing::uintValue; // mesh.m_nodesEdges[k].size();

        for (UInt jj = 0; jj < mesh.m_nodesNumEdges[k]; ++jj)
        // for (UInt jj = 0; jj < mesh.m_nodesEdges[k].size(); ++jj)
        {
            if (mesh.m_nodesEdges[k][jj] == L)
            {
                j = jj;
                break;
            }
        }

        if (j != constants::missing::uintValue)
        {
            // TODO use ShiftArray(i, mesh.m_XXX);
            // TODO check loop range
            for (UInt jj = j; jj < mesh.m_nodesNumEdges[k]; ++jj)
            {
                mesh.m_nodesEdges[k][jj] = mesh.m_nodesEdges[k][jj + 1];
            }
        }
        else
        {
            // Error?
            // call qnerror('cleanup_nod: link not found', ' ', ' ')
            return;
        }

        --mesh.m_nodesNumEdges[k];
    }

    mesh.GetEdge(L).first = constants::missing::uintValue;
    mesh.GetEdge(L).second = constants::missing::uintValue;
}

void meshkernel::CasulliDeRefinement::StoreMesh(Mesh2D& mesh, const UInt k,
                                                const std::vector<UInt>& kDirect,
                                                const std::vector<UInt>& kIndirect,
                                                std::vector<UInt>& savedNodes,
                                                std::vector<Edge>& savedEdges,
                                                std::vector<UInt>& savedCells)
{

    [[maybe_unused]] UInt maxCells = 1u + kDirect.size() + kIndirect.size();
    std::vector<UInt> allCellIds;
    allCellIds.reserve(maxCells);

    allCellIds.push_back(k);
    allCellIds.insert(allCellIds.end(), kDirect.begin(), kDirect.end());
    allCellIds.insert(allCellIds.end(), kIndirect.begin(), kIndirect.end());

    [[maybe_unused]] UInt numberOfNodes = std::accumulate(allCellIds.begin(), allCellIds.end(), 0,
                                                          [&grid = mesh](const UInt total, const UInt index)
                                                          {
                                                              return total + grid.m_numFacesNodes[index];
                                                          });

    [[maybe_unused]] UInt maxNodes = numberOfNodes;
    [[maybe_unused]] UInt maxLinks = numberOfNodes;

    savedNodes.reserve(maxNodes);
    savedEdges.reserve(maxLinks);
    savedCells.reserve(maxCells);

    // for (UInt i = 0; i < maxCells; ++)
    // {
    // }
}

void meshkernel::CasulliDeRefinement::InitialiseBoundaryNodes(Mesh2D& mesh, std::vector<NodeMask>& nodeMask)
{
    // Find nodes that lie on the boundary of the domain.
    for (UInt i = 0; i < mesh.GetNumEdges(); ++i)
    {
        UInt node1 = mesh.GetEdge(i).first;
        UInt node2 = mesh.GetEdge(i).second;

        if (mesh.m_edgesNumFaces[i] == 1)
        {
            if (nodeMask[node1] != NodeMask::Unassigned)
            {
                nodeMask[node1] = NodeMask::BoundaryNode;
            }

            if (nodeMask[node2] != NodeMask::Unassigned)
            {
                nodeMask[node2] = NodeMask::BoundaryNode;
            }
        }
    }
}

void meshkernel::CasulliDeRefinement::InitialiseCornerNodes(Mesh2D& mesh, std::vector<NodeMask>& nodeMask)
{
    for (UInt i = 0; i < mesh.GetNumNodes(); ++i)
    {

        for (UInt j = 0; j < mesh.m_nodesNumEdges[i]; ++j)
        {
            UInt edge1 = mesh.m_nodesEdges[i][j];

            if (mesh.m_edgesNumFaces[edge1] != 1)
            {
                continue;
            }

            UInt elementId = mesh.m_edgesFaces[edge1][0];
            UInt nodeCount = mesh.m_numFacesNodes[elementId];

            UInt faceEdgeIndex = 0;
            UInt edge2 = mesh.m_facesEdges[elementId][faceEdgeIndex];

            // Check the loop termination, especially the faceEdgeIndex < nodeCount - 1
            // Perhaps change to for loop checking the condition then break.
            while (((mesh.GetEdge(edge2).first != i && mesh.GetEdge(edge2).second != i) || edge2 == edge1) && faceEdgeIndex < nodeCount - 1)
            {
                ++faceEdgeIndex;
                edge2 = mesh.m_facesEdges[elementId][faceEdgeIndex];
            }

            if (mesh.m_edgesNumFaces[edge2] == 1)
            {
                if (nodeMask[i] > NodeMask::Unassigned)
                {
                    nodeMask[i] = NodeMask::CornerNode;
                    break;
                }
            }
        }
    }

    // Find included corner nodes
    for (UInt i = 0; i < mesh.GetNumNodes(); ++i)
    {
        if (nodeMask[i] != NodeMask::Unassigned && mesh.m_nodesTypes[i] == 3)
        {
            nodeMask[i] = NodeMask::CornerNode;
        }
    }
}

void meshkernel::CasulliDeRefinement::InitialiseFaceNodes(Mesh2D& mesh, std::vector<NodeMask>& nodeMask)
{

    std::vector<UInt> sharedFaces;
    std::vector<UInt> connectedNodes;
    std::vector<std::vector<UInt>> faceNodeMapping;

    for (UInt i = 0; i < mesh.GetNumNodes(); ++i)
    {
        if (nodeMask[i] == NodeMask::Unassigned)
        {
            continue;
        }

        if (mesh.m_nodesNumEdges[i] > 1)
        {
            FindPatchIds(mesh, i, sharedFaces, connectedNodes, faceNodeMapping);
        }
        else
        {
            nodeMask[i] = NodeMask::Unassigned;
            continue;
        }

        UInt elementCount = 0;

        for (UInt j = 0; j < sharedFaces.size(); ++j)
        {
            if (sharedFaces[j] != constants::missing::uintValue)
            {
                ++elementCount;
            }
        }

        if (elementCount == 0)
        {
            nodeMask[i] = NodeMask::Unassigned;
        }

        if (elementCount < mesh.m_nodesNumEdges[i] - 1 && nodeMask[i] == NodeMask::BoundaryNode)
        {
            nodeMask[i] = NodeMask::CornerNode;
        }

        if (elementCount > Mesh::m_maximumNumberOfEdgesPerNode && nodeMask[i] > NodeMask::Unassigned && nodeMask[i] < NodeMask::BoundaryNode)
        {
            nodeMask[i] = NodeMask::CornerNode;
        }
    }
}

std::vector<meshkernel::CasulliDeRefinement::NodeMask> meshkernel::CasulliDeRefinement::InitialiseNodeMask(Mesh2D& mesh, const Polygons& polygon)
{
    std::vector<NodeMask> nodeMask(10 * mesh.GetNumNodes(), NodeMask::Unassigned);

    // Find nodes that are inside the polygon.
    // If the polygon is empty then all nodes will be taken into account.

    for (UInt i = 0; i < mesh.GetNumNodes(); ++i)
    {
        auto [containsPoint, pointIndex] = polygon.IsPointInPolygons(mesh.Node(i));

        if (containsPoint)
        {
            nodeMask[i] = NodeMask::RegisteredNode;
        }
    }

    InitialiseBoundaryNodes(mesh, nodeMask);
    InitialiseCornerNodes(mesh, nodeMask);
    InitialiseFaceNodes(mesh, nodeMask);

    return nodeMask;
}

void meshkernel::CasulliDeRefinement::FindPatchIds(Mesh2D& mesh,
                                                   const UInt currentNode,
                                                   std::vector<UInt>& sharedFaces,
                                                   std::vector<UInt>& connectedNodes,
                                                   std::vector<std::vector<UInt>>& faceNodeMapping)
{
    sharedFaces.clear();
    connectedNodes.clear();
    faceNodeMapping.clear();

    if (currentNode >= mesh.GetNumNodes())
    {
        throw AlgorithmError("Node index out of range: {} >= {}", currentNode, mesh.GetNumNodes());
    }

    if (mesh.m_nodesNumEdges[currentNode] < 2)
    {
        return;
    }

    mesh.FindFacesConnectedToNode(currentNode, sharedFaces);

    // no shared face found
    if (sharedFaces.empty())
    {
        return;
    }

    mesh.GetConnectingNodes(currentNode, connectedNodes);
    mesh.FindNodesSharedByFaces(currentNode, sharedFaces, connectedNodes, faceNodeMapping);
}

bool meshkernel::CasulliDeRefinement::IsConvexCell(const Mesh2D& mesh, const UInt cell)
{
    constexpr double TOL = 0.000001;

    UInt ip1;
    UInt ip2;
    UInt N = mesh.m_numFacesNodes[cell];

    for (UInt i = 0; i < N; ++i)
    {
        ip1 = RotateIndex(i, N, true /* forward */);
        ip2 = RotateIndex(ip1, N, true /* forward */);

        UInt k1 = mesh.m_facesNodes[cell][i];
        UInt k2 = mesh.m_facesNodes[cell][ip1];
        UInt k3 = mesh.m_facesNodes[cell][ip2];

        double cosphi = NormalizedInnerProductTwoSegments(mesh.Node(k1), mesh.Node(k2), mesh.Node(k2), mesh.Node(k3), mesh.m_projection);
        if (std::abs(1.0 - std::abs(cosphi) < TOL))
        {
            return false;
        }
    }

    return true;
}
