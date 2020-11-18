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

#include <MeshKernel/CurvilinearGridFromPolygon.hpp>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Operations.cpp>
#include <MeshKernel/Polygons.hpp>

meshkernel::CurvilinearGridFromPolygon::CurvilinearGridFromPolygon(std::shared_ptr<Polygons> polygon) : m_polygon(polygon){};

void meshkernel::CurvilinearGridFromPolygon::Compute(int firstNode,
                                                     int secondNode,
                                                     int thirdNode,
                                                     bool useFourthSide,
                                                     CurvilinearGrid& curvilinearGrid) const
{
    if (m_polygon->IsEmpty())
    {
        throw std::invalid_argument("CurvilinearGridFromPolygon::CurvilinearGridFromPolygon: The polygon contains no nodes.");
    }

    const auto areNodesValid = firstNode != secondNode &&
                               secondNode != thirdNode &&
                               firstNode != thirdNode;
    if (!areNodesValid)
    {
        throw std::invalid_argument("CurvilinearGridFromPolygon::CurvilinearGridFromPolygon: Invalid nodes.");
    }

    // for the current polygon find the number of nodes
    const auto start = m_polygon->m_indices[0][0];
    const auto end = m_polygon->m_indices[0][1];
    const int numPolygonNodes = end - start + 1;

    // get rid of size and orientation first part
    int diffForward = secondNode - firstNode;
    if (diffForward < 0)
    {
        diffForward = diffForward + numPolygonNodes;
    }

    int diffBackward = firstNode - secondNode;
    if (diffBackward < 0)
    {
        diffBackward = diffBackward + numPolygonNodes;
    }

    int direction;
    int numMNodes;
    if (diffForward <= diffBackward)
    {
        direction = 1;
        numMNodes = diffForward + 1;
    }
    else
    {

        direction = -1;
        numMNodes = diffBackward + 1;
    }

    // get rid of size and orientation second part

    diffForward = thirdNode - secondNode;
    if (diffForward < 0)
    {
        diffForward = diffForward + numPolygonNodes;
    }

    diffBackward = secondNode - thirdNode;
    if (diffBackward < 0)
    {
        diffBackward = diffBackward + numPolygonNodes;
    }

    int numNNodes;
    if (direction == 1)
    {
        numNNodes = diffForward + 1;
    }
    else
    {
        numNNodes = diffBackward + 1;
    }

    // get the fourth node
    int fourthNode;
    fourthNode = thirdNode + direction * (numMNodes - 1);
    if (fourthNode < start)
    {
        fourthNode += numPolygonNodes;
    }
    if (fourthNode >= numPolygonNodes)
    {
        fourthNode -= numPolygonNodes;
    }

    int numRequiredPoints = 0;
    if (useFourthSide)
    {
        numRequiredPoints = 2 * (numMNodes - 1) + 2 * (numNNodes - 1);
    }
    else
    {
        numRequiredPoints = 1 + 2 * (numMNodes - 1) + (numNNodes - 1);
    }

    if (numRequiredPoints > numPolygonNodes)
    {
        throw std::invalid_argument("CurvilinearGridFromPolygon::CurvilinearGridFromPolygon: The polygon does not contain enough nodes to compute the curvilinear grid.");
    }

    int maximumNumberOfNodes = std::max(numNNodes, numMNodes);
    std::vector<Point> sideOne(maximumNumberOfNodes, {doubleMissingValue, doubleMissingValue});
    std::vector<Point> sideTwo(maximumNumberOfNodes, {doubleMissingValue, doubleMissingValue});
    std::vector<Point> sideThree(maximumNumberOfNodes, {doubleMissingValue, doubleMissingValue});
    std::vector<Point> sideFour(maximumNumberOfNodes, {doubleMissingValue, doubleMissingValue});

    // Fill boundary coordinates
    auto assignPolygonPointsToSegment = [this, start, end, numPolygonNodes](int nodeIndex, int numPointsSide, int dir, std::vector<Point>& sideToFill) {
        for (int i = 0; i < numPointsSide; i++)
        {
            sideToFill[i] = m_polygon->m_nodes[nodeIndex];
            nodeIndex = nodeIndex + dir;
            if (nodeIndex < start)
            {
                nodeIndex += numPolygonNodes;
            }
            if (nodeIndex > end)
            {
                nodeIndex -= numPolygonNodes;
            }
        }
    };

    if (useFourthSide)
    {
        assignPolygonPointsToSegment(firstNode, numNNodes, -direction, sideOne);
    }
    else
    {
        // Interpolate fourth side
        for (int i = 0; i < numNNodes; i++)
        {
            const double fac = double(i) / double(numNNodes - 1);
            sideOne[i] = m_polygon->m_nodes[firstNode] * (1.0 - fac) +
                         m_polygon->m_nodes[fourthNode] * fac;
        }
    }

    assignPolygonPointsToSegment(secondNode, numNNodes, direction, sideTwo);
    assignPolygonPointsToSegment(firstNode, numMNodes, direction, sideThree);
    assignPolygonPointsToSegment(fourthNode, numMNodes, -direction, sideFour);

    std::vector<std::vector<Point>> result;
    InterpolateTransfinite(sideOne,
                           sideTwo,
                           sideThree,
                           sideFour,
                           m_polygon->m_projection,
                           numMNodes - 1,
                           numNNodes - 1,
                           result);

    // Assign the points to the curvilinear grid
    curvilinearGrid.Set(numMNodes, numNNodes);
    for (int i = 0; i < numMNodes; i++)
    {
        for (int j = 0; j < numNNodes; j++)
        {
            curvilinearGrid.m_grid[i][j] = result[i][j];
        }
    }
}

void meshkernel::CurvilinearGridFromPolygon::Compute(int firstNode,
                                                     int secondNode,
                                                     int thirdNode,
                                                     CurvilinearGrid& curvilinearGrid) const
{
    if (m_polygon->IsEmpty())
    {
        throw std::invalid_argument("CurvilinearGridFromPolygon::Compute: The polygon contains no nodes.");
    }

    const auto areNodesValid = firstNode != secondNode &&
                               secondNode != thirdNode &&
                               firstNode != thirdNode;

    if (!areNodesValid)
    {
        throw std::invalid_argument("CurvilinearGridFromPolygon::Compute: Invalid nodes.");
    }

    // for the current polygon find the number of nodes
    const auto start = m_polygon->m_indices[0][0];
    const auto end = m_polygon->m_indices[0][1];
    const int numPolygonNodes = end - start + 1;

    // get rid of size and orientation first part
    int numPointsFirstSide = secondNode - firstNode;
    if (numPointsFirstSide < 0)
    {
        numPointsFirstSide = numPointsFirstSide + numPolygonNodes;
    }

    int numPointsSecondSide = thirdNode - secondNode;
    if (numPointsSecondSide < 0)
    {
        numPointsSecondSide = numPointsSecondSide + numPolygonNodes;
    }

    int numPointsThirdSide = numPolygonNodes - (numPointsFirstSide + numPointsSecondSide);
    int blockSize = (numPointsFirstSide + numPointsSecondSide + numPointsThirdSide) / 2;

    int n1 = blockSize - numPointsThirdSide;
    int n2 = blockSize - numPointsSecondSide;
    int n3 = blockSize - numPointsFirstSide;

    if (n1 < 1 || n2 < 1 || n3 < 1)
    {
        throw std::invalid_argument("CurvilinearGridFromPolygon::Compute: The block size is less than the number of points.");
    }

    // compute the midpoint

    int firstSideMiddlePoint = firstNode + n1;
    if (firstSideMiddlePoint >= numPolygonNodes)
    {
        firstSideMiddlePoint = firstSideMiddlePoint - numPolygonNodes;
    }
    int secondSideMiddlePoint = secondNode + n3;
    if (secondSideMiddlePoint >= numPolygonNodes)
    {
        secondSideMiddlePoint = secondSideMiddlePoint - numPolygonNodes;
    }
    int thirdSideMiddlePoint = thirdNode + n2;
    if (thirdSideMiddlePoint >= numPolygonNodes)
    {
        thirdSideMiddlePoint = thirdSideMiddlePoint - numPolygonNodes;
    }

    // set dimensions of blocks
    std::vector<int> numM{n1, n3, n2};
    std::vector<int> numN{n3, n2, n1};

    // set pointers of block corners
    std::vector<int> cornerPoints{firstNode, secondNode, thirdNode};
    std::vector<int> iLeft{thirdSideMiddlePoint, firstSideMiddlePoint, secondSideMiddlePoint};
    std::vector<int> iRight{firstSideMiddlePoint, secondSideMiddlePoint, thirdSideMiddlePoint};

    // compute triangle middle point
    const auto xia = double(n1) / double(numPointsFirstSide);
    const auto xib = double(n2) / double(numPointsSecondSide);
    const auto xic = double(n3) / double(numPointsThirdSide);

    auto triangleCenter = ((m_polygon->m_nodes[firstNode] * (1.0 - xia) + m_polygon->m_nodes[secondNode] * xia) * xic + m_polygon->m_nodes[thirdNode] * (1.0 - xic) +
                           (m_polygon->m_nodes[secondNode] * (1.0 - xib) + m_polygon->m_nodes[thirdNode] * xib) * xia + m_polygon->m_nodes[firstNode] * (1.0 - xia) +
                           (m_polygon->m_nodes[thirdNode] * (1.0 - xic) + m_polygon->m_nodes[firstNode] * xic) * xib + m_polygon->m_nodes[secondNode] * (1.0 - xib)) /
                          3.0;

    const auto maxM = *std::max_element(numM.begin(), numM.end());
    const auto maxN = *std::max_element(numN.begin(), numN.end());
    const auto maximumNumberOfNodes = std::max(maxM, maxN) + 1;
    std::vector<Point> sideOne(maximumNumberOfNodes, {doubleMissingValue, doubleMissingValue});
    std::vector<Point> sideTwo(maximumNumberOfNodes, {doubleMissingValue, doubleMissingValue});
    std::vector<Point> sideThree(maximumNumberOfNodes, {doubleMissingValue, doubleMissingValue});
    std::vector<Point> sideFour(maximumNumberOfNodes, {doubleMissingValue, doubleMissingValue});

    curvilinearGrid.Set(n1 + n3 + 1, n2 + n3 + 1);
    for (int t = 0; t < 3; ++t)
    {
        std::fill(sideOne.begin(), sideOne.end(), Point{doubleMissingValue, doubleMissingValue});
        std::fill(sideTwo.begin(), sideTwo.end(), Point{doubleMissingValue, doubleMissingValue});
        std::fill(sideThree.begin(), sideThree.end(), Point{doubleMissingValue, doubleMissingValue});
        std::fill(sideFour.begin(), sideFour.end(), Point{doubleMissingValue, doubleMissingValue});

        // backward
        int cornerIndex = cornerPoints[t];
        for (int i = 0; i < numN[t] + 1; ++i)
        {
            sideOne[i] = m_polygon->m_nodes[cornerIndex];
            cornerIndex -= 1;
            if (cornerIndex < start)
            {
                cornerIndex = cornerIndex + numPolygonNodes;
            }
            if (cornerIndex > end)
            {
                cornerIndex = cornerIndex - numPolygonNodes;
            }
        }

        // forward
        cornerIndex = cornerPoints[t];
        for (int i = 0; i < numM[t] + 1; ++i)
        {
            sideThree[i] = m_polygon->m_nodes[cornerIndex];
            cornerIndex += 1;
            if (cornerIndex < start)
            {
                cornerIndex = cornerIndex + numPolygonNodes;
            }
            if (cornerIndex > end)
            {
                cornerIndex = cornerIndex - numPolygonNodes;
            }
        }

        // fill side four
        for (int i = 0; i < numM[t] + 1; ++i)
        {
            double localXia = double(i) / double(numM[t]);
            sideFour[i] = m_polygon->m_nodes[iLeft[t]] * (1.0 - localXia) + triangleCenter * localXia;
        }

        // fill side two
        for (int i = 0; i < numN[t] + 1; ++i)
        {
            double localXia = double(i) / double(numN[t]);
            sideTwo[i] = m_polygon->m_nodes[iRight[t]] * (1.0 - localXia) + triangleCenter * localXia;
        }

        std::vector<std::vector<Point>> result;
        InterpolateTransfinite(sideOne,
                               sideTwo,
                               sideThree,
                               sideFour,
                               m_polygon->m_projection,
                               numM[t],
                               numN[t],
                               result);

        // add to grid
        if (t == 0)
        {
            for (int i = 0; i < result.size(); ++i)
            {
                for (int j = 0; j < result[0].size(); ++j)
                {
                    curvilinearGrid.m_grid[i][j] = result[i][j];
                }
            }
        }
        if (t == 1)
        {
            for (int i = 0; i < result.size(); ++i)
            {
                for (int j = 0; j < result[0].size(); ++j)
                {
                    int iIndex = n1 + n3 - i;
                    int jIndex = n2 + n3 - j;
                    curvilinearGrid.m_grid[iIndex][jIndex] = result[i][j];
                }
            }
        }
        if (t == 2)
        {
            for (int i = 0; i < result[0].size(); ++i)
            {
                for (int j = 0; j < result.size(); ++j)
                {
                    int jIndex = n2 + n3 - j;
                    curvilinearGrid.m_grid[i][jIndex] = result[j][i];
                }
            }
        }
    }
}
