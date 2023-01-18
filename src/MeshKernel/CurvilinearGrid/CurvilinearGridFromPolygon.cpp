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

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridFromPolygon.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Polygons.hpp>

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridFromPolygon;

CurvilinearGridFromPolygon::CurvilinearGridFromPolygon(std::shared_ptr<Polygons> polygon) : m_polygon(polygon) {}

CurvilinearGrid CurvilinearGridFromPolygon::Compute(size_t firstNode,
                                                    size_t secondNode,
                                                    size_t thirdNode,
                                                    bool useFourthSide) const
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
    auto const& [start, end] = m_polygon->OuterIndices(0);

    if (end <= start)
    {
        throw std::invalid_argument("CurvilinearGridFromPolygon::CurvilinearGridFromPolygon: Not enough points in polygon.");
    }

    const size_t numPolygonNodes = end - start + 1;

    // get rid of size and orientation first part
    size_t diffForward;
    if (firstNode > secondNode)
    {
        diffForward = secondNode + numPolygonNodes - firstNode;
    }
    else
    {
        diffForward = secondNode - firstNode;
    }

    size_t diffBackward;
    if (secondNode > firstNode)
    {
        diffBackward = firstNode + numPolygonNodes - secondNode;
    }
    else
    {
        diffBackward = firstNode - secondNode;
    }

    int direction;
    size_t numMNodes;
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
    if (secondNode > thirdNode)
    {
        diffForward = thirdNode + numPolygonNodes - secondNode;
    }
    else
    {
        diffForward = thirdNode - secondNode;
    }

    if (thirdNode > secondNode)
    {
        diffBackward = secondNode + numPolygonNodes - thirdNode;
    }
    else
    {
        diffBackward = secondNode - thirdNode;
    }

    size_t numNNodes;
    if (direction == 1)
    {
        numNNodes = diffForward + 1;
    }
    else
    {
        numNNodes = diffBackward + 1;
    }

    // get the fourth node
    auto fourthNode = thirdNode + direction * (numMNodes - 1);
    if (fourthNode < start)
    {
        fourthNode += numPolygonNodes;
    }
    if (fourthNode >= numPolygonNodes)
    {
        fourthNode -= numPolygonNodes;
    }

    size_t numRequiredPoints;
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

    const auto maximumNumberOfNodes = std::max(numNNodes, numMNodes);
    std::vector<Point> sideOne(maximumNumberOfNodes, {doubleMissingValue, doubleMissingValue});
    std::vector<Point> sideTwo(maximumNumberOfNodes, {doubleMissingValue, doubleMissingValue});
    std::vector<Point> sideThree(maximumNumberOfNodes, {doubleMissingValue, doubleMissingValue});
    std::vector<Point> sideFour(maximumNumberOfNodes, {doubleMissingValue, doubleMissingValue});

    // Fill boundary coordinates
    auto assignPolygonPointsToSegment = [this, &start, &end, &numPolygonNodes](size_t nodeIndex, size_t numPointsSide, int dir, std::vector<Point>& sideToFill)
    {
        for (size_t i = 0; i < numPointsSide; i++)
        {
            sideToFill[i] = m_polygon->Node(nodeIndex);

            if ((nodeIndex == 0 && dir == -1) || nodeIndex + dir < start)
            {
                nodeIndex = nodeIndex + numPolygonNodes + dir;
            }
            else if (nodeIndex + dir > end)
            {
                nodeIndex = nodeIndex + dir - numPolygonNodes;
            }
            else
            {
                nodeIndex = nodeIndex + dir;
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
        for (size_t i = 0; i < numNNodes; i++)
        {
            const double fac = static_cast<double>(i) / static_cast<double>(numNNodes - 1);
            sideOne[i] = m_polygon->Node(firstNode) * (1.0 - fac) +
                         m_polygon->Node(fourthNode) * fac;
        }
    }

    assignPolygonPointsToSegment(secondNode, numNNodes, direction, sideTwo);
    assignPolygonPointsToSegment(firstNode, numMNodes, direction, sideThree);
    assignPolygonPointsToSegment(fourthNode, numMNodes, -direction, sideFour);

    Projection const polygonProjection = m_polygon->GetProjection();

    const auto result = DiscretizeTransfinite(sideOne, sideTwo, sideThree, sideFour,
                                              polygonProjection, numMNodes - 1, numNNodes - 1);

    // Assign the points to the curvilinear grid
    std::vector<std::vector<Point>> gridNodes(numMNodes, std::vector<Point>(numNNodes));
    for (size_t i = 0; i < numMNodes; i++)
    {
        for (size_t j = 0; j < numNNodes; j++)
        {
            gridNodes[i][j] = result[i][j];
        }
    }
    return CurvilinearGrid(std::move(gridNodes), polygonProjection);
}

CurvilinearGrid CurvilinearGridFromPolygon::Compute(size_t firstNode,
                                                    size_t secondNode,
                                                    size_t thirdNode) const
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
    const auto& [start, end] = m_polygon->OuterIndices(0);

    if (end <= start)
    {
        throw std::invalid_argument("CurvilinearGridFromPolygon::Compute: Not enough points in polygon.");
    }

    const auto numPolygonNodes = end - start + 1;

    // get rid of size and orientation first part
    size_t numPointsFirstSide;
    if (firstNode > secondNode)
    {
        numPointsFirstSide = secondNode + numPolygonNodes - firstNode;
    }
    else
    {
        numPointsFirstSide = secondNode - firstNode;
    }

    size_t numPointsSecondSide;
    if (secondNode > thirdNode)
    {
        numPointsSecondSide = thirdNode + numPolygonNodes - secondNode;
    }
    else
    {
        numPointsSecondSide = thirdNode - secondNode;
    }

    const auto numPointsThirdSide = numPolygonNodes - (numPointsFirstSide + numPointsSecondSide);
    const auto blockSize = static_cast<size_t>(static_cast<double>(numPointsFirstSide + numPointsSecondSide + numPointsThirdSide) * 0.5);

    if (numPointsThirdSide >= blockSize || numPointsSecondSide >= blockSize || numPointsFirstSide >= blockSize)
    {
        throw std::invalid_argument("CurvilinearGridFromPolygon::Compute: The block size is less than the number of points.");
    }

    const auto n1 = blockSize - numPointsThirdSide;
    const auto n2 = blockSize - numPointsSecondSide;
    const auto n3 = blockSize - numPointsFirstSide;

    // compute the midpoint

    size_t firstSideMiddlePoint = firstNode + n1;
    if (firstSideMiddlePoint >= numPolygonNodes)
    {
        firstSideMiddlePoint = firstSideMiddlePoint - numPolygonNodes;
    }
    size_t secondSideMiddlePoint = secondNode + n3;
    if (secondSideMiddlePoint >= numPolygonNodes)
    {
        secondSideMiddlePoint = secondSideMiddlePoint - numPolygonNodes;
    }
    size_t thirdSideMiddlePoint = thirdNode + n2;
    if (thirdSideMiddlePoint >= numPolygonNodes)
    {
        thirdSideMiddlePoint = thirdSideMiddlePoint - numPolygonNodes;
    }

    // set dimensions of blocks
    std::vector<size_t> numM{n1, n3, n2};
    std::vector<size_t> numN{n3, n2, n1};

    // set pointers of block corners
    std::vector<size_t> cornerPoints{firstNode, secondNode, thirdNode};
    std::vector<size_t> iLeft{thirdSideMiddlePoint, firstSideMiddlePoint, secondSideMiddlePoint};
    std::vector<size_t> iRight{firstSideMiddlePoint, secondSideMiddlePoint, thirdSideMiddlePoint};

    // compute triangle middle point
    const auto xia = static_cast<double>(n1) / static_cast<double>(numPointsFirstSide);
    const auto xib = static_cast<double>(n2) / static_cast<double>(numPointsSecondSide);
    const auto xic = static_cast<double>(n3) / static_cast<double>(numPointsThirdSide);

    auto const& polygonNodes = m_polygon->Nodes();

    const auto triangleCenter = ((polygonNodes[firstNode] * (1.0 - xia) + polygonNodes[secondNode] * xia) * xic + polygonNodes[thirdNode] * (1.0 - xic) +
                                 (polygonNodes[secondNode] * (1.0 - xib) + polygonNodes[thirdNode] * xib) * xia + polygonNodes[firstNode] * (1.0 - xia) +
                                 (polygonNodes[thirdNode] * (1.0 - xic) + polygonNodes[firstNode] * xic) * xib + polygonNodes[secondNode] * (1.0 - xib)) *
                                oneThird;

    const auto maxM = *std::max_element(numM.begin(), numM.end());
    const auto maxN = *std::max_element(numN.begin(), numN.end());
    const auto maximumNumberOfNodes = std::max(maxM, maxN) + 1;
    std::vector<Point> sideOne(maximumNumberOfNodes, {doubleMissingValue, doubleMissingValue});
    std::vector<Point> sideTwo(maximumNumberOfNodes, {doubleMissingValue, doubleMissingValue});
    std::vector<Point> sideThree(maximumNumberOfNodes, {doubleMissingValue, doubleMissingValue});
    std::vector<Point> sideFour(maximumNumberOfNodes, {doubleMissingValue, doubleMissingValue});

    std::vector<std::vector<Point>> gridNodes(n1 + n3 + 1, std::vector<Point>(n2 + n3 + 1));

    Projection const polygonProjection = m_polygon->GetProjection();

    for (size_t t = 0; t < numNodesInTriangle; ++t)
    {
        std::fill(sideOne.begin(), sideOne.end(), Point{doubleMissingValue, doubleMissingValue});
        std::fill(sideTwo.begin(), sideTwo.end(), Point{doubleMissingValue, doubleMissingValue});
        std::fill(sideThree.begin(), sideThree.end(), Point{doubleMissingValue, doubleMissingValue});
        std::fill(sideFour.begin(), sideFour.end(), Point{doubleMissingValue, doubleMissingValue});

        // backward
        auto cornerIndex = cornerPoints[t];
        for (size_t i = 0; i < numN[t] + 1; ++i)
        {
            sideOne[i] = polygonNodes[cornerIndex];
            if (cornerIndex == 0 || cornerIndex < start)
            {
                cornerIndex = cornerIndex + numPolygonNodes - 1;
            }
            else if (cornerIndex > end)
            {
                cornerIndex = cornerIndex - numPolygonNodes - 1;
            }
            else
            {
                cornerIndex -= 1;
            }
        }

        // forward
        cornerIndex = cornerPoints[t];
        for (size_t i = 0; i < numM[t] + 1; ++i)
        {
            sideThree[i] = polygonNodes[cornerIndex];
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
        for (size_t i = 0; i < numM[t] + 1; ++i)
        {
            double localXia = static_cast<double>(i) / static_cast<double>(numM[t]);
            sideFour[i] = polygonNodes[iLeft[t]] * (1.0 - localXia) + triangleCenter * localXia;
        }

        // fill side two
        for (size_t i = 0; i < numN[t] + 1; ++i)
        {
            double localXia = static_cast<double>(i) / static_cast<double>(numN[t]);
            sideTwo[i] = polygonNodes[iRight[t]] * (1.0 - localXia) + triangleCenter * localXia;
        }

        const auto result = DiscretizeTransfinite(sideOne,
                                                  sideTwo,
                                                  sideThree,
                                                  sideFour,
                                                  polygonProjection,
                                                  numM[t],
                                                  numN[t]);

        // add to grid
        if (t == 0)
        {
            for (size_t i = 0; i < result.size(); ++i)
            {
                for (size_t j = 0; j < result[0].size(); ++j)
                {
                    gridNodes[i][j] = result[i][j];
                }
            }
        }
        if (t == 1)
        {
            for (size_t i = 0; i < result.size(); ++i)
            {
                for (size_t j = 0; j < result[0].size(); ++j)
                {
                    const auto iIndex = n1 + n3 - i;
                    const auto jIndex = n2 + n3 - j;
                    gridNodes[iIndex][jIndex] = result[i][j];
                }
            }
        }
        if (t == 2)
        {
            for (size_t i = 0; i < result[0].size(); ++i)
            {
                for (size_t j = 0; j < result.size(); ++j)
                {
                    const auto jIndex = n2 + n3 - j;
                    gridNodes[i][jIndex] = result[j][i];
                }
            }
        }
    }

    return CurvilinearGrid(std::move(gridNodes), polygonProjection);
}
