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

#include <vector>
#include <iostream>
#include <algorithm>
#include "Operations.cpp"
#include "Entities.hpp"
#include "Mesh.hpp"
#include "Polygons.hpp"
#include "CurvilinearGridFromPolygon.hpp"
#include "CurvilinearGrid.hpp"

MeshKernel::CurvilinearGridFromPolygon::CurvilinearGridFromPolygon() :
    m_polygon(nullptr)
{



}

MeshKernel::CurvilinearGridFromPolygon::CurvilinearGridFromPolygon(std::shared_ptr<Polygons> polygon) :
    m_polygon(polygon)
{

};

bool MeshKernel::CurvilinearGridFromPolygon::Compute( int firstNode, 
                                                      int secondNode, 
                                                      int thirdNode, 
                                                      bool useFourthSide, 
                                                      CurvilinearGrid& curvilinearGrid ) const
{
    if (m_polygon->m_indexses.size() < 0)
    {
        return true;
    }

    // for the current polygon find the number of nodes
    const auto start = m_polygon->m_indexses[0][0];
    const auto end = m_polygon->m_indexses[0][1];
    const int numPolygonNodes = end - start + 1;

    // get rid of size and orientation first part
    int  diffForward = secondNode - firstNode;
    if (diffForward < 0)
    {
        diffForward = diffForward + numPolygonNodes;
    }

    int  diffBackward = firstNode - secondNode;
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
        return false;
    }

    int maximumNumberOfNodes = std::max(numNNodes, numMNodes);
    std::vector<Point> sideOne(maximumNumberOfNodes, { doubleMissingValue, doubleMissingValue });
    std::vector<Point> sideTwo(maximumNumberOfNodes, { doubleMissingValue, doubleMissingValue });
    std::vector<Point> sideThree(maximumNumberOfNodes, { doubleMissingValue, doubleMissingValue });
    std::vector<Point> sideFour(maximumNumberOfNodes, { doubleMissingValue, doubleMissingValue });

    // Fill boundary coordinates
    auto assignPolygonPointsToSegment = [&](int nodeIndex, int numPointsSide, int direction, std::vector<Point>& sideToFill)
    {
        for (int i = 0; i < numPointsSide; i++)
        {
            sideToFill[i] = m_polygon->m_nodes[nodeIndex];
            nodeIndex = nodeIndex + direction;
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
            const double fac = double(i) / double((numNNodes - 1));
            sideOne[i] = sideOne[i] * (1.0 - fac) + sideOne[i] * fac;
        }
    }

    assignPolygonPointsToSegment(secondNode, numNNodes, direction, sideTwo);
    assignPolygonPointsToSegment(firstNode, numMNodes, direction, sideThree);
    assignPolygonPointsToSegment(fourthNode, numMNodes, -direction, sideFour);

    std::vector<std::vector<Point>> result;
    bool successful = InterpolateTransfinite( sideOne,
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

    return successful;
}