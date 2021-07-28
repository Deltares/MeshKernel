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

#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/Polygons.hpp>

void meshkernel::CurvilinearGrid::Delete(Polygons const& polygons, size_t polygonIndex)
{
    // no polygons available
    if (polygons.IsEmpty())
    {
        return;
    }

    // no grid, return
    if (m_grid.empty())
    {
        return;
    }

    const auto numN = m_grid.size() - 1;
    const auto numM = m_grid[0].size() - 1;

    std::vector<std::vector<bool>> nodeBasedMask(numN, std::vector<bool>(numM, false));
    std::vector<std::vector<bool>> faceBasedMask(numN - 1, std::vector<bool>(numM - 1, true));
    // Mark points inside a polygonIndex
    for (auto n = 0; n < numN; ++n)
    {
        for (auto m = 0; m < numM; ++m)
        {
            const auto isInPolygon = polygons.IsPointInPolygon(m_grid[n][m], polygonIndex);
            if (isInPolygon)
            {
                nodeBasedMask[n][m] = true;
            }
        }
    }

    // Mark faces when all nodes are inside
    for (auto n = 0; n < numN - 1; ++n)
    {
        for (auto m = 0; m < numM - 1; ++m)
        {
            if (!nodeBasedMask[n][m] ||
                !nodeBasedMask[n + 1][m] ||
                !nodeBasedMask[n][m + 1] ||
                !nodeBasedMask[n + 1][m + 1])
            {
                faceBasedMask[n][m] = false;
            }
        }
    }

    // Mark only the nodes of faces completely included in the polygonIndex
    std::fill(nodeBasedMask.begin(), nodeBasedMask.end(), std::vector<bool>(numM, false));
    for (auto n = 0; n < numN - 1; ++n)
    {
        for (auto m = 0; m < numM - 1; ++m)
        {
            if (faceBasedMask[n][m])
            {
                nodeBasedMask[n][m] = true;
                nodeBasedMask[n + 1][m] = true;
                nodeBasedMask[n][m + 1] = true;
                nodeBasedMask[n + 1][m + 1] = true;
            }
        }
    }

    // mark points inside a polygonIndex
    for (auto n = 0; n < numN; ++n)
    {
        for (auto m = 0; m < numM; ++m)
        {
            if (!nodeBasedMask[n][m])
            {
                m_grid[n][m].x = doubleMissingValue;
                m_grid[n][m].y = doubleMissingValue;
            }
        }
    }
}