//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2023.
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

#include <cmath>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Vector.hpp"

#include "MeshKernel/Utilities/LinearAlgebra.hpp"

#include "MeshKernel/CurvilinearGrid/CurvilinearGridSmoothness.hpp"

void meshkernel::CurvilinearGridSmoothness::Compute(const CurvilinearGrid& grid, const CurvilinearDirection direction, lin_alg::Matrix<double>& smoothness)
{
    lin_alg::ResizeAndFillMatrix(smoothness, grid.m_numM, grid.m_numN, false, constants::missing::doubleValue);

    if (direction == CurvilinearDirection::M)
    {
        for (UInt i = 1; i < grid.m_numM - 1; ++i)
        {
            for (UInt j = 0; j < grid.m_numN; ++j)
            {
                smoothness(i, j) = ComputeNodeSmoothness(grid.m_gridNodes(i - 1, j), grid.m_gridNodes(i, j), grid.m_gridNodes(i + 1, j));
            }
        }
    }
    else if (direction == CurvilinearDirection::N)
    {
        for (UInt i = 0; i < grid.m_numM; ++i)
        {
            for (UInt j = 1; j < grid.m_numN - 1; ++j)
            {
                smoothness(i, j) = ComputeNodeSmoothness(grid.m_gridNodes(i, j - 1), grid.m_gridNodes(i, j), grid.m_gridNodes(i, j + 1));
            }
        }
    }
    else
    {
        throw MeshKernelError("Unknown curvilinear direction values {} with integer value {}", CurvilinearDirectionToString(direction), static_cast<int>(direction));
    }
}

double meshkernel::CurvilinearGridSmoothness::ComputeNodeSmoothness(const Point& p0, const Point& p1, const Point& p2)
{
    double nodeSmoothness = constants::missing::doubleValue;

    if (p0.IsValid() && p1.IsValid() && p2.IsValid())
    {
        double diffX10 = p1.x - p0.x;
        double diffY10 = p1.y - p0.y;

        double diffX21 = p2.x - p1.x;
        double diffY21 = p2.y - p1.y;

        double lengthSquared10 = diffX10 * diffX10 + diffY10 * diffY10;
        double lengthSquared21 = diffX21 * diffX21 + diffY21 * diffY21;

        if (lengthSquared10 != 0.0 && lengthSquared21 != 0.0)
        {
            nodeSmoothness = std::sqrt(lengthSquared21 / lengthSquared10);

            if (nodeSmoothness < 1.0)
            {
                nodeSmoothness = 1.0 / nodeSmoothness;
            }
        }
        else
        {
            nodeSmoothness = 1.0;
        }
    }

    return nodeSmoothness;
}
