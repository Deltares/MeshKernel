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

#include <cmath>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Vector.hpp"

#include "MeshKernel/CurvilinearGrid/CurvilinearGridCurvature.hpp"

void meshkernel::CurvilinearGridCurvature::Compute(const CurvilinearGrid& grid, const CurvilinearDirection direction, lin_alg::Matrix<double>& curvature)
{
    lin_alg::ResizeAndFillMatrix(curvature, grid.NumM(), grid.NumN(), false, constants::missing::doubleValue);

    if (direction == CurvilinearDirection::M)
    {
        for (UInt i = 1; i < grid.NumM() - 1; ++i)
        {
            for (UInt j = 0; j < grid.NumN(); ++j)
            {
                curvature(i, j) = ComputeNodeCurvature(grid.m_gridNodes(i - 1, j), grid.m_gridNodes(i, j), grid.m_gridNodes(i + 1, j));
            }
        }
    }
    else if (direction == CurvilinearDirection::N)
    {
        for (UInt i = 0; i < grid.NumM(); ++i)
        {
            for (UInt j = 1; j < grid.NumN() - 1; ++j)
            {
                curvature(i, j) = ComputeNodeCurvature(grid.m_gridNodes(i, j - 1), grid.m_gridNodes(i, j), grid.m_gridNodes(i, j + 1));
            }
        }
    }
    else
    {
        throw MeshKernelError("Unknown curvilinear direction values {} with integer value {}", CurvilinearDirectionToString(direction), static_cast<int>(direction));
    }
}

double meshkernel::CurvilinearGridCurvature::ComputeNodeCurvature(const Point& p0, const Point& p1, const Point& p2)
{
    double nodeCurvature = constants::missing::doubleValue;

    if (p0.IsValid() && p1.IsValid() && p2.IsValid())
    {
        double diffX10 = p1.x - p0.x;
        double diffY10 = p1.y - p0.y;

        double diffX20 = p2.x - p0.x;
        double diffY20 = p2.y - p0.y;

        double diffX21 = p2.x - p1.x;
        double diffY21 = p2.y - p1.y;

        // Twice the area of the triangle formed by the points {p0, p1, p2}
        double area2 = diffX10 * diffY21 - diffY10 * diffX21;
        double distance01 = std::sqrt(diffX10 * diffX10 + diffY10 * diffY10);
        double radius = 999999.0;

        if (area2 != 0.0)
        {
            radius = distance01 * std::abs((diffX21 * diffX20 + diffY21 * diffY20) / area2);
        }

        if (radius == 0.0)
        {
            radius = 999999.0;
        }

        nodeCurvature = 1000.0 * std::abs(1.0 / radius);
    }

    return nodeCurvature;
}
