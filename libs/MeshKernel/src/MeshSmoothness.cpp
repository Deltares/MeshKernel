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

#include "MeshKernel/MeshSmoothness.hpp"
#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Exceptions.hpp"

std::vector<double> meshkernel::MeshSmoothness::Compute(const Mesh2D& mesh) const
{

    std::vector<double> smoothness(mesh.GetNumEdges(), constants::missing::doubleValue);
    Compute(mesh, smoothness);

    return smoothness;
}

void meshkernel::MeshSmoothness::Compute(const Mesh2D& mesh, std::span<double> smoothness) const
{

    if (smoothness.size() != mesh.GetNumEdges())
    {
        throw ConstraintError("array for smoothness values is not the correct size");
    }

    const UInt numEdges = mesh.GetNumEdges();

    for (UInt e = 0; e < numEdges; e++)
    {
        auto val = constants::missing::doubleValue;

        const auto [firstNode, secondNode] = mesh.GetEdge(e);

        const auto firstFaceIndex = mesh.m_edgesFaces[e][0];
        const auto secondFaceIndex = mesh.m_edgesFaces[e][1];

        if (firstNode != constants::missing::uintValue &&
            secondNode != constants::missing::uintValue &&
            firstFaceIndex != constants::missing::uintValue &&
            secondFaceIndex != constants::missing::uintValue && !mesh.IsEdgeOnBoundary(e))
        {
            const auto leftFaceArea = mesh.m_faceArea[firstFaceIndex];
            const auto rightFaceArea = mesh.m_faceArea[secondFaceIndex];

            if (leftFaceArea > m_minimumCellArea && rightFaceArea > m_minimumCellArea)
            {
                val = rightFaceArea / leftFaceArea;

                if (val < 1.0)
                {
                    val = 1.0 / val;
                }
            }
        }

        smoothness[e] = val;
    }
}
