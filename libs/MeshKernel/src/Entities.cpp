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

#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Exceptions.hpp"

std::tuple<int, double*, double*> meshkernel::ConvertFromNodesVector(const std::vector<Point>& nodes)
{

    double* xNodes = nullptr;
    double* yNodes = nullptr;

    try
    {
        if (nodes.empty())
        {
            return {0, nullptr, nullptr};
        }

        xNodes = new double[nodes.size()];
        yNodes = new double[nodes.size()];

        for (size_t i = 0; i < nodes.size(); ++i)
        {
            xNodes[i] = nodes[i].x;
            yNodes[i] = nodes[i].y;
        }

        return {static_cast<int>(nodes.size()), xNodes, yNodes};
    }
    catch (...)
    {
        if (xNodes != nullptr)
        {
            delete[] xNodes;
        }
        if (yNodes != nullptr)
        {
            delete[] yNodes;
        }

        throw;
    }
}
