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

#include "MeshKernel/Utilities/Utilities.hpp"
#include "MeshKernel/Exceptions.hpp"

void meshkernel::Print(const std::vector<Point>& nodes, const std::vector<Edge>& edges, std::ostream& out)
{
    out << "nullId = " << constants::missing::uintValue << ";" << std::endl;
    out << "nullValue = " << constants::missing::doubleValue << ";" << std::endl;
    out << "nodex = zeros ( " << nodes.size() << ", 1);" << std::endl;
    out << "nodey = zeros ( " << nodes.size() << ", 1);" << std::endl;
    out << "edges = zeros ( " << edges.size() << ", 2);" << std::endl;

    for (UInt i = 0; i < nodes.size(); ++i)
    {
        out << "nodex (" << i + 1 << " ) = " << nodes[i].x << ";" << std::endl;
    }

    for (UInt i = 0; i < nodes.size(); ++i)
    {
        out << "nodey (" << i + 1 << " ) = " << nodes[i].y << ";" << std::endl;
    }

    out << "edges = zeros ( " << edges.size() << ", 2 );" << std::endl;

    for (UInt i = 0; i < edges.size(); ++i)
    {
        out << "edges ( " << i + 1 << ", 1 ) = " << edges[i].first + 1 << ";" << std::endl;
        out << "edges ( " << i + 1 << ", 2 ) = " << edges[i].second + 1 << ";" << std::endl;
    }
}

void meshkernel::Print(const std::vector<double>& xNodes,
                       const std::vector<double>& yNodes,
                       const std::vector<int>& edges,
                       std::ostream& out)
{
    // xNodes and yNodes should be the same size
    if (xNodes.size() != yNodes.size())
    {
        throw ConstraintError("x-node and y-nodes are not the same size, {} /= {}", xNodes.size(), yNodes.size());
    }

    out << "nullId = " << constants::missing::uintValue << ";" << std::endl;
    out << "nullValue = " << constants::missing::doubleValue << ";" << std::endl;
    out << "nodex = zeros ( " << xNodes.size() << ", 1);" << std::endl;
    out << "nodey = zeros ( " << xNodes.size() << ", 1);" << std::endl;
    out << "edges = zeros ( " << edges.size() << ", 2);" << std::endl;

    for (UInt i = 0; i < xNodes.size(); ++i)
    {
        out << "nodex (" << i + 1 << " ) = " << xNodes[i] << ";" << std::endl;
    }

    for (UInt i = 0; i < xNodes.size(); ++i)
    {
        out << "nodey (" << i + 1 << " ) = " << yNodes[i] << ";" << std::endl;
    }

    out << "edges = zeros ( " << edges.size() / 2 << ", 2 );" << std::endl;

    for (UInt i = 0; i < edges.size() / 2; ++i)
    {
        out << "edges ( " << i + 1 << ", 1 ) = " << edges[2 * i] + 1 << ";" << std::endl;
        out << "edges ( " << i + 1 << ", 2 ) = " << edges[2 * i + 1] + 1 << ";" << std::endl;
    }
}
