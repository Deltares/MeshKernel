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

#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Point.hpp"

namespace meshkernel
{

    /// @brief Print the (simplified) graph in a form that can be loaded into matlab/octave.
    ///
    /// Only nodes and node connectivity need be printed to visualise the graph.
    void Print(const std::vector<Point>& nodes, const std::vector<Edge>& edges, std::ostream& out = std::cout);

    /// @brief Print the (simplified) graph in a form that can be loaded into matlab/octave.
    ///
    /// Only nodes and node connectivity need be printed to visualise the graph.
    void Print(const std::vector<double>& xNodes,
               const std::vector<double>& yNodes,
               const std::vector<int>& edges, std::ostream& out = std::cout);

    /// @brief Save the mesh data in a vtk file format
    ///
    /// @note saves only triangle and quadrilateral elements.
    void SaveVtk(const std::vector<Point>& nodes, const std::vector<std::vector<UInt>>& faces, const std::string& fileName);

} // namespace meshkernel
