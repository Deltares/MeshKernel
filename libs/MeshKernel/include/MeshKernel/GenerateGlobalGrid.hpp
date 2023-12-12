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

#pragma once
#include <array>
#include <utility>

#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Point.hpp"

namespace meshkernel
{

    class GenerateGlobalGrid
    {
    public:
        // TODO may be better to return mesh rather than void
        // nx and ny are number of elements in each direction
        static void Compute(const UInt nx, const UInt ny, Mesh2D& mesh);

    private:
        static double getDeltaY(const double y, const double deltaX);

        static void isNodeDB(const Mesh& mesh, const Point& x, UInt& kp);

        static void addMaze(Mesh& mesh, const std::array<Point, 8>& points, const double ySign, const UInt pointSize, const bool jafive);

        static void mergenodesinpolygon(Mesh2D& mesh);
    };

} // namespace meshkernel
