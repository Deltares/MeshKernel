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

#pragma once
#include <vector>
#include <memory>

#include "Entities.hpp"

namespace meshkernel
{
    class Mesh2D;

    
    /// \namespace CutCellNodeClasses
    /// @brief Contains the logic of the C++ static library
    namespace CutCellNodeClasses
    {
        const int inactiveFlag{0};
        const int virtualNodeFlag{1};
        const int innerNodeFlag{2};
    }

    /// @brief A class implementing some of the CutCell functionality
    class CutCell
    {

    public:
        /// @brief CutCell ctor
        /// @brief mesh a shared pointer to mesh
        explicit CutCell(std::shared_ptr<Mesh2D> mesh);

        /// @brief Gets the edges crossed by a segment
        /// @param[in] boundaryLines The boundary polyline
        /// @return A vector containing the polyline vertices
        [[nodiscard]] std::vector<int> ClassifyNodes(const std::vector<Point>& boundaryLines) const;

    private:

        // A reference to the mesh
        std::shared_ptr<Mesh2D> m_mesh; ///< Pointer to mesh
        
    };
} // namespace meshkernel
