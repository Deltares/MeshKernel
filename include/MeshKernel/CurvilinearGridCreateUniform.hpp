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

#include <MeshKernelApi/MakeMeshParameters.hpp>

#include <memory>

namespace meshkernel
{
    class Polygons;
    class CurvilinearGrid;

    /// @brief A class implementing the generation of a uniform curvilinear grid
    class CurvilinearGridCreateUniform
    {
    public:
        /// @brief Class constructor
        ///
        /// @param[in] MakeMeshParameters The structure containing the make grid parameters
        /// @param[in] polygons The polygon to account for
        CurvilinearGridCreateUniform(const meshkernelapi::MakeMeshParameters& MakeMeshParameters, std::shared_ptr<Polygons> polygons);

        /// @brief Compute an uniform curvilinear grid
        CurvilinearGrid Compute() const;

    private:
        meshkernelapi::MakeMeshParameters m_makeMeshParameters; ///< A copy of the structure containing the parameters used for making the grid
        std::shared_ptr<Polygons> m_polygons;                   ///< A pointer to the polygon to use
    };
} // namespace meshkernel
