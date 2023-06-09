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

#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Parameters.hpp>

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
        /// @param[in] projection The projection to use
        CurvilinearGridCreateUniform(Projection projection);

        /// @brief Compute an uniform curvilinear grid using the make mesh parameters
        /// @param[in] makeGridParameters The parameters to use for the curvilinear grid generation
        CurvilinearGrid Compute(const MakeGridParameters& makeGridParameters) const;

        /// @brief Compute an uniform curvilinear grid in one polygon
        /// @param[in] makeGridParameters The parameters to use for the curvilinear grid generation
        /// @param[in] polygons The input polygons
        /// @param[in] polygonIndex The polygon index
        CurvilinearGrid Compute(const MakeGridParameters& makeGridParameters,
                                std::shared_ptr<Polygons> polygons,
                                size_t polygonIndex) const;

    private:
        /// @brief Compute an uniform curvilinear grid on spherical coordinates.
        /// A correction is applied to the longitudinal discretization to preserve an aspect ratio ds/dy = 1 on real distances
        /// For preventing small edges, the correction stops when the distance is less than 2000 meters and the grid is generated around the poles
        /// This is a customized fix for GTSM models.
        /// @param[in] makeGridParameters The parameters to use for the curvilinear grid generation
        static std::vector<std::vector<Point>> computeSpherical(const MakeGridParameters& makeGridParameters);

        /// @brief Compute an uniform curvilinear grid on cartesian coordinates.
        /// @param[in] makeGridParameters The parameters to use for the curvilinear grid generation
        static std::vector<std::vector<Point>> computeCartesian(const MakeGridParameters& makeGridParameters);

        Projection m_projection; ///< The projection to use
    };
} // namespace meshkernel
