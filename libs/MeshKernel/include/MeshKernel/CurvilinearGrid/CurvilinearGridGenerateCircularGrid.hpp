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

#include <memory>

#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Parameters.hpp"
#include "MeshKernel/Point.hpp"
#include "MeshKernel/UndoActions/UndoAction.hpp"

namespace meshkernel
{

    class CurvilinearGridGenerateCircularGrid final
    {
    public:
        /// @brief Generate the circular grid
        CurvilinearGrid Compute(const MakeGridParameters& parameters, const Projection projection) const;

    private:
        ///@{
        /// @name Generation of graded rectangular grid

        /// @brief Compute the x/m-values
        std::vector<double> ComputeXValues(const MakeGridParameters& parameters) const;

        /// @brief Compute the y/n-values
        std::vector<double> ComputeYValues(const MakeGridParameters& parameters) const;

        /// @brief Generate the graded rectangular grid
        lin_alg::Matrix<Point> GenerateRectangularGrid(const MakeGridParameters& parameters) const;
        ///@}

        ///@{
        /// @name Generation of radial grid

        /// @brief Compute the radius-values
        std::vector<double> ComputeRadiusValues(const MakeGridParameters& parameters) const;

        /// @brief Compute the theta-values (angle values)
        std::vector<double> ComputeThetaValues(const MakeGridParameters& parameters) const;

        /// @brief Generate the circular grid
        lin_alg::Matrix<Point> GenerateCircularGrid(const MakeGridParameters& parameters) const;
        ///@}

        /// @brief Generate the grid points for the curvilinear grid.
        lin_alg::Matrix<Point> GenerateGridPoints(const MakeGridParameters& parameters) const;
    };

} // namespace meshkernel
