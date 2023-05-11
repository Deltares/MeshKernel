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

#include <memory>

#include <MeshKernel/CurvilinearGrid/CurvilinearGridAlgorithm.hpp>
#include <MeshKernel/Parameters.hpp>
#include <MeshKernel/Splines.hpp>

namespace meshkernel
{
    class CurvilinearGrid;

    /// @brief A class implementing the curvilinear grid orthogonalization algorithm
    class CurvilinearGridOrthogonalization : public CurvilinearGridAlgorithm
    {
    public:
        /// @brief Class constructor
        /// @param[in] grid                        The input curvilinear grid
        /// @param[in] orthogonalizationParameters The orthogonalization parameters
        CurvilinearGridOrthogonalization(std::shared_ptr<CurvilinearGrid> grid,
                                         const OrthogonalizationParameters& orthogonalizationParameters);

        /// @brief Orthogonalize the curvilinear grid (modifies the grid point by m_grid)
        CurvilinearGrid Compute() override;

    private:
        /// @brief Solve one orthogonalization iteration, using the method of successive over-relaxation SOR (ORTSOR)
        void Solve();

        /// @brief Project the m boundary nodes onto the original grid (BNDSMT)
        void ProjectHorizontalBoundaryGridNodes();

        /// @brief Project the n boundary nodes onto the original grid (BNDSMT)
        void ProjectVerticalBoundariesGridNodes();

        /// @brief Compute frozen grid points (FIXDDBOUNDARIES)
        void ComputeFrozenGridPoints();

        /// @brief Computes the matrix coefficients (ATPPAR)
        void ComputeCoefficients();

        /// @brief Computes the matrix coefficients for m-gridlines (SOMDIST)
        void ComputeHorizontalCoefficients();

        /// @brief Compute the matrix coefficients for n-gridlines (SOMDIST)
        void ComputeVerticalCoefficients();

        /// @brief Computes the invalid nodes on m boundary grid lines
        /// @return A vector with true if the node is an invalid boundary node, false otherwise
        [[nodiscard]] std::vector<std::vector<bool>> ComputeInvalidHorizontalBoundaryNodes() const;

        /// @brief Computes the invalid nodes on n boundary grid lines
        /// @return A vector with true if the node is an invalid boundary node, false otherwise
        [[nodiscard]] std::vector<std::vector<bool>> ComputeInvalidVerticalBoundaryNodes() const;

        OrthogonalizationParameters m_orthogonalizationParameters; ///< The orthogonalization parameters

        std::vector<std::vector<double>> m_a;   ///< The a term of the orthogonalization equation
        std::vector<std::vector<double>> m_b;   ///< The b term of the orthogonalization equation
        std::vector<std::vector<double>> m_c;   ///< The c term of the orthogonalization equation
        std::vector<std::vector<double>> m_d;   ///< The d term of the orthogonalization equation
        std::vector<std::vector<double>> m_e;   ///< The e term of the orthogonalization equation
        std::vector<std::vector<double>> m_atp; ///< The atp term of the orthogonalization equation

        std::vector<std::vector<bool>> m_isGridNodeFrozen; ///< A mask for setting some of the grid nodes frozen

        Splines m_splines; ///< The grid lines stored as splines
    };
} // namespace meshkernel
