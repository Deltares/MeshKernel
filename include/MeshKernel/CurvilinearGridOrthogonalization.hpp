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

#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Splines.hpp>
#include <MeshKernelApi/OrthogonalizationParameters.hpp>

namespace meshkernel
{
    class CurvilinearGrid;

    /// @brief A class implementing the curvilinear grid orthogonalization algorithm
    class CurvilinearGridOrthogonalization
    {
    public:
        /// @brief Class constructor
        /// @param[in] grid                        The input curvilinear grid
        /// @param[in] orthogonalizationParameters The orthogonalization parameters
        CurvilinearGridOrthogonalization(std::shared_ptr<CurvilinearGrid> grid,
                                         const meshkernelapi::OrthogonalizationParameters& orthogonalizationParameters);

        /// @brief Orthogonalize the curvilinear grid (modifies the grid point by m_grid)
        void Compute();

        /// @brief Sets the orthogonalization block
        /// @param[in] firstCornerPoint            The first point defining the orthogonalization bounding box
        /// @param[in] secondCornerPoint           The second point defining the orthogonalization bounding box
        void SetBlock(Point const& firstCornerPoint, Point const& secondCornerPoint);

        /// @brief Sets a line in the grid that will not move during the orthogonalization process
        /// @param[in] firstPoint The geometry list containing the first point of the line to freeze
        /// @param[in] secondPoint The geometry list containing the second point of the line to freeze
        void SetLineToFreeze(Point const& firstPoint, Point const& secondPoint);

    private:
        /// @brief Solve one orthogonalization iteration, using the method of successive over-relaxation SOR (ORTSOR)
        void Solve();

        /// @brief Project the m boundary nodes onto the original grid (BNDSMT)
        void ProjectHorizontalBoundaryGridNodes();

        /// @brief Project the n boundary nodes onto the original grid (BNDSMT)
        void ProjectVerticalBoundariesGridNodes();

        /// @brief Account for frozen lines (FIXDDBOUNDARIES)
        void AccountForFrozenLines() const;

        /// @brief Computes the matrix coefficients (ATPPAR)
        void ComputeCoefficients();

        /// @brief Computes the matrix coefficients for m-gridlines (SOMDIST)
        void ComputeHorizontalCoefficients();

        /// @brief Compute the matrix coefficients for n-gridlines (SOMDIST)
        void ComputeVerticalCoefficients();

        /// @brief Orders the m and n coordinates of the points on the grid
        /// @param firstPointM The m coordinate of the first point
        /// @param firstPointN The n coordinate of the first point
        /// @param secondPointM The m coordinate of the second point
        /// @param secondPointN The n coordinate of the second point
        /// @return The min m, min n, max m and max n of the input points
        std::tuple<size_t, size_t, size_t, size_t> OrderCoordinates(size_t firstPointM, size_t firstPointN, size_t secondPointM, size_t secondPointN);

        /// @brief Some nodes on m boundary grid lines
        [[nodiscard]] std::vector<std::vector<bool>> ComputeInvalidHorizontalBoundaryNodes() const;

        /// @brief Some nodes on n boundary grid lines
        [[nodiscard]] std::vector<std::vector<bool>> ComputeInvalidVerticalBoundaryNodes() const;

        std::shared_ptr<CurvilinearGrid> m_grid;                                  ///< A pointer to the curvilinear grid to modify
        meshkernelapi::OrthogonalizationParameters m_orthogonalizationParameters; ///< The orthogonalization parameters

        size_t m_minM; ///< The minimum m grid index of the orthogonalization bounding box
        size_t m_minN; ///< The minimum n grid index of the orthogonalization bounding box
        size_t m_maxM; ///< The maximum m grid index of the orthogonalization bounding box
        size_t m_maxN; ///< The maximum n grid index of the orthogonalization bounding box

        std::vector<std::vector<double>> m_a;   ///< The a term of the orthogonalization equation
        std::vector<std::vector<double>> m_b;   ///< The b term of the orthogonalization equation
        std::vector<std::vector<double>> m_c;   ///< The c term of the orthogonalization equation
        std::vector<std::vector<double>> m_d;   ///< The d term of the orthogonalization equation
        std::vector<std::vector<double>> m_e;   ///< The e term of the orthogonalization equation
        std::vector<std::vector<double>> m_atp; ///< The atp term of the orthogonalization equation

        std::vector<std::tuple<size_t, size_t, size_t, size_t>> m_frozenLines; ///< The frozen lines, expresse as start point and end point m and n corrdinates

        Splines m_splines; ///< The grid lines stored as splines
    };
} // namespace meshkernel
