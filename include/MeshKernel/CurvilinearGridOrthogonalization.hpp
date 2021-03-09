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

    /// @brief A class implementing the curvilinear grid derefinement algorithm
    class CurvilinearGridOrthogonalization
    {
    public:
        /// @brief Class constructor
        /// @param[in] grid The input curvilinear grid
        /// @param[in] firstPoint The first vertex of the segment defining the derefinement zone
        /// @param[in] secondPoint The second vertex of the segment defining the refinement zone
        CurvilinearGridOrthogonalization(std::shared_ptr<CurvilinearGrid> grid,
                                         const meshkernelapi::OrthogonalizationParameters& orthogonalizationParameters,
                                         const Point& firstPoint,
                                         const Point& secondPoint);

        /// @brief Refine the curvilinear grid
        void Compute();

    private:
        /// @brief Matrix solve (BNDSMT)
        void TreatBoundaryConditions();

        /// @brief Matrix solve (ORTSOR)
        void Solve();

        /// @brief Computes matrix coefficients (FIXDDBOUNDARIES)
        void FreezeBoundaries();

        /// @brief Computes matrix coefficients (ATPPAR)
        void ComputeHorizontalMatrixCoefficients();

        /// @brief Computes the horizontal matrix coefficients (SOMDIST)
        void ComputeHorizontalCoefficients();

        /// @brief Computes the vertical matrix coefficients (SOMDIST)
        void ComputeVerticalCoefficients();

        std::shared_ptr<CurvilinearGrid> m_grid;                                  ///< A pointer to the curvilinear grid to modify
        meshkernelapi::OrthogonalizationParameters m_orthogonalizationParameters; ///< The orthogonalization parameters
        Point m_firstPoint;                                                       ///< The first vertex of the segment defining the derefinement zone
        Point m_secondPoint;                                                      ///< The second vertex of the segment defining the derefinement zone

        // the node indices
        size_t m_minM;
        size_t m_minN;
        size_t m_maxM;
        size_t m_maxN;

        /// matrix coefficients
        std::vector<std::vector<double>> m_a;
        std::vector<std::vector<double>> m_b;
        std::vector<std::vector<double>> m_c;
        std::vector<std::vector<double>> m_d;
        std::vector<std::vector<double>> m_e;
        std::vector<std::vector<double>> m_atp;

        Splines m_splines;
    };
} // namespace meshkernel