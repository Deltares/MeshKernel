//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2020.
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
#include <MeshKernel/CurvilinearParametersNative.hpp>

namespace meshkernel
{
    class CurvilinearGrid;
    class Splines;

    class CurvilinearGridFromSplinesTransfinite
    {
    public:
        /// <summary>
        /// Ctor with splines and parameters
        /// </summary>
        /// <returns></returns>
        CurvilinearGridFromSplinesTransfinite(std::shared_ptr<Splines> splines, meshkernelapi::CurvilinearParametersNative curvilinearParametersNative);

        /// @brief Computes the adimensional intersections between splines.
        /// Also orders the m splines (the horizontal ones) before the n splines (the vertical ones)
        void ComputeIntersections();

        /// Computes the curvilinear grid from the splines using transfinite interpolation
        /// @param curvilinearGrid
        void Compute(CurvilinearGrid& curvilinearGrid);

        std::shared_ptr<Splines> m_splines; // A pointer to spline

    private:
        /// @brief Order the splines such that their index increases in m or n direction
        /// @param startFirst
        /// @param endFirst
        /// @param startSecond
        /// @param endSecond
        /// @returns Boolean to indicate that procedure has to be repeated
        [[nodiscard]] bool OrderSplines(int startFirst,
                                        int endFirst,
                                        int startSecond,
                                        int endSecond);

        /// @brief Swap the rows of a two dimensional vector
        /// @param v The input vector
        /// @param firstRow The first row
        /// @param secondRow The second row
        template <typename T>
        void SwapRows(std::vector<std::vector<T>>& v, int firstRow, int secondRow) const;

        /// @brief Swap the columns of a two dimensional vector (MAKESR)
        /// @tparam T The input vector
        /// @param v
        /// @param firstColumn The first column
        /// @param secondColumn The second column
        template <typename T>
        void SwapColumns(std::vector<std::vector<T>>& v, int firstColumn, int secondColumn) const;

        /// Compute the distances following an exponential increase
        /// @param[in] factor
        /// @param[in] leftDistance
        /// @param[in] rightDistance
        /// @param[out] distances
        void ComputeExponentialDistances(double factor,
                                         double leftDistance,
                                         double rightDistance,
                                         std::vector<double>& distances) const;

        /// @brief Computes the distances along the spline where to generate the points
        /// @param[in] numIntersections
        /// @param[in] numPoints
        /// @param[in] numDiscretizations
        /// @param[in] intersectionDistances
        /// @param[out] distances
        void ComputeDiscretizations(int numIntersections,
                                    int numPoints,
                                    int numDiscretizations,
                                    const std::vector<double>& intersectionDistances,
                                    std::vector<double>& distances) const;

        std::vector<int> m_splineType;                                          // The spline types (1 horizontal, -1 vertical)
        std::vector<std::vector<double>> m_splineIntersectionRatios;            // For each spline, stores the intersections in terms of total spline length
        std::vector<std::vector<int>> m_splineGroupIndexAndFromToIntersections; // For each spline: position in m or n group, from and to spline crossing indices (MN12)
        int m_numMSplines = -1;                                                 // The index of the last m spline
        int m_numNSplines = -1;                                                 // The index of the last m spline
        int m_numM = 0;                                                         // Number of m columns
        int m_numN = 0;                                                         // Number of n rows
    };

} // namespace meshkernel
