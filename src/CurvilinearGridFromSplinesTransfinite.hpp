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
#include "Entities.hpp"
#include "CurvilinearParametersNative.hpp"
#include "SplinesToCurvilinearParametersNative.hpp"

namespace GridGeom
{
    class CurvilinearGrid;
    class Splines;

    class CurvilinearGridFromSplinesTransfinite
    {
    public:

        /// <summary>
        /// Ctor
        /// </summary>
        /// <returns></returns>
        CurvilinearGridFromSplinesTransfinite();

        /// <summary>
        /// Ctor
        /// </summary>
        /// <returns></returns>
        CurvilinearGridFromSplinesTransfinite(Splines* splines);

        /// <summary>
        /// Computes the adimensional intersections between splines.
        /// Also orders the n splines (the horizontal ones) before the m splines (the vertical ones)
        /// </summary>
        /// <returns></returns>
        bool ComputeSplineIntersections();

        /// <summary>
        /// Computes the curvilinear grid from the splines
        /// </summary>
        /// <param name="curvilinearGrid"></param>
        /// <returns></returns>
        bool Compute(CurvilinearGrid& curvilinearGrid);
        
        Splines* m_splines;                                      // A pointer to spline

    private:

        bool OrderSplines(int startFirst,
            int endFirst,
            int startSecond,
            int endSecond);

        template<typename T>
        bool SwapRows(std::vector<std::vector<T>>& v, int firstRow, int secondRow);

        template<typename T>
        bool SwapColumns(std::vector<std::vector<T>>& v, int firstColumn, int secondColumn);


        std::vector<int>                 m_splineType;
        std::vector<std::vector<double>> m_splineIntersectionRatios;
        std::vector<std::vector<int>>    m_countSplineIntersections;  // spline index, 0 n splines intersections, 1 m spline intersections 
        int                              m_firstNSplines = -1;
        int                              m_numM = 0;
        int                              m_numN = 0;

    };

}


