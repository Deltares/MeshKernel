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

namespace MeshKernel
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
        /// Ctor with splines and parameters
        /// </summary>
        /// <returns></returns>
        CurvilinearGridFromSplinesTransfinite(std::shared_ptr<Splines> splines, MeshKernelApi::CurvilinearParametersNative curvilinearParametersNative);


        /// <summary>
        /// Computes the adimensional intersections between splines.
        /// Also orders the m splines (the horizontal ones) before the n splines (the vertical ones)
        /// </summary>
        /// <returns>If the method succeeded</returns>
        bool ComputeIntersections();

        /// <summary>
        /// Computes the curvilinear grid from the splines using transfinite interpolation
        /// </summary>
        /// <param name="curvilinearGrid"></param>
        /// <returns>If the method succeeded</returns>
        bool Compute(CurvilinearGrid& curvilinearGrid);

        
        std::shared_ptr<Splines> m_splines;                                      // A pointer to spline

    private:

        /// <summary>
        /// Order the splines such that their index increases in m or n direction
        /// </summary>
        /// <param name="startFirst"></param>
        /// <param name="endFirst"></param>
        /// <param name="startSecond"></param>
        /// <param name="endSecond"></param>
        /// <returns>If the method succeeded</returns>
        bool OrderSplines(int startFirst,
                          int endFirst,
                          int startSecond,
                          int endSecond);

        /// <summary>
        /// Swap the rows of a two dimensional vector
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="v">The input vector</param>
        /// <param name="firstRow">The first row</param>
        /// <param name="secondRow">The second row</param>
        /// <returnsIf the method succeeded></returns>
        template<typename T>
        bool SwapRows(std::vector<std::vector<T>>& v, int firstRow, int secondRow);

        /// <summary>
        /// Swap the columns of a two dimensional vector (MAKESR)
        /// </summary>
        /// <typeparam name="T">The input vector</typeparam>
        /// <param name="v"></param>
        /// <param name="firstColumn">The first column</param>
        /// <param name="secondColumn">The second column</param>
        /// <returns>If the method succeeded</returns>
        template<typename T>
        bool SwapColumns(std::vector<std::vector<T>>& v, int firstColumn, int secondColumn);


        /// <summary>
        /// Compute the distances following an exponential increase
        /// </summary>
        /// <param name="factor"></param>
        /// <param name="leftDistance"></param>
        /// <param name="rightDistance"></param>
        /// <param name="distances"></param>
        /// <returns></returns>
        bool ComputeExponentialDistances(double factor, 
                                           double leftDistance,
                                           double rightDistance,
                                           std::vector<double>& distances) const;

        /// <summary>
        /// Computes the distances along the spline where to generate the points
        /// </summary>
        /// <param name="numIntersections"></param>
        /// <param name="numPoints"></param>
        /// <param name="numDiscretizations"></param>
        /// <param name="intersectionDistances"></param>
        /// <param name="distances"></param>
        /// <returns></returns>
        bool ComputeDiscretizations( int numIntersections,
                                     int numPoints,
                                     int numDiscretizations,
                                     const std::vector<double>& intersectionDistances,
                                     std::vector<double>& distances ) const;


        std::vector<int>                 m_splineType;                              // The spline types (1 horizontal, -1 vertical)
        std::vector<std::vector<double>> m_splineIntersectionRatios;                // For each spline, stores the intersections in terms of total spline length
        std::vector<std::vector<int>>    m_splineGroupIndexAndFromToIntersections;  // For each spline: position in m or n group, from and to spline crossing indexses (MN12)
        int                              m_numMSplines = -1;                        // The index of the last m spline
        int                              m_numNSplines = -1;                        // The index of the last m spline
        int                              m_numM = 0;                                // Number of m columns
        int                              m_numN = 0;                                // Number of n rows

    };

}


