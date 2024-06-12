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

#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"
#include "MeshKernel/Parameters.hpp"
#include "MeshKernel/Splines.hpp"

namespace meshkernel
{

    // TODO rename, there exists a class CurvilinearGridGridFromSplines
    // Or rename the other
    class CurvilinearGridSplineToGrid
    {
    public:
        void Compute(const Splines& splines,
                     const CurvilinearParameters& curvilinearParameters,
                     CurvilinearGrid& grid) const;

        // CurvilinearGrid Compute(const Splines& splines,
        //                         const CurvilinearParameters& curvilinearParameters) const;

    private:
        /// @brief Maximum number of spline points allowed when doubling of the spline points.
        static const UInt MaximumNumberOfSplinePoints = 10; // TODO make higher for production code

        /// @brief The maximum number of
        static const UInt MaximumCumulativeUnlabeledSplineCount = 1000;

        //
#undef USE_EIGEN

#ifdef USE_EIGEN
        template <typename T, int storage = Eigen::RowMajor>
        using EigenMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, storage>;
#else
        template <typename Type>
        using DoubleVector = std::vector<Type>;

        template <typename Type>
        using EigenMatrix = std::vector<DoubleVector<Type>>;

        template <typename T>
        void SwapColumns(std::vector<std::vector<T>>& v, UInt firstColumn, UInt secondColumn) const
        {
            for (UInt i = 0; i < v.size(); i++)
            {
                if (firstColumn >= v[i].size() || secondColumn >= v[i].size())
                {
                    continue;
                }

                std::swap(v[i][firstColumn], v[i][secondColumn]);
            }
        }

#endif

        using AnotherMatrix = std::vector<std::array<int, 3>>;

        UInt longestSplineLength(const Splines& splines) const;

        /// @brief Double the number of spline support points for all splines.
        void DoubleSplinePoints(Splines& splines) const;

        void ComputeSplineIntersections(Splines& splines,
                                        EigenMatrix<double>& splineIntersections,
                                        UInt& numMSplines) const;

        void splrgf(Splines& splines,
                    const EigenMatrix<double>& splineIntersections,
                    const AnotherMatrix& mn12,
                    CurvilinearGrid& grid,
                    const UInt numMSplines,
                    const UInt mFac,
                    const UInt nFac) const;

        void GenerateGridPoints(const Splines& splines,
                                const UInt whichSpline,
                                const UInt mFac,
                                std::vector<double>& intersectionPoints,
                                std::vector<Point>& gridPoints) const;

        void makessq(const std::vector<double>& fixedPoints,
                     const UInt mFac,
                     std::vector<double>& ssq) const;

        void makesr(const double ar,
                    const double s0,
                    const double s1,
                    std::vector<double>& sr) const;

        void getdis(const Splines& splines,
                    const UInt whichSpline,
                    double& tValue,
                    double& sValue) const;

        bool CheckSplines(const Splines& splines) const;

        void OrderSplines(Splines& splines,
                          const UInt numMSplines,
                          EigenMatrix<double>& splineIntersections) const;

        bool SortSplines(Splines& splines,
                         const UInt outerStartIndex,
                         const UInt outerEndIndex,
                         const UInt innerStartIndex,
                         const UInt innerEndIndex,
                         EigenMatrix<double>& splineIntersections,
                         bool& jaChange) const;

        Point GetXy(const Splines& splines,
                    const UInt whichSpline,
                    const std::vector<double>& intersectionPoints,
                    const double ssq) const;

        // great name (paktij)
        std::vector<double> CompressRow(const EigenMatrix<double>& splineIntersections,
                                        const UInt whichRow) const;

        void determineIntersection(Splines& splines,
                                   const UInt splineI,
                                   const UInt splineJ,
                                   UInt& numberTimesCrossing,
                                   double& crossProductOfIntersection,
                                   double& firstNormalisedIntersectionLength,
                                   double& secondNormalisedIntersectionLength) const;

        bool ComputeInteractions(Splines& splines,
                                 std::vector<int>& splineType,
                                 EigenMatrix<double>& splineIntersections) const;

        /// @brief One or more splines remain unlabeled.
        /// @return true if one or more splines is unlabeled, false if all splines have been labeled.
        bool SplinesRemainUnlabeled(const std::vector<int>& splineType, UInt& unlabledSplineCount) const;

        void SortInteractionsOnSplineType(Splines& splines,
                                          std::vector<int>& splineType,
                                          EigenMatrix<double>& splineIntersections) const;

        UInt GetNumberOfSplinesInDirectionM(const std::vector<int>& splineType) const;

        void ComputeSplineInteractions(const Splines& splines,
                                       const UInt numMSplines,
                                       const EigenMatrix<double>& splineIntersections,
                                       AnotherMatrix& splineInteraction) const;

        void assignBoundaryPoint(const UInt loopIndex,
                                 const UInt boundaryIndex,
                                 const UInt mnFac,
                                 std::vector<Point>& startBoundaryPoints,
                                 std::vector<Point>& endBoundaryPoints,
                                 const Point gridNode) const;
    };

} // namespace meshkernel
