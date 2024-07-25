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
        /// @brief Generate the curvilinear grid from the set of splines
        CurvilinearGrid Compute(const Splines& splines,
                                const CurvilinearParameters& curvilinearParameters) const;

    private:
        /// @brief Maximum number of spline points allowed when doubling of the spline points.
        static constexpr UInt MaximumNumberOfSplinePoints = 1000;

        /// @brief The maximum number of checks for unlabeled splines
        static constexpr UInt MaximumCumulativeUnlabeledSplineCount = 1000;

        /// @brief The maximum refinement factor.
        static constexpr int MaximumRefinementFactor = 1000;

        /// @brief Increase factor when adding mode spling support points.
        ///
        /// If value other than 2 is used, then the function IncreaseSplinePoints
        /// needs to be updated to reflect this change
        static constexpr UInt SplineIncreaseFactor = 2;

        /// @brief Array of doubles
        using DoubleVector = std::vector<double>;

        using VectorOfDoubleVectors = std::vector<DoubleVector>;

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

        /// @brief Bounded array of integers.
        using ArrayOfThree = std::array<int, 3>;

        /// @brief Vector of ArrayOfThree.
        using VectorOfThreeInts = std::vector<ArrayOfThree>;

        /// @brief Increase the number of spline support points for all splines.
        void IncreaseSplinePoints(Splines& splines) const;

        /// @brief Compute the intersection of all splines and order into splines along m- and n-grid-lines
        void ComputeSplineIntersections(Splines& splines,
                                        VectorOfDoubleVectors& splineIntersections,
                                        UInt& numMSplines) const;

        /// @brief Generate the grid points
        void GenerateGrid(const Splines& splines,
                          const VectorOfDoubleVectors& splineIntersections,
                          const VectorOfThreeInts& splineInteraction,
                          const UInt numMSplines,
                          const UInt mRefinement,
                          const UInt nRefinement,
                          CurvilinearGrid& grid) const;

        /// Generate points for a m- or n-grid-line
        void GenerateGridPoints(const Splines& splines,
                                const UInt whichSpline,
                                const UInt refinementFactor,
                                DoubleVector& intersectionPoints,
                                std::vector<Point>& gridPoints) const;

        /// @brief Assign the points along all splines.
        void GenerateGridPointsAlongSpline(const Splines& splines,
                                           const VectorOfDoubleVectors& splineIntersections,
                                           const VectorOfThreeInts& splineInteraction,
                                           const UInt numMSplines,
                                           const UInt mRefinement,
                                           const UInt nRefinement,
                                           lin_alg::Matrix<Point>& gridNodes) const;

        /// @brief Assign the interior points to a all patches.
        ///
        /// A patch is bounded by the intersection of four splines
        void FillPatchesWithPoints(const Splines& splines,
                                   const UInt numMSplines,
                                   const UInt mRefinement,
                                   const UInt nRefinement,
                                   lin_alg::Matrix<Point>& gridNodes) const;

        /// @brief Assign the interior points to a patch.
        ///
        /// A patch is bounded by the intersection of four splines
        void AssignPatchGridPoints(const UInt i, const UInt j, const UInt nRefinement, const UInt mRefinement,
                                   const lin_alg::Matrix<Point>& patchNodes,
                                   lin_alg::Matrix<Point>& gridNodes) const;

        /// @brief Ensure that the normalised distances are increasing with array index
        void PrepareNormalisedDistances(const UInt intervalStart,
                                        const UInt intervalEnd,
                                        DoubleVector& normalisedDistances) const;

        // makessq
        void ComputeDiscretisedDistances(const DoubleVector& splineIntervalLengths,
                                         const UInt refinementFactor,
                                         DoubleVector& discretisedDistances) const;

        void ComputeExponentialDistances(const double factor,
                                         const double intervalStart,
                                         const double intervalEnd,
                                         DoubleVector& intervalDiscretisation) const;

        void ComputeSplineIntervalLength(const Splines& splines,
                                         const UInt whichSpline,
                                         double& normalisedDistance,
                                         double& intervalLength) const;

        /// @brief Check all spline are valid.
        ///
        /// Ensure that there are no splines with only a single point
        bool AreSplinesValid(const Splines& splines) const;

        /// @brief Order the splines and spline-intersections into m- and n-grid-lines
        void OrderSplines(Splines& splines,
                          const UInt numMSplines,
                          VectorOfDoubleVectors& splineIntersections) const;

        /// @brief Sort a section of the splines, either along m- or n-grid-lines
        bool SortSplines(Splines& splines,
                         const UInt outerStartIndex,
                         const UInt outerEndIndex,
                         const UInt innerStartIndex,
                         const UInt innerEndIndex,
                         VectorOfDoubleVectors& splineIntersections,
                         bool& splinesSwapped) const;

        /// @brief Compute a point at a distance along the spline
        Point ComputePoint(const Splines& splines,
                           const UInt whichSpline,
                           const DoubleVector& intersectionPoints,
                           const double ssq) const;

        /// @brief Get only the non-zero values from a row of the intersection matrix (paktij)
        DoubleVector CompressRow(const VectorOfDoubleVectors& splineIntersections,
                                 const UInt whichRow) const;

        /// @brief Determine the intersection point or points of two splines
        void DetermineIntersection(const Splines& splines,
                                   const UInt splineI,
                                   const UInt splineJ,
                                   UInt& numberTimesCrossing,
                                   double& crossProductOfIntersection,
                                   double& firstNormalisedIntersectionLength,
                                   double& secondNormalisedIntersectionLength) const;

        /// @brief Compute intersections of all splines and update the spline-type if necessary.
        bool ComputeInteractions(Splines& splines,
                                 std::vector<int>& splineType,
                                 VectorOfDoubleVectors& splineIntersections) const;

        /// @brief Compute the intersection point between two splines and determine if there is a problem with the intersection of the splines.
        ///
        /// May throw AlgorithmError if:
        /// 1. a doubling of the spline points is required but the doubling will exceed the maximum number of points allowed.
        /// 2. The splines intersect more than once.
        bool ComputeAndCheckIntersection(Splines& splines,
                                         const UInt splineI,
                                         const UInt splineJ,
                                         std::vector<int>& splineType,
                                         VectorOfDoubleVectors& splineIntersections) const;

        /// @brief One or more splines remain unlabeled.
        /// @return true if one or more splines is unlabeled, false if all splines have been labeled.
        bool SplinesRemainUnlabeled(const std::vector<int>& splineType, UInt& unlabledSplineCount) const;

        /// @brief Sort spline type and spline intersections into m- and n-grid-lines
        void SortInteractionsOnSplineType(Splines& splines,
                                          std::vector<int>& splineType,
                                          VectorOfDoubleVectors& splineIntersections) const;

        /// @brief Determine if splines are along m- or n-grid lines
        void DetermineSplineOrientation(const Splines& splines,
                                        const UInt numMSplines,
                                        const VectorOfDoubleVectors& splineIntersections,
                                        VectorOfThreeInts& splineInteraction) const;

        /// @brief Assign a single point along one or other boundary arrays
        void AssignBoundaryPoint(const UInt loopIndex,
                                 const UInt boundaryIndex,
                                 const UInt refinementFactor,
                                 std::vector<Point>& startBoundaryPoints,
                                 std::vector<Point>& endBoundaryPoints,
                                 const Point gridNode) const;
    };

} // namespace meshkernel
