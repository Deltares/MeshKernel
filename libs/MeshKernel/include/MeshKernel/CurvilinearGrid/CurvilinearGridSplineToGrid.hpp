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
#include "MeshKernel/Splines.hpp"

namespace meshkernel
{

    // TODO rename, there exists a class CurvilinearGridGridFromSplines
    // Or rename the other
    class CurvilinearGridSplineToGrid
    {
    public:
        void Compute(const Splines& splines, CurvilinearGrid& grid) const;

    private:
        using DoubleVector = std::vector<double>;
        using DoubleMatrix = std::vector<DoubleVector>;

        UInt longestSplineLength(const Splines& splines) const;

        void sectr(Splines& splines,
                   lin_alg::Matrix<double>& splineIntersections,
                   lin_alg::Matrix<int>& mn12) const;

        void splrgf(Splines& splines,
                    const lin_alg::Matrix<double>& splineIntersections,
                    const lin_alg::Matrix<int>& mn12,
                    CurvilinearGrid& grid) const;

        void makespl(const Splines& splines,
                     const UInt whichSpline,
                     const UInt mFac,
                     const std::vector<double>& intersectionPoints) const;

        void makessq(const std::vector<double>& fixedPoints,
                     const UInt mFac,
                     std::vector<double>& ssq) const;

        void makesr (const double ar,
                     const double s0,
                     const double s1,
                     std::vector<double>& sr) const;

        void getdis(const Splines& splines,
                    const UInt whichSpline,
                    double& tValue,
                    double& sValue) const;

        bool checkSplines(const Splines& splines) const;

        // great name
        std::vector<double> paktij(const lin_alg::Matrix<double>& splineIntersections,
                                   const UInt whichRow) const;

        void determineIntersection(Splines& splines,
                                   const UInt i,
                                   const UInt j,
                                   UInt& numberTimesCrossing,
                                   double& crossProductOfIntersection,
                                   double& firstNormalisedIntersectionLength,
                                   double& secondNormalisedIntersectionLength) const;
    };

} // namespace meshkernel
