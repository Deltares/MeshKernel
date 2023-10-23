//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2023.
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

#include <array>
#include <vector>

#include "MeshKernel/Hessian.hpp"
#include "MeshKernel/Point.hpp"
#include "MeshKernel/Utilities/LinearAlgebra.hpp"
#include "MeshKernel/Utilities/RTree.hpp"

namespace meshkernel
{

#if 0

    class RidgeRefinement final
    {
    public:
        /// @brief Intended to be the computation of the Hessian
        void Compute(const std::vector<Point>& rawSamplePoints,
                     const std::vector<double>& rawSampleData,
                     const Projection projection,
                     const UInt numX, const UInt numY) const;

        /// @brief Tidy the sample points.
        ///
        /// Remove invlalid points and duplicate points.
        /// From (tidysample.f90)
        void TidySamples(std::vector<Point>& samplePoints, std::vector<double>& sampleData) const;

    private:
        /// @brief Remove invalid points
        void GetValidPoints(const std::vector<Point>& samplePoints, std::vector<Point>& validSamplePoints, std::vector<UInt>& iperm) const;

        /// @brief Smooth sample data
        ///
        /// From (smooth_samples.f90)
        void SmoothSamples(const std::vector<double>& sampleData,
                           const UInt numberOfSmoothingIterations,
                           Hessian& hessian) const;

        /// @brief Compute the gradient in a control volume defined by the polygon (0-R-1-L)
        ///
        /// From (comp_grad.f90)
        void ComputeGradient(const std::vector<Point>& samplePoints,
                             const std::vector<double>& sampleData,
                             const Projection projection,
                             const Hessian& hessian,
                             const UInt ip0,
                             const UInt ip1,
                             const UInt ip0L,
                             const UInt ip0R,
                             const UInt ip1L,
                             const UInt ip1R,
                             meshkernel::Vector& gradient,
                             meshkernel::Vector& S,
                             double& dareaL,
                             double& dareaR) const;

        /// @brief Compute the gradient of the sample data along an edge
        ///
        /// From (comp_samplegradi.f90)
        void ComputeSampleGradient(const std::vector<Point>& samplePoints,
                                   const std::vector<double>& sampleData,
                                   const Projection projection,
                                   const Hessian& hessian,
                                   const UInt direction,
                                   const UInt i,
                                   const UInt j,
                                   meshkernel::Vector& gradient,
                                   meshkernel::Vector& sn,
                                   double& dareaL,
                                   double& dareaR) const;

        /// @brief Compute the Hessian
        ///
        /// From (comp_samplehessian.f90)
        void ComputeHessian(const std::vector<Point>& samplePoints,
                            const std::vector<double>& sampleData,
                            const Projection projection,
                            Hessian& hessian) const;

        /// @brief Remove in-valid and duplicate sample points and data.
        void RemoveDuplicates(std::vector<Point>& samplePoints, std::vector<double>& sampleData, std::vector<UInt>& sampleIndices /*ipsam*/) const;

        /// @brief Prepare the sample data for the Hessian computation.
        ///
        /// From (prepare_samplehessian.f90)
        void PrepareSampleForHessian(const std::vector<Point>& samplePoints,
                                     const std::vector<double>& sampleData,
                                     const Projection projection,
                                     Hessian& hessian) const;
    };
#endif

#if 0
    class HessianCalculator final
    {
    public:
        /// @brief Intended to be the computation of the Hessian
        void Compute(const std::vector<Sample>& rawSamplePoints,
                     const Projection projection,
                     const UInt numX,
                     const UInt numY,
                     Hessian& hessian) const;

        /// @brief Intended to be the computation of the Hessian
        void Compute(const std::vector<Sample>& rawSamplePoints,
                     const Projection projection,
                     const UInt numX,
                     const UInt numY,
                     std::vector<Sample>& hessianSamples) const;

    private:
        /// @brief Smooth sample data
        ///
        /// From (smooth_samples.f90)
        void SmoothSamples(const std::vector<Sample>& sampleData,
                           const UInt numberOfSmoothingIterations,
                           Hessian& hessian) const;

        /// @brief Compute the gradient in a control volume defined by the polygon (0-R-1-L)
        ///
        /// From (comp_grad.f90)
        void ComputeGradient(const std::vector<Sample>& samplePoints,
                             const Projection projection,
                             const Hessian& hessian,
                             const UInt ip0,
                             const UInt ip1,
                             const UInt ip0L,
                             const UInt ip0R,
                             const UInt ip1L,
                             const UInt ip1R,
                             meshkernel::Vector& gradient,
                             meshkernel::Vector& S,
                             double& dareaL,
                             double& dareaR) const;

        /// @brief Compute the gradient of the sample data along an edge
        ///
        /// From (comp_samplegradi.f90)
        void ComputeSampleGradient(const std::vector<Sample>& samplePoints,
                                   const Projection projection,
                                   const Hessian& hessian,
                                   const UInt direction,
                                   const UInt i,
                                   const UInt j,
                                   meshkernel::Vector& gradient,
                                   meshkernel::Vector& sn,
                                   double& dareaL,
                                   double& dareaR) const;

        /// @brief Compute the Hessian
        ///
        /// From (comp_samplehessian.f90)
        void ComputeHessian(const std::vector<Sample>& samplePoints,
                            const Projection projection,
                            Hessian& hessian) const;

        /// @brief Prepare the sample data for the Hessian computation.
        ///
        /// From (prepare_samplehessian.f90)
        void PrepareSampleForHessian(const std::vector<Sample>& samplePoints,
                                     const Projection projection,
                                     Hessian& hessian) const;
    };
#endif

    class HessianSampleCalculator final
    {
    public:
        /// @brief Intended to be the computation of the Hessian
        static void Compute(const std::vector<Sample>& rawSamplePoints,
                            const Projection projection,
                            const UInt numX,
                            const UInt numY,
                            std::vector<Sample>& hessian);

    private:
        static meshkernel::UInt get1DIndex(const UInt numX, const UInt numY, const UInt dim1, const UInt dim2)
        {
            (void)numY;
            (void)numX;
            // return numY * dim1 + dim2;
            return dim1 + numX * dim2;
        }

        /// @brief Smooth sample data
        ///
        /// From (smooth_samples.f90)
        static void SmoothSamples(const UInt numX,
                                  const UInt numY,
                                  const UInt numberOfSmoothingIterations,
                                  std::vector<Sample>& sampleData);

        /// @brief Compute the gradient in a control volume defined by the polygon (0-R-1-L)
        ///
        /// From (comp_grad.f90)
        static void ComputeGradient(const std::vector<Sample>& sampleData,
                                    const Projection projection,
                                    const UInt ip0,
                                    const UInt ip1,
                                    const UInt ip0L,
                                    const UInt ip0R,
                                    const UInt ip1L,
                                    const UInt ip1R,
                                    meshkernel::Vector& gradient,
                                    meshkernel::Vector& S,
                                    double& dareaL,
                                    double& dareaR);

        /// @brief Compute the gradient of the sample data along an edge
        ///
        /// From (comp_samplegradi.f90)
        static void ComputeSampleGradient(const std::vector<Sample>& samplePoints,
                                          const UInt numX,
                                          const UInt numY,
                                          const Projection projection,
                                          const UInt direction,
                                          const UInt i,
                                          const UInt j,
                                          meshkernel::Vector& gradient,
                                          meshkernel::Vector& sn,
                                          double& dareaL,
                                          double& dareaR);

        /// @brief Compute the Hessian
        ///
        /// From (comp_samplehessian.f90)
        static void ComputeHessian(const std::vector<Sample>& samplePoints,
                                   const UInt numX,
                                   const UInt numY,
                                   const Projection projection,
                                   std::vector<Sample>& hessian);
    };

} // namespace meshkernel
