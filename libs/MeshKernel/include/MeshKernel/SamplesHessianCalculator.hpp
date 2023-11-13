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

#include <vector>

#include "MeshKernel/Hessian.hpp"
#include "MeshKernel/Point.hpp"
#include "MeshKernel/Utilities/LinearAlgebra.hpp"

namespace meshkernel
{
    /// @brief A class implementing the computation of the real component the local hessian eigenvalues
    class SamplesHessianCalculator final
    {
    public:
        /// @brief Computation of the samples Hessian
        /// @param[in] rawSamplePoints The raw sample points
        /// @param[in] projection The mesh projection
        /// @param[in] numberOfSmoothingIterations The number of smoothing iterations for the samples
        /// @param[in] numX The number of samples x coordinates
        /// @param[in] numY The number of samples y coordinates
        /// @returns A vector of samples containing the real component the local hessian eigenvalues
        static std::vector<Sample> ComputeSamplesHessian(const std::vector<Sample>& rawSamplePoints,
                                                         const Projection projection,
                                                         UInt numberOfSmoothingIterations,
                                                         const UInt numX,
                                                         const UInt numY);

    private:
        /// @brief Smooth sample data
        ///
        /// From (smooth_samples.f90)
        /// @param[in] sampleData The sample points
        /// @param[in] numberOfSmoothingIterations  The number of smoothing iterations for the samples
        /// @param[in] hessian The hessian matrix where to store the result of the hessian calculations
        static void SmoothSamples(const std::vector<Sample>& sampleData,
                                  const UInt numberOfSmoothingIterations,
                                  Hessian& hessian);

        /// @brief Compute the gradient in a control volume defined by the polygon (0-R-1-L, comp_grad.f90)
        static void ComputeGradient(const std::vector<Sample>& samplePoints,
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
                                    double& dareaR);

        /// @brief Compute the gradient of the sample data along an edge (comp_samplegradi.f90)
        static void ComputeSampleGradient(const std::vector<Sample>& samplePoints,
                                          const Projection projection,
                                          const Hessian& hessian,
                                          const UInt direction,
                                          const UInt i,
                                          const UInt j,
                                          meshkernel::Vector& gradient,
                                          meshkernel::Vector& sn,
                                          double& dareaL,
                                          double& dareaR);

        /// @brief Compute the Hessian From (comp_samplehessian.f90)
        static void ComputeHessian(const std::vector<Sample>& samplePoints,
                                   const Projection projection,
                                   Hessian& hessian);

        /// @brief Prepare the sample data for the Hessian computation (prepare_samplehessian.f90)
        static void PrepareSampleForHessian(const std::vector<Sample>& samplePoints,
                                            const Projection projection,
                                            UInt numberOfSmoothingIterations,
                                            Hessian& hessian);
    };

} // namespace meshkernel
