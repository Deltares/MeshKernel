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

#include "MeshKernel/Point.hpp"
#include "MeshKernel/Utilities/LinearAlgebra.hpp"
#include "MeshKernel/Utilities/RTree.hpp"

// namespace lin_alg
// {
//     /// @brief Row major dynamic matrix
//     /// @tparam T Data type
//     template <class T>
//     using MatrixRowMajor = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

//     /// @brief Column major dynamic matrix
//     /// @tparam T Data type
//     template <class T>
//     using MatrixColMajor = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

namespace meshkernel
{

    using HessianDimension = std::array<UInt, 3>;

    class Hessian
    {
    public:
        Hessian() = default;

        Hessian(const UInt dim1, const UInt dim2, const UInt dim3);

        void resize(const UInt dim1, const UInt dim2, const UInt dim3);

        /// @brief Get the dimension for each dimension?
        ///
        /// @param [in] dim For which dimension is the size required, dim in range [0,2]
        UInt size(const UInt dim) const;

        const HessianDimension& size() const;

        /// @brief Get the 1-dimension index of
        UInt get1DIndex(const UInt dim2, const UInt dim3) const;

        double operator()(const UInt dim1, const UInt dim2, const UInt dim3) const;

        double& operator()(const UInt dim1, const UInt dim2, const UInt dim3);

        /// @brief Access the matrix in 'dim1' as though it were a 1 dimensional array
        ///
        /// dim2 = i + size(1) * j?
        double operator()(const UInt dim1, const UInt dim2) const;

        const lin_alg::MatrixColMajor<double>& getMatrix(const UInt dim) const;

        lin_alg::MatrixColMajor<double>& getMatrix(const UInt dim);

        /// @brief Set all entries to zero.
        void zero();

    private:
        // Since the size of the first index will be 5, and most accesses are vary the first index fastest
        // try: std::vector<lin_alg::MatrixColMajor<std::array<double,5>>> m_hessian
        std::vector<lin_alg::MatrixColMajor<double>> m_hessian;
        HessianDimension m_dimensions{0, 0, 0};
    };

    class RidgeRefinement final
    {
    public:
        void Compute(const std::vector<Point>& rawSamplePoints,
                     const std::vector<double>& rawSampleData,
                     const Projection projection,
                     const UInt numX, const UInt numY) const;

        // void Compute(Mesh2D& mesh,
        //              const std::vector<Point>& rawSamplePoints,
        //              const std::vector<double>& rawSampleData,
        //              const Projection projection,
        //              const UInt numX, const UInt numY) const;

    private:
        void GetValidPoints(const std::vector<Point>& samplePoints, std::vector<Point>& validSamplePoints, std::vector<UInt>& iperm) const;

        void smoothSamples(const std::vector<double>& sampleData,
                           const UInt numberOfSmoothingIterations,
                           Hessian& hessian) const;

        void computeGradient(const std::vector<Point>& samplePoints,
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

        void computeSampleGradient(const std::vector<Point>& samplePoints,
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

        void computeHessian(const std::vector<Point>& samplePoints,
                            const std::vector<double>& sampleData,
                            const Projection projection,
                            Hessian& hessian) const;

        /// @brief Remove in-valid and duplicate sample points and data.
        void RemoveDuplicates(std::vector<Point>& samplePoints, std::vector<double>& sampleData, std::vector<UInt>& sampleIndices /*ipsam*/) const;

        void prepareSampleForHessian(const std::vector<Point>& samplePoints,
                                     const std::vector<double>& sampleData,
                                     const Projection projection,
                                     Hessian& hessian) const;

        void TidySamples(std::vector<Point>& samplePoints, std::vector<double>& sampleData) const;
    };

} // namespace meshkernel

inline meshkernel::UInt meshkernel::Hessian::size(const UInt dim) const
{
    return m_dimensions[dim];
}

inline const meshkernel::HessianDimension& meshkernel::Hessian::size() const
{
    return m_dimensions;
}

inline meshkernel::UInt meshkernel::Hessian::get1DIndex(const UInt dim2, const UInt dim3) const
{
    return dim2 + size(1) * dim3;
}

inline double meshkernel::Hessian::operator()(const UInt dim1, const UInt dim2, const UInt dim3) const
{
    return m_hessian[dim1](dim2, dim3);
}

inline double& meshkernel::Hessian::operator()(const UInt dim1, const UInt dim2, const UInt dim3)
{
    return m_hessian[dim1](dim2, dim3);
}

inline double meshkernel::Hessian::operator()(const UInt dim1, const UInt dim2) const
{
    return m_hessian[dim1](dim2);
}

inline const lin_alg::MatrixColMajor<double>& meshkernel::Hessian::getMatrix(const UInt dim) const
{
    return m_hessian[dim];
}

inline lin_alg::MatrixColMajor<double>& meshkernel::Hessian::getMatrix(const UInt dim)
{
    return m_hessian[dim];
}
