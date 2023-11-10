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

#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Utilities/LinearAlgebra.hpp"

namespace meshkernel
{

    /// @brief Array containing dimensions of the hessian
    using HessianDimension = std::array<UInt, 3>;

    /// @brief Define column major orientation
    using MatrixColMajor = lin_alg::Matrix<double, Eigen::ColMajor>;

    /// @brief The hessian values
    ///
    /// Implemented as an array of matrices.
    /// Not sure what is the best implementation yet for performance.
    class Hessian
    {
    public:
        /// @brief Default constructor
        Hessian() = default;

        /// @brief Constructor taking 3 parameters
        Hessian(const UInt dim1, const UInt dim2, const UInt dim3);

        /// @brief Resize taking 3 parameters
        void resize(const UInt dim1, const UInt dim2, const UInt dim3);

        /// @brief Get the dimension for each dimension
        ///
        /// @param [in] dim For which dimension is the size required, dim in range [0,2]
        UInt size(const UInt dim) const;

        /// @brief Get all the Hessian dimensions
        const HessianDimension& size() const;

        /// @brief Get the 1-dimension index of
        UInt get1DIndex(const UInt dim2, const UInt dim3) const;

        /// @brief Get the value of the hessian
        double operator()(const UInt dim1, const UInt dim2, const UInt dim3) const;

        /// @brief Get the value of the hessian
        double& operator()(const UInt dim1, const UInt dim2, const UInt dim3);

        /// @brief Access the matrix in 'dim1' as though it were a 1 dimensional array
        ///
        /// dim2 = i + size(1) * j?
        double operator()(const UInt dim1, const UInt dim2) const;

        /// @brief Get the matrix for a dimension
        const MatrixColMajor& getMatrix(const UInt dim) const;

        /// @brief Get the matrix for a dimension
        MatrixColMajor& getMatrix(const UInt dim);

        /// @brief Set all entries to zero.
        void zero();

    private:
        // Since the size of the first index will be 2, and most accesses are vary the first index fastest
        // try: std::vector<lin_alg::MatrixColMajor<std::array<double,2>>> m_hessian
        std::vector<MatrixColMajor> m_hessian;
        HessianDimension m_dimensions{0, 0, 0};
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
    // return m_dimensions[2] * dim2 + dim3;
    return dim2 + m_dimensions[1] * dim3;
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

inline const meshkernel::MatrixColMajor& meshkernel::Hessian::getMatrix(const UInt dim) const
{
    return m_hessian[dim];
}

inline meshkernel::MatrixColMajor& meshkernel::Hessian::getMatrix(const UInt dim)
{
    return m_hessian[dim];
}
