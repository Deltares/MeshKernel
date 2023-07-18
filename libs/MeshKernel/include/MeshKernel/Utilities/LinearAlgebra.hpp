#pragma once

#include "MeshKernel/Definitions.hpp"

#include <Eigen/Core>

namespace meshkernel
{

    namespace lin_alg
    {
        /// @brief Row major dynamic matrix
        /// @tparam T Data type
        template <class T>
        using MatrixRowMajor = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

        /// @brief Column major dynamic matrix
        /// @tparam T Data type
        template <class T>
        using MatrixColMajor = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

        /// @brief Resizes and fills a two dimensional vector
        /// @tparam T           The type of the matrix elements
        /// @param[in,out] matrix   The matrix
        /// @param[in]     rows     The number of rows
        /// @param[in]     cols     The number of columns
        /// @param[in]     preserve Optional parameter that determines if the existing elements in
        ///                         the matrix should be presrved or overwritten by the provided
        ///                         or deflault (if not provided) fill value
        /// @param[in]     fill_value The optional the fill value (defaults to the template class constructor)
        template <class T>
        inline static void ResizeAndFillMatrix(MatrixRowMajor<T>& matrix,
                                               UInt rows,
                                               UInt cols,
                                               bool preserve = false,
                                               T const& fill_value = {})
        {
            UInt const rows_old = static_cast<UInt>(matrix.rows());
            UInt const cols_old = static_cast<UInt>(matrix.cols());
            matrix.resize(rows, cols);
            if (!preserve)
            {
                matrix.fill(fill_value);
            }
            else
            {
                for (UInt i = rows_old; i < matrix.rows(); ++i)
                {
                    for (UInt j = cols_old; j < matrix.cols(); ++j)
                    {
                        matrix(i, j) = fill_value;
                    }
                }
            }
        }
    } // namespace lin_alg

} // namespace meshkernel
