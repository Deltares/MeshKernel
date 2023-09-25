#pragma once

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

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

    /// @brief A dynamically sized column vector
    /// @tparam T Data type
    template <class T>
    using ColumnVector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

} // namespace lin_alg
