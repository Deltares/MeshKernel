#pragma once

#include <Eigen/Core>

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
} // namespace lin_alg
