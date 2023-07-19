#pragma once

#include <Eigen/Core>

#include <stdexcept>

namespace meshkernel
{

    namespace lin_alg
    {
        /// @brief Row major dynamic matrix
        /// @tparam T Data type
        template <typename T>
        using MatrixRowMajor = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

        /// @brief Column major dynamic matrix
        /// @tparam T Data type
        template <class T>
        using MatrixColMajor = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

        // @brief Dynamic matrix
        /// @tparam T Data type
        template <typename T, Eigen::StorageOptions storage_option>
        using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, storage_option>;

        /// @brief Checks if a matrix is empty by inspecting its size
        /// @tparam T      The type of the matrix elements
        /// @param matrix  The matrix
        /// @return        true if the matrix is empty, false otherwise
        template <class T>
        [[nodiscard]] inline static bool MatrixIsEmpty(MatrixRowMajor<T> const& matrix)
        {
            return matrix.size() == 0;
        }

        /// @brief Resizes and fills a two dimensional vector
        /// @tparam T               The type of the matrix elements
        /// @param[in,out] matrix   The matrix
        /// @param[in]     rows     The number of rows
        /// @param[in]     cols     The number of columns
        /// @param[in]     preserve Optional parameter that determines if the existing elements in
        ///                         the matrix should be preserved or overwritten by the provided
        ///                         or default (if not provided) fill value
        /// @param[in]     fill_value The optional the fill value (defaults to the template class constructor)
        template <typename T>
        inline static bool ResizeAndFillMatrix(MatrixRowMajor<T>& matrix,
                                               Eigen::Index rows,
                                               Eigen::Index cols,
                                               bool preserve = false,
                                               T const& fill_value = {})
        {
            Eigen::Index const rows_old = matrix.rows();
            Eigen::Index const cols_old = matrix.cols();
            bool const resize = (rows != rows_old) || (cols != cols_old);
            if (!resize)
            {
                if (!preserve)
                {
                    matrix.fill(fill_value);
                }
            }
            else
            {
                if (!preserve)
                {
                    matrix.resize(rows, cols);
                    matrix.fill(fill_value);
                }
                else
                {
                    matrix.conservativeResize(rows, cols);
                    if (rows > rows_old)
                    {
                        matrix.bottomRows(rows - rows_old).setConstant(fill_value);
                    }
                    if (cols > cols_old)
                    {
                        matrix.topRightCorner(rows_old, cols - cols_old).setConstant(fill_value);
                    }
                }
            }
            return resize;
        }

        /// @brief Erases a range of rows from a matrix given the begin and end row indices
        /// @tparam T The type of the matrix elements
        /// @param matrix The matrix
        /// @param col_begin The begin row index
        /// @param col_end The end row index
        template <typename T>
        inline static void EraseRows(MatrixRowMajor<T>& matrix,
                                     Eigen::Index row_begin,
                                     Eigen::Index row_end)
        {
            if (row_begin < 0 ||
                row_end < row_begin ||
                row_end > matrix.rows() - 1)
            {
                throw std::invalid_argument("Invalid range");
            }

            Eigen::Index const rows_to_remove = row_end - row_begin + 1;
            Eigen::Index const rows_new = matrix.rows() - rows_to_remove;
            Eigen::Index const rows_to_move_up = rows_new - row_begin;
            Eigen::Index const cols = matrix.cols();

            // move the block of rows to keep to the position of the block of rows to remove
            matrix.block(row_begin, 0, rows_to_move_up, cols) =
                matrix.block(row_end + 1, 0, rows_to_move_up, cols);

            // resize and preserve matrix contents
            matrix.conservativeResize(rows_new, cols);
        }

        /// @brief Erases a row from a matrix given the row index
        /// @tparam T The type of the matrix elements
        /// @param[in,out] matrix The matrix
        /// @param[in] row The row index
        template <typename T>
        inline static void EraseRow(MatrixRowMajor<T>& matrix, Eigen::Index row)
        {
            EraseRows(matrix, row, row);
        }

        /// @brief Erases a range of columns from a matrix given the begin and end column indices
        /// @tparam T The type of the matrix elements
        /// @param matrix The matrix
        /// @param col_begin The begin column index
        /// @param col_end The end column index
        template <typename T>
        inline static void EraseCols(MatrixRowMajor<T>& matrix,
                                     Eigen::Index col_begin,
                                     Eigen::Index col_end)
        {
            if (col_begin < 0 ||
                col_end < col_begin ||
                col_end > matrix.cols() - 1)
            {
                throw std::invalid_argument("Invalid range");
            }

            Eigen::Index const rows = matrix.rows();
            Eigen::Index const cols_to_remove = col_end - col_begin + 1;
            Eigen::Index const cols_new = matrix.cols() - cols_to_remove;
            Eigen::Index const cols_to_move_left = cols_new - col_begin;

            // move the block of columns to keep to the position of the block of columns to remove
            matrix.block(0, col_begin, rows, cols_to_move_left) =
                matrix.block(0, col_end + 1, rows, cols_to_move_left);

            // resize and preserve matrix contents
            matrix.conservativeResize(rows, cols_new);
        }

        /// @brief Erases a column from a matrix given the column index
        /// @tparam T The type of the matrix elements
        /// @param[in,out] matrix The matrix
        /// @param[in] row The column index
        template <typename T>
        inline static void EraseCol(MatrixRowMajor<T>& matrix, Eigen::Index col)
        {
            EraseCols(matrix, col, col);
        }

    } // namespace lin_alg

} // namespace meshkernel
