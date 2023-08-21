#pragma once

#include <Eigen/Core>

#include <stdexcept>
#include <vector>

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

    /// @brief  Dynamic vector
    /// @tparam T Data type
    template <typename T>
    using Vector = Eigen::Matrix<T, 1, Eigen::Dynamic, Eigen::RowMajor>;

    /// @brief     Checks if a matrix is empty by inspecting its size
    /// @tparam    T      The type of the matrix elements
    /// @param[in] matrix The matrix
    /// @return    True if the matrix is empty, false otherwise
    template <typename T>
    [[nodiscard]] inline static bool MatrixIsEmpty(MatrixRowMajor<T> const& matrix)
    {
        return matrix.size() == 0;
    }

    /// @brief         Resizes and fills a matrix
    /// @tparam        T        The data type of the matrix elements
    /// @param[in,out] matrix   The matrix
    /// @param[in]     rows     The number of rows
    /// @param[in]     cols     The number of columns
    /// @param[in]     preserve Optional parameter that determines if the existing elements in
    ///                         the matrix should be preserved or overwritten by the provided
    ///                         or default (if not provided) fill value
    /// @param[in]     fill_value The optional fill value (defaults to the template class constructor)
    /// @return        True if the input matrix has been resized, false otherwise
    template <typename T>
    inline static bool ResizeAndFillMatrix(MatrixRowMajor<T>& matrix,
                                           Eigen::Index rows,
                                           Eigen::Index cols,
                                           bool preserve = false,
                                           T const& fill_value = {})
    {
        if (rows < 0 || cols < 0)
        {
            throw std::invalid_argument("Invalid range");
        }

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

    /// @brief         Erases a range of rows from a matrix given the begin and end row indices
    /// @tparam        T         The type of the matrix elements
    /// @param[in,out] matrix    The matrix
    /// @param[in]     col_begin The begin row index
    /// @param[in]     col_end   The end row index
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

    /// @brief        Erases a row from a matrix given the row index
    /// @tparam        T      The type of the matrix elements
    /// @param[in,out] matrix The matrix
    /// @param[in]     row    The row index
    template <typename T>
    inline static void EraseRow(MatrixRowMajor<T>& matrix, Eigen::Index row)
    {
        EraseRows(matrix, row, row);
    }

    /// @brief         Erases a range of columns from a matrix given the begin and end column indices
    /// @tparam        T         The type of the matrix elements
    /// @param[in,out] matrix    The matrix
    /// @param[in]     col_begin The begin column index
    /// @param[in]     col_end   The end column index
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

    /// @brief         Erases a column from a matrix given the column index
    /// @tparam        T      The type of the matrix elements
    /// @param[in,out] matrix The matrix
    /// @param[in]     col    The column index
    template <typename T>
    inline static void EraseCol(MatrixRowMajor<T>& matrix, Eigen::Index col)
    {
        EraseCols(matrix, col, col);
    }

    /// @brief         Inserts an Eigen vector in a matrix at a given row index
    /// @tparam        T      The type of the matrix elements
    /// @param[in,out] matrix The matrix
    /// @param[in]     vector The vector to insert
    /// @param[in]     row    The row index where the vector is to be inserted in the matrix
    template <typename T>
    inline static void InsertRow(MatrixRowMajor<T>& matrix,
                                 Vector<T> const& vector,
                                 Eigen::Index row)
    {
        if (row < 0 || row > matrix.rows())
        {
            throw std::invalid_argument("Invalid range");
        }

        Eigen::Index const rows_new = matrix.rows() + 1;
        Eigen::Index const rows_to_move_down = matrix.rows() - row;
        Eigen::Index const cols = matrix.cols();

        auto const block_to_move_down = matrix.block(row, 0, rows_to_move_down, cols).eval();
        matrix.conservativeResize(rows_new, cols);
        matrix.block(row, 0, 1, cols) = vector;
        matrix.block(row + 1, 0, rows_to_move_down, cols) = block_to_move_down;
    }

    /// @brief     Transforms a STL vector to an Eigen row-major vector
    /// @tparam    T      The type of the matrix elements
    /// @param[in] vector The STL vector
    /// @return    The Eigen row-major vector
    template <typename T>
    inline static Vector<T> STLVectorToVector(std::vector<T> const& stl_vector)
    {
        // Vector<T> vector(1, stl_vector.size());
        // for (Eigen::Index i = 0; i < vector.cols(); ++i)
        // {
        //     vector(0, i) = stl_vector[i];
        // }
        // return vector;

        return Eigen::Map<Vector<T> const>(stl_vector.data(), 1, stl_vector.size());
    }

    /// @brief         Inserts a STL vector in a matrix at a given row index
    /// @tparam        T      The type of the matrix elements
    /// @param[in,out] matrix     The matrix
    /// @param[in]     stl_vector The STL vector to insert
    /// @param[in]     row    The row index where the vector is to be inserted in the matrix
    template <typename T>
    inline static void InsertRow(MatrixRowMajor<T>& matrix,
                                 std::vector<T> const& stl_vector,
                                 Eigen::Index row)
    {
        InsertRow(matrix, STLVectorToVector(stl_vector), row);
    }

    /// @brief         Inserts an Eigen vector in a matrix at a given column index
    /// @tparam        T      The type of the matrix elements
    /// @param[in,out] matrix The matrix
    /// @param[in]     vector The vector to insert
    /// @param[in]     col    The column index where the vector is to be inserted in the matrix
    template <typename T>
    inline static void InsertCol(MatrixRowMajor<T>& matrix,
                                 Vector<T> const& vector,
                                 Eigen::Index col)
    {
        if (col < 0 || col > matrix.cols())
        {
            throw std::invalid_argument("Invalid range");
        }

        Eigen::Index const rows = matrix.rows();
        Eigen::Index const cols_new = matrix.cols() + 1;
        Eigen::Index const cols_to_move_right = matrix.cols() - col;

        auto const block_to_move_right = matrix.block(0, col, rows, cols_to_move_right).eval();
        matrix.conservativeResize(rows, cols_new);
        matrix.block(0, col, rows, 1) = vector.transpose();
        matrix.block(0, col + 1, rows, cols_to_move_right) = block_to_move_right;
    }

    /// @brief         Inserts a STL vector in a matrix at a given column index
    /// @tparam        T      The type of the matrix elements
    /// @param[in,out] matrix     The matrix
    /// @param[in]     stl_vector The STL vector to insert
    /// @param[in]     col    The column index where the vector is to be inserted in the matrix
    template <typename T>
    inline static void InsertCol(MatrixRowMajor<T>& matrix,
                                 std::vector<T> const& stl_vector,
                                 Eigen::Index col)
    {
        InsertCol(matrix, STLVectorToVector(stl_vector), col);
    }

} // namespace lin_alg
