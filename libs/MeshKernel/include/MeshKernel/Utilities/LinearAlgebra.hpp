#pragma once

#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Exceptions.hpp"

#include <Eigen/Core>

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace lin_alg
{

    /// @brief  Dynamic matrix
    /// @tparam T       Matrix data type
    /// @tparam storage Matrix storage option
    template <typename T, int storage = Eigen::RowMajor>
    using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, storage>;

    template <typename T, int storage = Eigen::RowMajor>
    using MatrixRow = Eigen::Block<Matrix<T, storage>, 1, Eigen::Dynamic, true>;

    /// @brief  Dynamic column vector
    /// @tparam T Data type
    template <typename T>
    using ColVector = Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::ColMajor>;

    /// @brief  Dynamic row vector
    /// @tparam T Data type
    template <typename T>
    using RowVector = Eigen::Matrix<T, 1, Eigen::Dynamic, Eigen::RowMajor>;

    /// @brief     Checks if a matrix is empty by inspecting its size
    /// @tparam    T       Matrix data type
    /// @tparam    storage Matrix storage option
    /// @param[in] matrix  The matrix
    /// @return    True if the matrix is empty, false otherwise
    template <typename T, int storage>
    [[nodiscard]] inline static bool MatrixIsEmpty(Matrix<T, storage> const& matrix)
    {
        return matrix.size() == 0;
    }

    /// @brief         Resizes and fills a matrix
    /// @tparam        T        The data type of the matrix elements
    /// @tparam        storage  Matrix storage option
    /// @param[in,out] matrix   The matrix
    /// @param[in]     rows     The number of rows
    /// @param[in]     cols     The number of columns
    /// @param[in]     preserve Optional parameter that determines if the existing elements in
    ///                         the matrix should be preserved or overwritten by the provided
    ///                         or default (if not provided) fill value
    /// @param[in]     fill_value The optional fill value (defaults to the template class constructor)
    /// @return        True if the input matrix has been resized, false otherwise
    template <typename T, int storage>
    inline static bool ResizeAndFillMatrix(Matrix<T, storage>& matrix,
                                           Eigen::Index rows,
                                           Eigen::Index cols,
                                           bool preserve = false,
                                           T const& fill_value = {})
    {
        if (rows < 0 || cols < 0)
        {
            throw meshkernel::LinearAlgebraError("Invalid dimensions: rows = {}, cols = {}. Dimensions must be strictly positive.",
                                                 rows,
                                                 cols);
        }

        Eigen::Index const rows_old = matrix.rows();
        Eigen::Index const cols_old = matrix.cols();

        bool const resize = (rows != rows_old) || (cols != cols_old);
        if (!resize && !preserve)
        {
            matrix.fill(fill_value);
        }

        if (resize && !preserve)
        {
            matrix.resize(rows, cols);
            matrix.fill(fill_value);
        }

        if (resize && preserve)
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

        return resize;
    }

    /// @brief         Erases a range of rows from a matrix given the begin and end row indices
    /// @tparam        T         Matrix data type
    /// @tparam        storage   Matrix storage option
    /// @param[in,out] matrix    The matrix
    /// @param[in]     row_begin The begin row index
    /// @param[in]     row_end   The end row index
    template <typename T, int storage>
    inline static void EraseRows(Matrix<T, storage>& matrix,
                                 Eigen::Index row_begin,
                                 Eigen::Index row_end)
    {
        if (row_begin < 0 ||
            row_end < row_begin ||
            row_end > matrix.rows() - 1)
        {
            throw meshkernel::LinearAlgebraError("Invalid range: row_begin_index = {}, row_end_index = {}, max_row_index = {}",
                                                 row_begin,
                                                 row_end,
                                                 matrix.rows() - 1);
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
    /// @tparam        T       Matrix data type
    /// @tparam        storage Matrix storage option
    /// @param[in,out] matrix The matrix
    /// @param[in]     row    The row index
    template <typename T, int storage>
    inline static void EraseRow(Matrix<T, storage>& matrix, Eigen::Index row)
    {
        EraseRows(matrix, row, row);
    }

    /// @brief         Erases a range of columns from a matrix given the begin and end column indices
    /// @tparam        T         Matrix data type
    /// @tparam        storage   Matrix storage option
    /// @param[in,out] matrix    The matrix
    /// @param[in]     col_begin The begin column index
    /// @param[in]     col_end   The end column index
    template <typename T, int storage>
    inline static void EraseCols(Matrix<T, storage>& matrix,
                                 Eigen::Index col_begin,
                                 Eigen::Index col_end)
    {
        if (col_begin < 0 ||
            col_end < col_begin ||
            col_end > matrix.cols() - 1)
        {
            throw meshkernel::LinearAlgebraError("Invalid range: col_begin_index = {}, col_end_index = {}, max_col_index = {}",
                                                 col_begin,
                                                 col_end,
                                                 matrix.cols() - 1);
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
    /// @tparam        T       Matrix data type
    /// @tparam        storage Matrix storage option
    /// @param[in,out] matrix  The matrix
    /// @param[in]     col     The column index
    template <typename T, int storage>
    inline static void EraseCol(Matrix<T, storage>& matrix, Eigen::Index col)
    {
        EraseCols(matrix, col, col);
    }

    /// @brief         Inserts an Eigen vector in a matrix at a given row index
    /// @tparam        T       Matrix data type
    /// @tparam        storage Matrix storage option
    /// @param[in,out] matrix  The matrix
    /// @param[in]     vector  The vector to insert
    /// @param[in]     row     The row index where the vector is to be inserted in the matrix
    template <typename T, int storage>
    inline static void InsertRow(Matrix<T, storage>& matrix,
                                 RowVector<T> const& vector,
                                 Eigen::Index row)
    {
        if (row < 0 || row > matrix.rows())
        {
            throw meshkernel::LinearAlgebraError("Invalid range: cannot insert at index = {}.", row);
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
    /// @tparam    T       Matrix data type
    /// @tparam    storage Matrix storage option
    /// @param[in] vector  The STL vector
    /// @return    The Eigen row-major vector
    template <typename T>
    inline static RowVector<T> STLVectorToRowVector(std::vector<T> const& stl_vector)
    {
        return Eigen::Map<RowVector<T> const>(stl_vector.data(), 1, stl_vector.size());
    }

    /// @brief         Inserts a STL vector in a matrix at a given row index
    /// @tparam        T          Matrix data type
    /// @tparam        storage    Matrix storage option
    /// @param[in,out] matrix     The matrix
    /// @param[in]     stl_vector The STL vector to insert
    /// @param[in]     row        The row index where the vector is to be inserted in the matrix
    template <typename T, int storage>
    inline static void InsertRow(Matrix<T, storage>& matrix,
                                 std::vector<T> const& stl_vector,
                                 Eigen::Index row)
    {
        InsertRow(matrix, STLVectorToRowVector(stl_vector), row);
    }

    /// @brief         Inserts an Eigen vector in a matrix at a given column index
    /// @tparam        T       Matrix data type
    /// @tparam        storage Matrix storage option
    /// @param[in,out] matrix  The matrix
    /// @param[in]     vector  The vector to insert
    /// @param[in]     col     The column index where the vector is to be inserted in the matrix
    template <typename T, int storage>
    inline static void InsertCol(Matrix<T, storage>& matrix,
                                 ColVector<T> const& vector,
                                 Eigen::Index col)
    {
        if (col < 0 || col > matrix.cols())
        {
            throw meshkernel::LinearAlgebraError("Invalid range: cannot insert at index = {}.", col);
        }

        Eigen::Index const rows = matrix.rows();
        Eigen::Index const cols_new = matrix.cols() + 1;
        Eigen::Index const cols_to_move_right = matrix.cols() - col;

        auto const block_to_move_right = matrix.block(0, col, rows, cols_to_move_right).eval();
        matrix.conservativeResize(rows, cols_new);
        matrix.block(0, col, rows, 1) = vector;
        matrix.block(0, col + 1, rows, cols_to_move_right) = block_to_move_right;
    }

    /// @brief     Transforms a STL vector to an Eigen column-major vector
    /// @tparam    T       Matrix data type
    /// @tparam    storage Matrix storage option
    /// @param[in] vector  The STL vector
    /// @return    The Eigen row-major vector
    template <typename T>
    inline static ColVector<T> STLVectorToColVector(std::vector<T> const& stl_vector)
    {
        return Eigen::Map<ColVector<T> const>(stl_vector.data(), stl_vector.size(), 1);
    }

    /// @brief         Inserts a STL vector in a matrix at a given column index
    /// @tparam        T          Matrix data type
    /// @tparam        storage    Matrix storage option
    /// @param[in,out] matrix     The matrix
    /// @param[in]     stl_vector The STL vector to insert
    /// @param[in]     col        The column index where the vector is to be inserted in the matrix
    template <typename T, int storage>
    inline static void InsertCol(Matrix<T, storage>& matrix,
                                 std::vector<T> const& stl_vector,
                                 Eigen::Index col)
    {
        InsertCol(matrix, STLVectorToColVector(stl_vector), col);
    }

    /// @brief         Swaps two columns of a matrix
    /// @tparam        T       Matrix data type
    /// @tparam        storage Matrix storage option
    /// @param[in,out] matrix The matrix
    /// @param[in]     col_1 Index of first column to swap
    /// @param[in]     col_2 Index of second column to swap
    template <typename T, int storage>
    void SwapColumns(Matrix<T, storage>& matrix,
                     Eigen::Index col_1,
                     Eigen::Index col_2)
    {
        if (col_1 < 0 ||
            col_1 > matrix.cols() - 1 ||
            col_2 < 0 ||
            col_2 > matrix.cols() - 1)
        {
            throw meshkernel::LinearAlgebraError("Invalid range");
        }

        // do nothing if column indices are identical
        if (col_1 == col_2)
        {
            return;
        }

        // swap
        matrix.col(col_1).swap(matrix.col(col_2));
    }

    /// @brief         Swaps two rows of a matrix
    /// @tparam        T       Matrix data type
    /// @tparam        storage Matrix storage option
    /// @param[in,out] matrix The matrix
    /// @param[in]     row_1 Index of first row to swap
    /// @param[in]     row_2 Index of second row to swap
    template <typename T, int storage>
    void SwapRows(Matrix<T, storage>& matrix,
                  Eigen::Index row_1,
                  Eigen::Index row_2)
    {
        if (row_1 < 0 ||
            row_1 > matrix.rows() - 1 ||
            row_2 < 0 ||
            row_2 > matrix.rows() - 1)
        {
            throw meshkernel::LinearAlgebraError("Invalid range");
        }

        // do nothing if row indices are identical
        if (row_1 == row_2)
        {
            return;
        }

        // swap
        matrix.row(row_1).swap(matrix.row(row_2));
    }

    /// @brief Returns the indices of a sorted matrix row without modifying the row
    /// @tparam    T       Matrix data type
    /// @tparam    storage Matrix storage option
    /// @param[in] row The matrix row
    /// @return    The indices of the sorted matrix row
    template <typename T, int storage>
    [[nodiscard]] RowVector<Eigen::Index> SortRow(MatrixRow<T, storage> const row)
    {
        RowVector<Eigen::Index> indices(row.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::ranges::stable_sort(indices,
                                 [&row](Eigen::Index row_index_1, Eigen::Index row_index_2)
                                 { return row[row_index_1] < row[row_index_2]; });
        return indices;
    }

    /// @brief Reorders a row of a matrix according to a row vector
    /// @tparam        T       Matrix data type
    /// @tparam        storage Matrix storage option
    /// @param[in,out] row The matrix row to reorder
    /// @param[in]     order The row vector containing the reordering indices
    template <typename T, int storage>
    void ReorderRow(MatrixRow<T, storage> row,
                    RowVector<Eigen::Index> const& order)
    {
        if (order.size() != row.size())
        {
            throw meshkernel::LinearAlgebraError("The matrix row and the order vector are not of the same size.");
        }

        std::vector<bool> swapped(row.size(), false);
        for (Eigen::Index i = 0; i < row.size(); ++i)
        {
            if (swapped[i])
            {
                continue;
            }
            swapped[i] = true;
            Eigen::Index last_j = i;
            Eigen::Index j = order[i];
            while (i != j)
            {
                std::swap(row[last_j], row[j]);
                swapped[j] = true;
                last_j = j;
                j = order[j];
            }
        }
    }

    /// @brief Copies a row from a matrix to an STL vector
    /// @tparam T         Matrix data type
    /// @tparam storage   Matrix storage option
    /// @param matrix[in] The matrix
    /// @param row[in]    The index of the row to copy
    /// @return  The matrix row as an STL vector
    template <typename T, int storage>
    [[nodiscard]] inline static std::vector<T> MatrixRowToSTLVector(Matrix<T, storage> const& matrix,
                                                                    Eigen::Index row)
    {
        if (row < 0 || row > matrix.rows() - 1)
        {
            throw meshkernel::LinearAlgebraError("Invalid range");
        }

        std::vector<T> vector(matrix.cols());
        Eigen::Map<Matrix<T, storage>>(vector.data(), 1, matrix.cols()) = matrix.row(row);
        return vector;
    }

    /// @brief Copies a column from a matrix to an STL vector
    /// @tparam T         Matrix data type
    /// @tparam storage   Matrix storage option
    /// @param matrix[in] The matrix
    /// @param col[in]    The index of the column to copy
    /// @return  The matrix column as an STL vector
    template <typename T, int storage>
    [[nodiscard]] inline static std::vector<T> MatrixColToSTLVector(Matrix<T, storage> const& matrix,
                                                                    Eigen::Index col)
    {
        if (col < 0 || col > matrix.cols() - 1)
        {
            throw meshkernel::LinearAlgebraError("Invalid range");
        }

        std::vector<T> vector(matrix.rows());
        Eigen::Map<Matrix<T, storage>>(vector.data(), matrix.rows(), 1) = matrix.col(col);
        return vector;
    }

    /// @brief Copies an STL vector of vectors to a row-major matrix
    /// @tparam T Matrix data type
    /// @param vector[in] vector of vectors to copy
    /// @return Row-major matrix
    template <typename T>
    [[nodiscard]] inline static Matrix<T> STLVectorOfVectorsToMatrix(std::vector<std::vector<T>> const& vector)
    {
        Matrix<T> matrix(vector.size(), vector.front().size());
        for (size_t i = 0; i < vector.size(); ++i)
        {
            matrix.row(i) = RowVector<T>::Map(vector[i].data(), 1, vector[i].size());
        }
        return matrix;
    }

} // namespace lin_alg
