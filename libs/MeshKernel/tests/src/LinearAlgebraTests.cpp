#include <gtest/gtest.h>

#include <MeshKernel/Utilities/LinearAlgebra.hpp>

TEST(LinearAlgebra, Empty)
{
    {
        lin_alg::Matrix<int> matrix;
        EXPECT_TRUE(lin_alg::MatrixIsEmpty(matrix));
    }

    {
        lin_alg::Matrix<int> matrix(3, 4);
        EXPECT_FALSE(lin_alg::MatrixIsEmpty(matrix));
    }
}

TEST(LinearAlgebra, Resize)
{
    lin_alg::Matrix<int> matrix;

    EXPECT_THROW(lin_alg::ResizeAndFillMatrix(matrix, -1, -1), std::invalid_argument);

    {
        Eigen::Index constexpr rows = 3;
        Eigen::Index constexpr cols = 4;

        EXPECT_TRUE(lin_alg::ResizeAndFillMatrix(matrix, rows, cols)); // 3x4
        EXPECT_EQ(matrix.rows(), rows);
        EXPECT_EQ(matrix.cols(), cols);

        // fill value was not set, expect default int initialisation
        lin_alg::Matrix<int> expected_matrix(rows, cols);
        expected_matrix.setZero();
        EXPECT_EQ(matrix, expected_matrix);

        // resize with the same dimensions
        EXPECT_FALSE(lin_alg::ResizeAndFillMatrix(matrix, rows, cols));
    }

    // Add extra row and extra col and fill with ones without preserving
    {
        Eigen::Index rows_new = matrix.rows() + 1;
        Eigen::Index cols_new = matrix.cols() + 1;
        EXPECT_TRUE(lin_alg::ResizeAndFillMatrix(matrix, rows_new, cols_new, false, 1)); // 4x5
        EXPECT_EQ(matrix.rows(), rows_new);
        EXPECT_EQ(matrix.cols(), cols_new);
        lin_alg::Matrix<int> expected_matrix(rows_new, cols_new); // 4x5
        expected_matrix.setOnes();
        EXPECT_EQ(matrix, expected_matrix);
    }

    // Add extra row and extra col and fill with 2s with preserving
    {
        Eigen::Index rows_new = matrix.rows() + 1;
        Eigen::Index cols_new = matrix.cols() + 1;
        EXPECT_TRUE(lin_alg::ResizeAndFillMatrix(matrix, rows_new, cols_new, true, 2)); // 5x6
        EXPECT_EQ(matrix.rows(), rows_new);
        EXPECT_EQ(matrix.cols(), cols_new);

        lin_alg::Matrix<int> expected_matrix(rows_new, cols_new); // 5x6
        expected_matrix.setOnes();
        expected_matrix.bottomRows(1).setConstant(2);
        expected_matrix.rightCols(1).setConstant(2);
        EXPECT_EQ(matrix, expected_matrix);
    }

    // shrink by one row and one column with preserving
    {
        Eigen::Index rows_new = matrix.rows() - 1;
        Eigen::Index cols_new = matrix.cols() - 1;
        EXPECT_TRUE(lin_alg::ResizeAndFillMatrix(matrix, rows_new, cols_new, true)); // 4x5
        EXPECT_EQ(matrix.rows(), rows_new);
        EXPECT_EQ(matrix.cols(), cols_new);
        {
            lin_alg::Matrix<int> expected_matrix(rows_new, cols_new); // 4x5
            expected_matrix.setOnes();
            EXPECT_EQ(matrix, expected_matrix);
        }
    }

    // add 3 extra cols, preserve old matrix and fill 2ith 666
    {
        Eigen::Index rows_new = matrix.rows();
        Eigen::Index cols_new = matrix.cols() + 3;
        EXPECT_TRUE(lin_alg::ResizeAndFillMatrix(matrix, rows_new, cols_new, true, 666)); // 4x8
        EXPECT_EQ(matrix.rows(), rows_new);
        EXPECT_EQ(matrix.cols(), cols_new);
        {
            lin_alg::Matrix<int> expected_matrix(rows_new, cols_new); // 4x8
            expected_matrix.setOnes();
            expected_matrix.rightCols(3).setConstant(666);
            EXPECT_EQ(matrix, expected_matrix);
        }
    }

    // add 4 extra rows, preserve old matrix and fill with -666
    {
        Eigen::Index rows_new = matrix.rows() + 4;
        Eigen::Index cols_new = matrix.cols();
        EXPECT_TRUE(lin_alg::ResizeAndFillMatrix(matrix, rows_new, cols_new, true, -666)); // 8x8
        EXPECT_EQ(matrix.rows(), rows_new);
        EXPECT_EQ(matrix.cols(), cols_new);
        {
            lin_alg::Matrix<int> expected_matrix(rows_new, cols_new); // 8x5
            expected_matrix.setOnes();
            expected_matrix.rightCols(3).setConstant(666);
            expected_matrix.bottomRows(4).setConstant(-666);
            EXPECT_EQ(matrix, expected_matrix);
        }
    }
}

TEST(LinearAlgebra, EraseOneRow)
{
    lin_alg::Matrix<int> matrix(5, 3);
    matrix << 1, 2, 3,
        4, 5, 6,
        7, 8, 9,
        10, 11, 12,
        13, 14, 15;

    lin_alg::EraseRow(matrix, 2);
    EXPECT_EQ(matrix.rows(), 4);
    {
        lin_alg::Matrix<int> expected_matrix(4, 3);
        expected_matrix << 1, 2, 3,
            4, 5, 6,
            10, 11, 12,
            13, 14, 15;
        EXPECT_EQ(matrix, expected_matrix);
    }

    lin_alg::EraseRow(matrix, 1);
    EXPECT_EQ(matrix.rows(), 3);
    {
        lin_alg::Matrix<int> expected_matrix(3, 3);
        expected_matrix << 1, 2, 3,
            10, 11, 12,
            13, 14, 15;
        EXPECT_EQ(matrix, expected_matrix);
    }

    lin_alg::EraseRow(matrix, 0);
    EXPECT_EQ(matrix.rows(), 2);
    {
        lin_alg::Matrix<int> expected_matrix(2, 3);
        expected_matrix << 10, 11, 12,
            13, 14, 15;
        EXPECT_EQ(matrix, expected_matrix);
    }
}

TEST(LinearAlgebra, EraseMultipleRows)
{
    Eigen::Index constexpr cols = 3;
    lin_alg::Matrix<int> matrix(6, cols);
    matrix << 1, 2, 3,
        4, 5, 6,
        7, 8, 9,
        10, 11, 12,
        13, 14, 15,
        16, 17, 18;

    lin_alg::EraseRows(matrix, 2, 3); // remove rows 2 and 3
    EXPECT_EQ(matrix.rows(), 4);      // 6 -2 = 4 rows left
    EXPECT_EQ(matrix.cols(), cols);   // unchanged
    lin_alg::Matrix<int> expected_matrix(4, cols);
    expected_matrix << 1, 2, 3,
        4, 5, 6,
        13, 14, 15,
        16, 17, 18;
    EXPECT_EQ(matrix, expected_matrix);
}

TEST(LinearAlgebra, EraseOneCol)
{
    Eigen::Index constexpr rows = 3;
    lin_alg::Matrix<int> matrix(rows, 5);
    matrix << 1, 2, 3, 4, 5,
        6, 7, 8, 9, 10,
        11, 12, 13, 14, 15;

    lin_alg::EraseCol(matrix, 2);
    EXPECT_EQ(matrix.cols(), 4);
    {
        lin_alg::Matrix<int> expected_matrix(rows, 4);
        expected_matrix << 1, 2, 4, 5,
            6, 7, 9, 10,
            11, 12, 14, 15;
        EXPECT_EQ(matrix, expected_matrix);
    }

    lin_alg::EraseCol(matrix, 1);
    EXPECT_EQ(matrix.cols(), 3);
    {
        lin_alg::Matrix<int> expected_matrix(rows, 3);
        expected_matrix << 1, 4, 5,
            6, 9, 10,
            11, 14, 15;
        EXPECT_EQ(matrix, expected_matrix);
    }

    lin_alg::EraseCol(matrix, 0);
    EXPECT_EQ(matrix.cols(), 2);
    {
        lin_alg::Matrix<int> expected_matrix(rows, 2);
        expected_matrix << 4, 5,
            9, 10,
            14, 15;
        EXPECT_EQ(matrix, expected_matrix);
    }
}

TEST(LinearAlgebra, EraseMultipleCols)
{
    Eigen::Index constexpr rows = 3;
    lin_alg::Matrix<int> matrix(rows, 6);
    matrix << 1, 2, 3, 4, 5, 6,
        7, 8, 9, 10, 11, 12,
        13, 14, 15, 16, 17, 18;
    lin_alg::EraseCols(matrix, 2, 3); // remove cols 2 and 3
    EXPECT_EQ(matrix.rows(), rows);   // unchanged
    EXPECT_EQ(matrix.cols(), 4);      // 6 - 2 = 4 cols left
    lin_alg::Matrix<int> expected_matrix(rows, 4);
    expected_matrix << 1, 2, 5, 6,
        7, 8, 11, 12,
        13, 14, 17, 18;
    EXPECT_EQ(matrix, expected_matrix);
}

TEST(LinearAlgebra, InsertRowFromRowVector)
{
    Eigen::Index rows = 1;
    Eigen::Index constexpr cols = 3;
    lin_alg::Matrix<int> matrix(rows, cols);
    matrix << 7, 8, 9;

    // insert zeroth row
    {
        lin_alg::RowVector<int> vector(cols);
        vector << 1, 2, 3;
        // row
        // 0 : 1, 2, 3, <- insert at 0
        // 1 : 7, 8, 9
        lin_alg::InsertRow(matrix, vector, 0);
        rows++;
        EXPECT_EQ(matrix.rows(), rows); // 1+1 rows
        EXPECT_EQ(matrix.cols(), cols); // unchanged
        lin_alg::Matrix<int> expected_matrix(rows, cols);
        expected_matrix << 1, 2, 3,
            7, 8, 9;
        EXPECT_EQ(matrix, expected_matrix);
    }

    // insert at second row
    {
        lin_alg::RowVector<int> vector(cols);
        vector << 4, 5, 6;
        // row
        // 0 : 1, 2, 3,
        // 1 : 4, 5, 6, <- insert at 1
        // 2 : 7, 8, 9
        lin_alg::InsertRow(matrix, vector, 1);
        rows++;
        EXPECT_EQ(matrix.rows(), rows); // 2+1 rows
        EXPECT_EQ(matrix.cols(), cols); // unchanged
        lin_alg::Matrix<int> expected_matrix(rows, cols);
        expected_matrix << 1, 2, 3,
            4, 5, 6,
            7, 8, 9;
        EXPECT_EQ(matrix, expected_matrix);
    }

    // insert at third row
    {
        lin_alg::RowVector<int> vector(cols);
        vector << 10, 11, 12;
        // row
        // 0 : 1, 2, 3,
        // 1 : 4, 5, 6,
        // 2 : 7, 8, 9,
        // 3 : 10, 11, 12  <- insert at 3
        lin_alg::InsertRow(matrix, vector, 3);
        rows++;
        EXPECT_EQ(matrix.rows(), rows); // 3+1 rows
        EXPECT_EQ(matrix.cols(), cols); // unchanged
        lin_alg::Matrix<int> expected_matrix(rows, cols);
        expected_matrix << 1, 2, 3,
            4, 5, 6,
            7, 8, 9,
            10, 11, 12;
        EXPECT_EQ(matrix, expected_matrix);
    }

    // throw out of range exceptions
    {
        lin_alg::RowVector<int> vector(cols);
        vector << 666, 666, 666;

        // at the top side of the matrix
        // row
        // 0 : 1, 2, 3,
        // 1 : 4, 5, 6,
        // 2 : 7, 8, 9,
        // 3 : 10, 11, 12
        // 4 : can insert here but won't
        // 5 : cannot insert here, should throw
        EXPECT_THROW(lin_alg::InsertRow(matrix, vector, 5), std::invalid_argument);

        // at the bottom side of the matrix
        // row
        // -1: cannot insert here, should throw
        // 0 : 1, 2, 3,
        // 1 : 4, 5, 6,
        // 2 : 7, 8, 9,
        // 3 : 10, 11, 12
        EXPECT_THROW(lin_alg::InsertRow(matrix, vector, -1), std::invalid_argument);
    }
}

TEST(LinearAlgebra, InsertRowFromSTLVector)
{
    Eigen::Index rows = 1;
    Eigen::Index constexpr cols = 3;
    lin_alg::Matrix<int> matrix(rows, cols);
    matrix << 7, 8, 9;

    // insert zeroth row
    {
        std::vector<int> vector{1, 2, 3};
        // row
        // 0 : 1, 2, 3, <- insert at 0
        // 1 : 7, 8, 9
        lin_alg::InsertRow(matrix, vector, 0);
        rows++;
        EXPECT_EQ(matrix.rows(), rows); // 1+1 rows
        EXPECT_EQ(matrix.cols(), cols); // unchanged
        lin_alg::Matrix<int> expected_matrix(rows, cols);
        expected_matrix << 1, 2, 3,
            7, 8, 9;
        EXPECT_EQ(matrix, expected_matrix);
    }

    // insert at first row
    {
        std::vector<int> vector{4, 5, 6};
        // row
        // 0 : 1, 2, 3,
        // 1 : 4, 5, 6, <- insert at 1
        // 2 : 7, 8, 9
        lin_alg::InsertRow(matrix, vector, 1);
        rows++;
        EXPECT_EQ(matrix.rows(), rows); // 2+1 rows
        EXPECT_EQ(matrix.cols(), cols); // unchanged
        lin_alg::Matrix<int> expected_matrix(rows, cols);
        expected_matrix << 1, 2, 3,
            4, 5, 6,
            7, 8, 9;
        EXPECT_EQ(matrix, expected_matrix);
    }

    // insert at third row
    {
        std::vector<int> vector{10, 11, 12};
        // row
        // 0 : 1, 2, 3,
        // 1 : 4, 5, 6,
        // 2 : 7, 8, 9,
        // 3 : 10, 11, 12  <- insert at 3
        lin_alg::InsertRow(matrix, vector, 3);
        rows++;
        EXPECT_EQ(matrix.rows(), rows); // 3+1 rows
        EXPECT_EQ(matrix.cols(), cols); // unchanged
        lin_alg::Matrix<int> expected_matrix(rows, cols);
        expected_matrix << 1, 2, 3,
            4, 5, 6,
            7, 8, 9,
            10, 11, 12;
        EXPECT_EQ(matrix, expected_matrix);
    }

    // throw out of range exceptions
    {
        std::vector<int> vector{666, 666, 666};

        // at the bottom side of the matrix
        // row
        // 0 : 1, 2, 3,
        // 1 : 4, 5, 6,
        // 2 : 7, 8, 9,
        // 3 : 10, 11, 12
        // 4 : can insert here but won't
        // 5 : cannot insert here, should throw
        EXPECT_THROW(lin_alg::InsertRow(matrix, vector, 5), std::invalid_argument);

        // at the top side of the matrix
        // row
        // -1: cannot insert here, should throw
        // 0 : 1, 2, 3,
        // 1 : 4, 5, 6,
        // 2 : 7, 8, 9,
        // 3 : 10, 11, 12
        EXPECT_THROW(lin_alg::InsertRow(matrix, vector, -1), std::invalid_argument);
    }
}

TEST(LinearAlgebra, InsertColFromColVector)
{
    Eigen::Index constexpr rows = 3;
    Eigen::Index cols = 1;
    lin_alg::Matrix<int> matrix(rows, 1);
    matrix << 7,
        8,
        9;

    // insert zeroth column
    {
        lin_alg::ColVector<int> vector(rows);
        vector << 1, 2, 3;
        lin_alg::InsertCol(matrix, vector, 0);
        cols++;
        EXPECT_EQ(matrix.rows(), rows); // unchanged
        EXPECT_EQ(matrix.cols(), cols); // 1+1 cols
        lin_alg::Matrix<int> expected_matrix(rows, cols);
        expected_matrix << 1, 7,
            2, 8,
            3, 9;
        EXPECT_EQ(matrix, expected_matrix);
    }

    // insert first column
    {
        lin_alg::ColVector<int> vector(rows);
        vector << 4, 5, 6;
        lin_alg::InsertCol(matrix, vector, 1);
        cols++;
        EXPECT_EQ(matrix.rows(), rows); // unchanged
        EXPECT_EQ(matrix.cols(), cols); // 1+2 cols
        lin_alg::Matrix<int> expected_matrix(rows, cols);
        expected_matrix << 1, 4, 7,
            2, 5, 8,
            3, 6, 9;
        EXPECT_EQ(matrix, expected_matrix);
    }

    // insert third column
    {
        lin_alg::ColVector<int> vector(rows);
        vector << 10, 11, 12;
        lin_alg::InsertCol(matrix, vector, 3);
        cols++;
        EXPECT_EQ(matrix.rows(), rows); // unchanged
        EXPECT_EQ(matrix.cols(), cols); // 1+3 cols
        lin_alg::Matrix<int> expected_matrix(rows, cols);
        expected_matrix << 1, 4, 7, 10,
            2, 5, 8, 11,
            3, 6, 9, 12;
        EXPECT_EQ(matrix, expected_matrix);
    }

    // throw out of range exceptions
    {
        lin_alg::ColVector<int> vector(rows);
        vector << 666, 666, 666;
        // at the right side of the matrix
        EXPECT_THROW(lin_alg::InsertCol(matrix, vector, 5), std::invalid_argument);
        // at the left side of the matrix
        EXPECT_THROW(lin_alg::InsertCol(matrix, vector, -1), std::invalid_argument);
    }
}

TEST(LinearAlgebra, InsertColFromSTLVector)
{
    Eigen::Index constexpr rows = 3;
    Eigen::Index cols = 1;
    lin_alg::Matrix<int> matrix(rows, 1);
    matrix << 7,
        8,
        9;

    // insert zeroth column
    {
        std::vector<int> vector{1, 2, 3};
        lin_alg::InsertCol(matrix, vector, 0);
        cols++;
        EXPECT_EQ(matrix.rows(), rows); // unchanged
        EXPECT_EQ(matrix.cols(), cols); // 1+1 cols
        lin_alg::Matrix<int> expected_matrix(rows, cols);
        expected_matrix << 1, 7,
            2, 8,
            3, 9;
        EXPECT_EQ(matrix, expected_matrix);
    }

    // insert first column
    {
        std::vector<int> vector{4, 5, 6};
        lin_alg::InsertCol(matrix, vector, 1);
        cols++;
        EXPECT_EQ(matrix.rows(), rows); // unchanged
        EXPECT_EQ(matrix.cols(), cols); // 1+2 cols
        lin_alg::Matrix<int> expected_matrix(rows, cols);
        expected_matrix << 1, 4, 7,
            2, 5, 8,
            3, 6, 9;
        EXPECT_EQ(matrix, expected_matrix);
    }

    // insert third column
    {
        std::vector<int> vector{10, 11, 12};
        lin_alg::InsertCol(matrix, vector, 3);
        cols++;
        EXPECT_EQ(matrix.rows(), rows); // unchanged
        EXPECT_EQ(matrix.cols(), cols); // 1+3 cols
        lin_alg::Matrix<int> expected_matrix(rows, cols);
        expected_matrix << 1, 4, 7, 10,
            2, 5, 8, 11,
            3, 6, 9, 12;
        EXPECT_EQ(matrix, expected_matrix);
    }

    // throw out of range exceptions
    {
        std::vector<int> vector{666, 666, 666};
        // at the right side of the matrix
        EXPECT_THROW(lin_alg::InsertCol(matrix, vector, 5), std::invalid_argument);
        // at the left side of the matrix
        EXPECT_THROW(lin_alg::InsertCol(matrix, vector, -1), std::invalid_argument);
    }
}

TEST(LinearAlgebra, STLVectorToEigenVector)
{
    std::vector<int> stl_vec{1, 2, 3, 4, 5};

    auto const row_vec = lin_alg::STLVectorToRowVector(stl_vec);
    for (Eigen::Index i = 0; i < row_vec.size(); ++i)
    {
        EXPECT_EQ(row_vec(i), stl_vec[i]);
    }

    auto const col_vec = lin_alg::STLVectorToColVector(stl_vec);
    for (Eigen::Index i = 0; i < col_vec.size(); ++i)
    {
        EXPECT_EQ(col_vec(i), stl_vec[i]);
    }

    EXPECT_EQ(row_vec, col_vec.transpose());
}

TEST(LinearAlgebra, SwapRows)
{
    lin_alg::Matrix<int> matrix(4, 3);
    matrix << 1, 2, 3,
        4, 5, 6,
        7, 8, 9,
        10, 11, 12;

    lin_alg::SwapRows(matrix, 1, 2);

    lin_alg::Matrix<int> expected_matrix(4, 3);
    expected_matrix << 1, 2, 3,
        7, 8, 9,
        4, 5, 6,
        10, 11, 12;

    EXPECT_EQ(matrix, expected_matrix);

    // nothing happens if indices are identical
    lin_alg::SwapRows(matrix, 0, 0);
    EXPECT_EQ(matrix, expected_matrix);
    lin_alg::SwapRows(matrix, 3, 3);
    EXPECT_EQ(matrix, expected_matrix);

    // go out of range
    EXPECT_THROW(lin_alg::SwapRows(matrix, 0, 4), std::invalid_argument);
    EXPECT_THROW(lin_alg::SwapRows(matrix, -1, 3), std::invalid_argument);
    EXPECT_THROW(lin_alg::SwapRows(matrix, -1, 4), std::invalid_argument);
}

TEST(LinearAlgebra, SwapColumns)
{
    lin_alg::Matrix<int> matrix(4, 3);
    matrix << 1, 2, 3,
        4, 5, 6,
        7, 8, 9,
        10, 11, 12;

    lin_alg::SwapColumns(matrix, 1, 2);

    lin_alg::Matrix<int> expected_matrix(4, 3);
    expected_matrix << 1, 3, 2,
        4, 6, 5,
        7, 9, 8,
        10, 12, 11;

    EXPECT_EQ(matrix, expected_matrix);

    // nothing happens if indices are identical
    lin_alg::SwapColumns(matrix, 0, 0);
    EXPECT_EQ(matrix, expected_matrix);
    lin_alg::SwapColumns(matrix, 2, 2);
    EXPECT_EQ(matrix, expected_matrix);

    // go out of range
    EXPECT_THROW(lin_alg::SwapColumns(matrix, 0, 3), std::invalid_argument);
    EXPECT_THROW(lin_alg::SwapColumns(matrix, -1, 2), std::invalid_argument);
    EXPECT_THROW(lin_alg::SwapColumns(matrix, -1, 3), std::invalid_argument);
}
