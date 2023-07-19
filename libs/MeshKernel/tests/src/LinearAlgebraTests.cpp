#include <gtest/gtest.h>

#include <MeshKernel/Utilities/LinearAlgebra.hpp>

using namespace meshkernel;

TEST(LinearAlgebra, Empty)
{
    {
        lin_alg::MatrixRowMajor<int> matrix;
        EXPECT_TRUE(lin_alg::MatrixIsEmpty(matrix));
    }

    {
        lin_alg::MatrixRowMajor<int> matrix(3, 4);
        EXPECT_FALSE(lin_alg::MatrixIsEmpty(matrix));
    }
}

TEST(LinearAlgebra, Resize)
{
    lin_alg::MatrixRowMajor<int> matrix;

    {
        Eigen::Index constexpr rows = 3;
        Eigen::Index constexpr cols = 4;

        EXPECT_TRUE(lin_alg::ResizeAndFillMatrix(matrix, rows, cols)); // 3x4
        EXPECT_EQ(matrix.rows(), rows);
        EXPECT_EQ(matrix.cols(), cols);
        // fill value was not set, expect default int initialisation

        lin_alg::MatrixRowMajor<int> expected_matrix;
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
        lin_alg::MatrixRowMajor<int> expected_matrix(rows_new, cols_new); // 4x5
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

        lin_alg::MatrixRowMajor<int> expected_matrix(rows_new, cols_new); // 5x6
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
            lin_alg::MatrixRowMajor<int> expected_matrix(rows_new, cols_new); // 4x5
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
            lin_alg::MatrixRowMajor<int> expected_matrix(rows_new, cols_new); // 4x8
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
            lin_alg::MatrixRowMajor<int> expected_matrix(rows_new, cols_new); // 8x5
            expected_matrix.setOnes();
            expected_matrix.rightCols(3).setConstant(666);
            expected_matrix.bottomRows(4).setConstant(-666);
            EXPECT_EQ(matrix, expected_matrix);
        }
    }
}

TEST(LinearAlgebra, EraseOneRow)
{
    lin_alg::MatrixRowMajor<int> matrix(5, 3);
    matrix << 1, 2, 3,
        4, 5, 6,
        7, 8, 9,
        10, 11, 12,
        13, 14, 15;

    lin_alg::EraseRow(matrix, 2);
    EXPECT_EQ(matrix.rows(), 4);
    {
        lin_alg::MatrixRowMajor<int> expected_matrix(4, 3);
        expected_matrix << 1, 2, 3,
            4, 5, 6,
            10, 11, 12,
            13, 14, 15;
        EXPECT_EQ(matrix, expected_matrix);
    }

    lin_alg::EraseRow(matrix, 1);
    EXPECT_EQ(matrix.rows(), 3);
    {
        lin_alg::MatrixRowMajor<int> expected_matrix(3, 3);
        expected_matrix << 1, 2, 3,
            10, 11, 12,
            13, 14, 15;
        EXPECT_EQ(matrix, expected_matrix);
    }

    lin_alg::EraseRow(matrix, 0);
    EXPECT_EQ(matrix.rows(), 2);
    {
        lin_alg::MatrixRowMajor<int> expected_matrix(2, 3);
        expected_matrix << 10, 11, 12,
            13, 14, 15;
        EXPECT_EQ(matrix, expected_matrix);
    }
}

TEST(LinearAlgebra, EraseMultipleRows)
{
    Eigen::Index constexpr cols = 3;
    lin_alg::MatrixRowMajor<int> matrix(6, cols);
    matrix << 1, 2, 3,
        4, 5, 6,
        7, 8, 9,
        10, 11, 12,
        13, 14, 15,
        16, 17, 18;

    lin_alg::EraseRows(matrix, 2, 3); // remove rows 2 and 3
    EXPECT_EQ(matrix.rows(), 4);      // 6 -2 = 4 rows left
    EXPECT_EQ(matrix.cols(), cols);   // unchanged
    lin_alg::MatrixRowMajor<int> expected_matrix(3, 3);
    expected_matrix << 1, 2, 3,
        4, 5, 6,
        13, 14, 15,
        16, 17, 18;
    EXPECT_EQ(matrix, expected_matrix);
}

TEST(LinearAlgebra, EraseOneCol)
{
    lin_alg::MatrixRowMajor<int> matrix(3, 5);
    matrix << 1, 2, 3, 4, 5,
        6, 7, 8, 9, 10,
        11, 12, 13, 14, 15;

    lin_alg::EraseCol(matrix, 2);
    EXPECT_EQ(matrix.cols(), 4);
    {
        lin_alg::MatrixRowMajor<int> expected_matrix(3, 4);
        expected_matrix << 1, 2, 4, 5,
            6, 7, 9, 10,
            11, 12, 14, 15;
        EXPECT_EQ(matrix, expected_matrix);
    }

    lin_alg::EraseCol(matrix, 1);
    EXPECT_EQ(matrix.cols(), 3);
    {
        lin_alg::MatrixRowMajor<int> expected_matrix(3, 3);
        expected_matrix << 1, 4, 5,
            6, 9, 10,
            11, 14, 15;
        EXPECT_EQ(matrix, expected_matrix);
    }

    lin_alg::EraseCol(matrix, 0);
    EXPECT_EQ(matrix.cols(), 2);
    {
        lin_alg::MatrixRowMajor<int> expected_matrix(3, 2);
        expected_matrix << 4, 5,
            9, 10,
            14, 15;
        EXPECT_EQ(matrix, expected_matrix);
    }
}

TEST(LinearAlgebra, EraseMultipleCols)
{
    Eigen::Index constexpr rows = 3;
    lin_alg::MatrixRowMajor<int> matrix(rows, 6);
    matrix << 1, 2, 3, 4, 5, 6,
        7, 8, 9, 10, 11, 12,
        13, 14, 15, 16, 17, 18;

    lin_alg::EraseCols(matrix, 2, 3); // remove cols 2 and 3
    EXPECT_EQ(matrix.rows(), rows);   // unchanged
    EXPECT_EQ(matrix.cols(), 4);      // 6 -2 = 4 cols left
    lin_alg::MatrixRowMajor<int> expected_matrix(3, 3);
    expected_matrix << 1, 2, 5, 6,
        7, 8, 11, 12,
        13, 14, 17, 18;
    EXPECT_EQ(matrix, expected_matrix);
}

TEST(LinearAlgebra, NewTemplate)
{
    lin_alg::Matrix<int, Eigen::RowMajor> matrix(3, 4);
    /*matrix << 1, 2, 3, 4,
        5, 6, 7, 8,
        9, 10, 11, 12;*/

    matrix.setRandom();
}