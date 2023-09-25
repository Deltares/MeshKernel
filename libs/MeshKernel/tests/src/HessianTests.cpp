#include <gtest/gtest.h>

#include <vector>

#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Point.hpp"
#include "MeshKernel/RidgeRefinement.hpp"

namespace mk = meshkernel;

std::vector<std::vector<mk::Point>> GenerateGridPoints(const mk::UInt rows, const mk::UInt cols)
{

    std::vector<std::vector<mk::Point>> meshPoints(cols, std::vector<mk::Point>(rows));

    mk::Point origin = mk::Point(0.0, 0.0);

    double hx = 1.0;
    double hy = 1.0;

    mk::Point current = origin;

    for (mk::UInt ix = 0; ix < cols; ++ix)
    {
        current.y = origin.y;

        for (mk::UInt jy = 0; jy < rows; ++jy)
        {
            meshPoints[ix][jy] = current;
            current.y += hy;
        }

        current.x += hx;
    }

    return meshPoints;
}

TEST(HessianTests, StructuredGridTest)
{

    const mk::UInt rows = 10;
    const mk::UInt cols = 10;

    std::vector<std::vector<mk::Point>> gridPoints(GenerateGridPoints(rows, cols));

    std::vector<mk::Point> samplePoints(rows * cols);
    std::vector<double> sampleData(rows * cols);

    mk::UInt count = 0;
    const double invMax = 1.0 / static_cast<double>(RAND_MAX);

    for (mk::UInt j = 0; j < cols; ++j)
    {
        for (mk::UInt i = 0; i < rows; ++i)
        {
            samplePoints[count] = gridPoints[i][j];
            sampleData[count] = rand() * invMax;
            ++count;
        }
    }

    mk::RidgeRefinement rigdeRefinement;
    rigdeRefinement.Compute(samplePoints, sampleData, mk::Projection::cartesian, rows, cols);
}
