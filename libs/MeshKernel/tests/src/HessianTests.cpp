#include <gtest/gtest.h>

#include <algorithm>
#include <iomanip>
#include <random>
#include <vector>

#include "MeshKernel/AveragingInterpolation.hpp"
#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/HessianCalculator.hpp"
#include "MeshKernel/MeshRefinement.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/Point.hpp"
#include "MeshKernel/RidgeRefinement.hpp"

#include "TestUtils/Definitions.hpp"
#include "TestUtils/MakeMeshes.hpp"
#include "TestUtils/SampleFileReader.hpp"

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

#if 0

TEST(HessianTests, TidyPointsNoDuplicatesTest)
{
    // Test that the TidySamples does not change a data set that has no duplicates.

    const mk::UInt size = 100;

    std::vector<mk::Point> samplePoints(size);
    std::vector<mk::Point> generated;
    std::vector<double> sampleData(size);

    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_real_distribution<double> distribution(0.0, 100.0);

    for (mk::UInt i = 0; i < size; ++i)
    {
        samplePoints[i] = mk::Point(distribution(generator), distribution(generator));
        sampleData[i] = distribution(generator);
    }

    std::vector<mk::Point> expectedSamplePoints = samplePoints;
    std::vector<double> expectedSampleData = sampleData;

    mk::RidgeRefinement ridgeRefinement;
    ridgeRefinement.TidySamples(samplePoints, sampleData);

    for (mk::UInt i = 0; i < size; ++i)
    {
        EXPECT_EQ(expectedSamplePoints[i].x, samplePoints[i].x);
        EXPECT_EQ(expectedSamplePoints[i].y, samplePoints[i].y);
        EXPECT_EQ(expectedSampleData[i], sampleData[i]);
    }
}

TEST(HessianTests, TidyPointsTest)
{

    const mk::UInt size = 100;

    std::vector<mk::Point> samplePoints(size);
    std::vector<mk::Point> expectedSamplePoints;
    std::vector<mk::Point> generated;
    std::vector<double> sampleData(size);

    const double invMax = 1.0 / static_cast<double>(RAND_MAX);

    for (mk::UInt i = 0; i < size; ++i)
    {
        // Use the simple pseudo random number generator to have repeatable values between tests when debugging.
        samplePoints[i] = mk::Point(rand() * invMax, rand() * invMax);
        sampleData[i] = rand() * invMax;

        if (i > 0 && i % 10 == 0)
        {
            samplePoints[i - 1] = samplePoints[i];
            sampleData[i - 1] = sampleData[i];
        }
        else if (i > 10 && i % 21 == 0)
        {
            samplePoints[i - 10] = samplePoints[i];
            sampleData[i - 10] = sampleData[i];
        }
        else
        {
            expectedSamplePoints.push_back(samplePoints[i]);
        }
    }

    mk::RidgeRefinement ridgeRefinement;
    ridgeRefinement.TidySamples(samplePoints, sampleData);

    EXPECT_EQ(expectedSamplePoints.size(), samplePoints.size());

#if 0
    for (const auto& point : expectedSamplePoints)
    {
        EXPECT_TRUE(std::find(samplePoints.begin(), samplePoints.end(), point) != samplePoints.end());
    }
#endif
}

#endif

std::vector<meshkernel::Sample> generateSampleData(meshkernel::UInt nx, meshkernel::UInt ny, double deltaX, double deltaY)
{
    meshkernel::UInt start = 0;
    meshkernel::UInt size = (nx - start) * (nx - start);
    std::vector<meshkernel::Sample> sampleData(size);

    [[maybe_unused]] double centreX = (static_cast<double>((nx - 1) / 2) * deltaX);
    [[maybe_unused]] double centreY = (static_cast<double>((ny - 1) / 2) * deltaY);

    std::cout << "centre: " << centreX << "  " << centreY << std::endl;

    [[maybe_unused]] double xInit = 0.0; // 0.5 * deltaX;
    [[maybe_unused]] double yInit = 0.0; // 0.5 * deltaY;

    [[maybe_unused]] double scale = (ny / 4.0) * deltaY;
    [[maybe_unused]] double x = xInit; // 0.0;
    [[maybe_unused]] double y = yInit; // 0.0;
    [[maybe_unused]] double r = (nx / 5) * deltaX;
    [[maybe_unused]] double maxx = (nx - 1) * deltaX;

    meshkernel::UInt count = 0;

    // Does not seem to matter if this is defined or not
    //
#define Y_FIRST_SAMPLE

    for (meshkernel::UInt i = start; i < nx; ++i)
    {
#ifdef Y_FIRST_SAMPLE
        y = yInit;
#else
        x = xInit;
#endif

        for (meshkernel::UInt j = start; j < ny; ++j)
        {

#if 0
            // Gaussian bump, in the centre of the grid
            double centre = (x - centreX) * (x - centreX) + (y - centreY) * (y - centreY);
            double sample = 100.0 * std::exp(-0.025 * centre);
#elif 0
            // Gaussian wave, along the centre line of the grid
            double centre = (x - centreX) * (x - centreX);
            double sample = 100.0 * std::exp(-0.025 * centre);

#elif 1
            // sine wave ridge in x-direction
            double sinx = std::sin(x / maxx * M_PI * 4.0);
            double xxx = scale * sinx + centreY;
            double sample = 10 * (std::atan(20.0 * (xxx - y)) + M_PI / 2.0);
#else

            // A arctan function
            double centre = (x - centreX) * (x - centreX) + (y - centreY) * (y - centreY);
            double sample = 10 * (std::atan(20.0 * (r * r - centre)) + M_PI / 2.0);
#endif

            sampleData[count] = {x, y, sample};

#ifdef Y_FIRST_SAMPLE
            y += deltaY;
#else
            x += deltaX;
#endif

            ++count;
        }

#ifdef Y_FIRST_SAMPLE
        x += deltaX;
#else
        y += deltaY;
#endif
    }

    return sampleData;
}

TEST(HessianTests, CheckHessian)
{
    meshkernel::UInt nx = 41; // 101;
    meshkernel::UInt ny = 21; // 101;

    double deltaX = 10.0;
    double deltaY = 10.0;

    double dimX = (nx - 1) * deltaX;
    double dimY = (ny - 1) * deltaY;

    std::shared_ptr<meshkernel::Mesh2D> mesh = MakeRectangularMeshForTesting(nx, ny, dimX, dimY, meshkernel::Projection::cartesian);

    meshkernel::UInt superSample = 2;

    meshkernel::UInt sampleNx = (nx - 1) * superSample + 1;
    meshkernel::UInt sampleNy = (ny - 1) * superSample + 1;

    double sampleDeltaX = deltaX / static_cast<double>(superSample);
    double sampleDeltaY = deltaY / static_cast<double>(superSample);

    [[maybe_unused]] double radius = 2.0 * sampleDeltaX;

    std::vector<meshkernel::Sample> sampleData = generateSampleData(sampleNx, sampleNy, sampleDeltaX, sampleDeltaY);

#if 0
    std::shared_ptr<meshkernel::HessianAveragingInterpolation> hessian = std::make_shared<meshkernel::HessianAveragingInterpolation>(*mesh,
                                                                                                                                     sampleData,
                                                                                                                                     sampleNx,
                                                                                                                                     sampleNy,
                                                                                                                                     meshkernel::AveragingInterpolation::Method::SimpleAveraging,
                                                                                                                                     meshkernel::Mesh::Location::Faces,
                                                                                                                                     radius,
                                                                                                                                     false,
                                                                                                                                     false,
                                                                                                                                     1);

    std::vector<meshkernel::Sample> samples = hessian->hessianSamples();
#elif 1
    std::vector<meshkernel::Sample> samples;
    meshkernel::HessianCalculator::Compute(sampleData, mesh->m_projection, sampleNx, sampleNy, samples);
#else
    std::vector<meshkernel::Sample> samples(sampleData.size());
    meshkernel::HessianCalculator::Compute(sampleData, mesh->m_projection, sampleNx, sampleNy, samples);
#endif
    const auto interpolator = std::make_shared<meshkernel::AveragingInterpolation>(*mesh,
                                                                                   samples,
                                                                                   // meshkernel::AveragingInterpolation::Method::Max,
                                                                                   meshkernel::AveragingInterpolation::Method::MinAbsValue,
                                                                                   // meshkernel::AveragingInterpolation::Method::InverseWeightedDistance,
                                                                                   // meshkernel::AveragingInterpolation::Method::Closest,
                                                                                   // meshkernel::AveragingInterpolation::Method::SimpleAveraging,
                                                                                   meshkernel::Mesh::Location::Faces,
                                                                                   1.0,
                                                                                   false,
                                                                                   false,
                                                                                   1);

    meshkernel::MeshRefinementParameters meshRefinementParameters;

    meshRefinementParameters.max_num_refinement_iterations = 2;
    meshRefinementParameters.max_num_refinement_iterations = 1;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 2.0;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 3;
    meshRefinementParameters.smoothing_iterations = 0;

    meshkernel::MeshRefinement meshRefinement(mesh, interpolator, meshRefinementParameters, false);

    // meshkernel::Print(mesh->m_nodes, mesh->m_edges);

    // Execute
    meshRefinement.Compute();

    std::cout << std::endl;
    std::cout << "--------------------------------" << std::endl;
    std::cout << std::endl;

    meshkernel::Print(mesh->m_nodes, mesh->m_edges);
}
