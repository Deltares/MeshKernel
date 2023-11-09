#include "SampleFileReader.hpp"

#include <gtest/gtest.h>

#include <vector>

#include "MeshKernel/AveragingInterpolation.hpp"
#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/HessianCalculator.hpp"
#include "MeshKernel/MeshRefinement.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/Point.hpp"
#include <TestUtils/Definitions.hpp>

#include "TestUtils/MakeMeshes.hpp"
#include "TestUtils/SampleFileWriter.hpp"

#include <fstream>

namespace mk = meshkernel;

std::vector<std::vector<mk::Point>> GenerateGridPoints(const mk::UInt rows, const mk::UInt cols)
{

    std::vector meshPoints(cols, std::vector<mk::Point>(rows));

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

std::vector<meshkernel::Sample> generateSampleData(meshkernel::UInt nx, meshkernel::UInt ny, double deltaX, double deltaY)
{
    meshkernel::UInt start = 0;
    meshkernel::UInt size = (nx - start) * (ny - start);
    std::vector<meshkernel::Sample> sampleData(size);

    std::vector sampleDataMatrix(ny, std::vector<double>(nx));

    [[maybe_unused]] double centreX = (static_cast<double>((nx - 1) / 2) * deltaX);
    [[maybe_unused]] double centreY = (static_cast<double>((ny - 1) / 2) * deltaY);

    [[maybe_unused]] double scale = (ny / 4.0) * deltaY;

    [[maybe_unused]] double r = (nx / 5) * deltaX;
    [[maybe_unused]] double maxx = (nx - 1) * deltaX;

    for (int i = ny - 1; i >= 0; --i)
    {

        for (meshkernel::UInt j = start; j < nx; ++j)
        {
            const double y = deltaY * i;
            const double x = deltaX * j;

#if 0
            // Gaussian bump, in the centre of the grid
            double centre = (x - centreX) * (x - centreX) + (y - centreY) * (y - centreY);
            double sample = 100.0 * std::exp(-0.025 * centre);
#elif 0
            // Gaussian wave, along the centre line of the grid
            double centre = (x - centreX) * (x - centreX) + (y - centreY) * (y - centreY);
            double factor = std::max(1e-6, std::exp(-0.00025 * centre));
            double sample = 100.0 * factor;

#elif 0
            // sine wave ridge in x-direction
            double sinx = std::sin(x / maxx * M_PI * 4.0);
            double xxx = scale * sinx + centreY;
            double sample = 10 * (std::atan(20.0 * (xxx - y)) + M_PI / 2.0);

#else

            // A arctan function
            double centre = (x - centreX) * (x - centreX) + (y - centreY) * (y - centreY);
            double sample = 10 * (std::atan(20.0 * (r * r - centre)) + M_PI / 2.0);
#endif

            sampleDataMatrix[(ny - 1) - i][j] = sample;
        }
    }

    meshkernel::UInt count = 0;
    for (meshkernel::UInt j = start; j < nx; ++j)
    {
        for (int i = ny - 1; i >= 0; --i)
        {
            const double y = deltaY * i;
            const double x = deltaX * j;
            sampleData[count] = {x, y, sampleDataMatrix[(ny - 1) - i][j]};
            count++;
        }
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

    const double sampleDeltaX = deltaX / static_cast<double>(superSample);
    const double sampleDeltaY = deltaY / static_cast<double>(superSample);

    const auto sampleData = generateSampleData(sampleNx, sampleNy, sampleDeltaX, sampleDeltaY);

    meshkernel::UInt numberOfSmoothingIterations = 0;

    auto samples = meshkernel::HessianCalculator::ComputeHessianSamples(sampleData, mesh->m_projection, numberOfSmoothingIterations, sampleNx, sampleNy);

    const auto interpolator = std::make_shared<meshkernel::AveragingInterpolation>(*mesh,
                                                                                   samples,
                                                                                   meshkernel::AveragingInterpolation::Method::Max,
                                                                                   meshkernel::Mesh::Location::Faces,
                                                                                   1.0,
                                                                                   false,
                                                                                   false,
                                                                                   1);

    meshkernel::MeshRefinementParameters meshRefinementParameters;

    meshRefinementParameters.max_num_refinement_iterations = 3;
    meshRefinementParameters.refine_intersected = 0;
    meshRefinementParameters.use_mass_center_when_refining = 0;
    meshRefinementParameters.min_edge_size = 2.0;
    meshRefinementParameters.account_for_samples_outside = 0;
    meshRefinementParameters.connect_hanging_nodes = 1;
    meshRefinementParameters.refinement_type = 3;
    meshRefinementParameters.smoothing_iterations = 0;

    meshkernel::MeshRefinement meshRefinement(mesh, interpolator, meshRefinementParameters, false);

    // Execute
    meshRefinement.Compute();

    // std::ofstream outputFile;

    // const auto resultFolder = TEST_FOLDER + "/result.m";
    // outputFile.open(resultFolder);

    // meshkernel::Print(mesh->m_nodes, mesh->m_edges, outputFile);
    // outputFile.close();
}
