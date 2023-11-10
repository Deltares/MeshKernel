#include "Definitions.hpp"
#include "SampleFileReader.hpp"

#include <gtest/gtest.h>
#include <vector>

#include "MeshKernel/AveragingInterpolation.hpp"
#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/HessianCalculator.hpp"
#include "MeshKernel/MeshRefinement.hpp"
#include "MeshKernel/Operations.hpp"
#include "TestUtils/MakeMeshes.hpp"

namespace mk = meshkernel;

enum class RidgeRefinementTestCase
{
    GaussianBump = 1,
    GaussianWave = 2,
    RidgeXdirection = 3,
    ArctanFunction = 4,
};

std::vector<meshkernel::Sample> generateSampleData(RidgeRefinementTestCase testcase,
                                                   meshkernel::UInt nx = 10,
                                                   meshkernel::UInt ny = 10,
                                                   double deltaX = 10.0,
                                                   double deltaY = 10.0)
{
    meshkernel::UInt start = 0;
    meshkernel::UInt size = (nx - start) * (ny - start);
    std::vector<meshkernel::Sample> sampleData(size);

    std::vector sampleDataMatrix(ny, std::vector<double>(nx));

    const double centreX = static_cast<double>((nx - 1) / 2) * deltaX;
    const double centreY = static_cast<double>((ny - 1) / 2) * deltaY;

    const double scale = ny / 4.0 * deltaY;

    const double r = nx / 5 * deltaX;
    const double maxx = (nx - 1) * deltaX;

    std::function<double(double, double)> generateSample;
    switch (testcase)
    {
    case RidgeRefinementTestCase::GaussianBump:
        generateSample = [&](double x, double y)
        {
            const double centre = (x - centreX) * (x - centreX) + (y - centreY) * (y - centreY);
            return 100.0 * std::exp(-0.025 * centre);
        };
        break;
    case RidgeRefinementTestCase::GaussianWave:
        generateSample = [&](double x, double y)
        {
            const double centre = (x - centreX) * (x - centreX) + (y - centreY) * (y - centreY);
            const double factor = std::max(1e-6, std::exp(-0.00025 * centre));
            return 100.0 * factor;
        };
        break;
    case RidgeRefinementTestCase::RidgeXdirection:
        generateSample = [&](double x, double y)
        {
            const double sinx = std::sin(x / maxx * M_PI * 4.0);
            const double xxx = scale * sinx + centreY;
            return 10 * (std::atan(20.0 * (xxx - y)) + M_PI / 2.0);
        };
        break;
    case RidgeRefinementTestCase::ArctanFunction:
        generateSample = [&](double x, double y)
        {
            const double centre = (x - centreX) * (x - centreX) + (y - centreY) * (y - centreY);
            return 10 * (std::atan(20.0 * (r * r - centre)) + M_PI / 2.0);
        };
        break;
    default:
        throw std::invalid_argument("invalid ridge refinement test case");
    }

    for (int i = ny - 1; i >= 0; --i)
    {

        for (meshkernel::UInt j = start; j < nx; ++j)
        {
            const double y = deltaY * i;
            const double x = deltaX * j;
            sampleDataMatrix[ny - 1 - i][j] = generateSample(x, y);
        }
    }

    meshkernel::UInt count = 0;
    for (meshkernel::UInt j = start; j < nx; ++j)
    {
        for (int i = ny - 1; i >= 0; --i)
        {
            const double y = deltaY * i;
            const double x = deltaX * j;
            sampleData[count] = {x, y, sampleDataMatrix[ny - 1 - i][j]};
            count++;
        }
    }

    return sampleData;
}

class RidgeRefinementTestCasesParameters : public ::testing::TestWithParam<std::tuple<RidgeRefinementTestCase, meshkernel::UInt, meshkernel::UInt>>
{
public:
    [[nodiscard]] static std::vector<std::tuple<RidgeRefinementTestCase, meshkernel::UInt, meshkernel::UInt>> GetData()
    {
        return std::vector{
            std::make_tuple<RidgeRefinementTestCase, meshkernel::UInt, meshkernel::UInt>(RidgeRefinementTestCase::GaussianBump, 1165, 2344),
            std::make_tuple<RidgeRefinementTestCase, meshkernel::UInt, meshkernel::UInt>(RidgeRefinementTestCase::GaussianWave, 5297, 10784),
            std::make_tuple<RidgeRefinementTestCase, meshkernel::UInt, meshkernel::UInt>(RidgeRefinementTestCase::RidgeXdirection, 2618, 5694),
            std::make_tuple<RidgeRefinementTestCase, meshkernel::UInt, meshkernel::UInt>(RidgeRefinementTestCase::ArctanFunction, 2309, 5028)};
    }
};

TEST_P(RidgeRefinementTestCasesParameters, expectedResults)
{
    // Prepare
    const auto [testCase, numNodes, numEdges] = GetParam();

    meshkernel::UInt nx = 41;
    meshkernel::UInt ny = 21;

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

    const auto sampleData = generateSampleData(testCase, sampleNx, sampleNy, sampleDeltaX, sampleDeltaY);

    auto samples = meshkernel::HessianCalculator::ComputeHessianSamples(sampleData, mesh->m_projection, 0, sampleNx, sampleNy);

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

    // Assert
    ASSERT_EQ(numNodes, mesh->GetNumNodes());
    ASSERT_EQ(numEdges, mesh->GetNumEdges());
}

INSTANTIATE_TEST_SUITE_P(RidgeRefinementTestCases,
                         RidgeRefinementTestCasesParameters,
                         ::testing::ValuesIn(RidgeRefinementTestCasesParameters::GetData()));
