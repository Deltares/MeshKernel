#include "SampleGenerator.hpp"

std::vector<meshkernel::Sample>
generateSampleData(FunctionTestCase testcase,
                   meshkernel::UInt nx,
                   meshkernel::UInt ny,
                   double deltaX,
                   double deltaY)
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
    case FunctionTestCase::GaussianBump:
        generateSample = [&](double x, double y)
        {
            const double centre = (x - centreX) * (x - centreX) + (y - centreY) * (y - centreY);
            return 100.0 * std::exp(-0.025 * centre);
        };
        break;
    case FunctionTestCase::GaussianWave:
        generateSample = [&](double x, double y)
        {
            const double centre = (x - centreX) * (x - centreX) + (y - centreY) * (y - centreY);
            const double factor = std::max(1e-6, std::exp(-0.00025 * centre));
            return 100.0 * factor;
        };
        break;
    case FunctionTestCase::RidgeXDirection:
        generateSample = [&](double x, double y)
        {
            const double sinx = std::sin(x / maxx * M_PI * 4.0);
            const double xxx = scale * sinx + centreY;
            return 10 * (std::atan(20.0 * (xxx - y)) + M_PI / 2.0);
        };
        break;
    case FunctionTestCase::ArctanFunction:
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
