#include "SplineReader.hpp"

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Point.hpp"

meshkernel::Splines LoadSplines(const std::string& fileName)
{
    std::ifstream splineFile;
    splineFile.open(fileName.c_str());

    meshkernel::Splines splines(meshkernel::Projection::cartesian);
    std::string line;

    std::vector<meshkernel::Point> splinePoints;

    while (std::getline(splineFile, line))
    {
        if (size_t found = line.find("L00"); found != std::string::npos)
        {
            std::getline(splineFile, line);
            std::istringstream sizes(line);

            meshkernel::UInt numPoints = 0;
            meshkernel::UInt numDim = 0;

            sizes >> numPoints;
            sizes >> numDim;

            splinePoints.clear();
            splinePoints.reserve(numPoints);

            for (meshkernel::UInt i = 0; i < numPoints; ++i)
            {
                std::getline(splineFile, line);
                std::istringstream values(line);
                double x;
                double y;
                values >> x;
                values >> y;
                splinePoints.emplace_back(meshkernel::Point(x, y));
            }

            splines.AddSpline(splinePoints);
        }
    }

    splineFile.close();
    return splines;
}

std::vector<meshkernel::Point> LoadSplinePoints(const std::string& fileName)
{
    std::ifstream splineFile;
    splineFile.open(fileName.c_str());

    std::vector<meshkernel::Point> splinePoints;
    std::string line;

    bool firstSpline = true;

    while (std::getline(splineFile, line))
    {
        if (size_t found = line.find("L00"); found != std::string::npos)
        {
            std::getline(splineFile, line);
            std::istringstream sizes(line);

            meshkernel::UInt numPoints = 0;
            meshkernel::UInt numDim = 0;

            sizes >> numPoints;
            sizes >> numDim;

            if (!firstSpline)
            {
                splinePoints.emplace_back(meshkernel::Point(meshkernel::constants::missing::doubleValue, meshkernel::constants::missing::doubleValue));
            }

            for (meshkernel::UInt i = 0; i < numPoints; ++i)
            {
                std::getline(splineFile, line);
                std::istringstream values(line);
                double x;
                double y;
                values >> x;
                values >> y;
                splinePoints.emplace_back(meshkernel::Point(x, y));
            }

            firstSpline = false;
        }
    }

    splineFile.close();
    return splinePoints;
}
