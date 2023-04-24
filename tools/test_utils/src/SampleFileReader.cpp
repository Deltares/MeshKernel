#include "MeshKernelApi/GriddedSamples.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <MeshKernel/Entities.hpp>
#include <TestUtils/SampleFileReader.hpp>

std::vector<meshkernel::Sample> ReadSampleFile(std::string const& filePath)
{
    std::vector<meshkernel::Sample> samples;

    // read sample file
    std::string line;
    std::ifstream infile(filePath.c_str());
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        double sampleX;
        double sampleY;
        double sampleValue;

        iss >> sampleX;
        iss >> sampleY;
        iss >> sampleValue;

        samples.push_back({sampleX, sampleY, sampleValue});
    }

    return samples;
}

std::tuple<int, int, double, double, double, double, std::vector<double>> ReadAscFile(const std::string& filePath)
{
    // read sample file
    std::string line;
    std::ifstream infile(filePath.c_str());

    int ncols = 0;
    int nrows = 0;
    double xllcenter = 0.0;
    double yllcenter = 0.0;
    double cellsize = 0.0;
    double nodata_value;

    int numlines = 0;
    std::vector<std::vector<double>> rows;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        if (numlines == 0)
        {
            std::string s;
            iss >> s >> ncols;
            numlines++;
            continue;
        }
        if (numlines == 1)
        {
            std::string s;
            iss >> s >> nrows;
            numlines++;
            continue;
        }
        if (numlines == 2)
        {
            std::string s;
            iss >> s >> xllcenter;
            numlines++;
            continue;
        }
        if (numlines == 3)
        {
            std::string s;
            iss >> s >> yllcenter;
            numlines++;
            continue;
        }
        if (numlines == 4)
        {
            std::string s;
            iss >> s >> cellsize;
            numlines++;
            continue;
        }
        if (numlines == 5)
        {
            std::string s;
            iss >> s >> nodata_value;
            numlines++;
            continue;
        }

        rows.push_back(std::vector<double>());
        for (auto i = 0; i < ncols + 1; ++i)
        {
            double value;
            iss >> value;
            rows.back().push_back(value);
        }
        numlines++;
    }

    std::reverse(rows.begin(), rows.end());
    std::vector<double> values;
    for (auto i = 0u; i < rows.size(); ++i)
    {
        for (auto j = 0u; j < rows[i].size(); ++j)
        {
            values.push_back(rows[i][j]);
        }
    }

    return {ncols, nrows, xllcenter, yllcenter, cellsize, nodata_value, values};
}
