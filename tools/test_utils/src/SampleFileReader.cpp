#include "SampleFileReader.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

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

    int numX = 0;
    int numY = 0;
    double xllCenter = 0.0;
    double yllCenter = 0.0;
    double cellSize = 0.0;
    double nodataValue;

    int numlines = 0;
    std::vector<std::vector<double>> rows;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        if (numlines == 0)
        {
            std::string s;
            iss >> s >> numX;
            numlines++;
            continue;
        }
        if (numlines == 1)
        {
            std::string s;
            iss >> s >> numY;
            numlines++;
            continue;
        }
        if (numlines == 2)
        {
            std::string s;
            iss >> s >> xllCenter;
            numlines++;
            continue;
        }
        if (numlines == 3)
        {
            std::string s;
            iss >> s >> yllCenter;
            numlines++;
            continue;
        }
        if (numlines == 4)
        {
            std::string s;
            iss >> s >> cellSize;
            numlines++;
            continue;
        }
        if (numlines == 5)
        {
            std::string s;
            iss >> s >> nodataValue;
            numlines++;
            continue;
        }

        rows.push_back(std::vector<double>());
        for (auto i = 0; i < numX; ++i)
        {
            double value;
            iss >> value;
            rows.back().push_back(value);
        }
        numlines++;
    }

    std::reverse(rows.begin(), rows.end());
    std::vector<double> values;
    for (size_t i = 0; i < rows.size(); ++i)
    {
        for (auto j = 0u; j < rows[i].size(); ++j)
        {
            values.push_back(rows[i][j]);
        }
    }

    return {numX, numY, xllCenter, yllCenter, cellSize, nodataValue, values};
}
