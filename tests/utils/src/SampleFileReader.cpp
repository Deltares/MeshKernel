#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include <MeshKernel/Entities.hpp>
#include <TestUtils/SampleFileReader.hpp>

std::vector<meshkernel::Sample> ReadSampleFile(std::string filePath)
{
    std::vector<meshkernel::Sample> samples;

    // read sample file
    std::string line;
    std::ifstream infile(filePath);
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
