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
