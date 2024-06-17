#include "MeshReaders.hpp"

#include <filesystem>
#include <fstream>
#include <sstream>

std::vector<meshkernel::Point> ReadPointFile(std::string const& filePath)
{
    std::vector<meshkernel::Point> points;

    // read point file
    std::string line;
    std::ifstream infile(filePath.c_str());

    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        double pointX;
        double pointY;

        iss >> pointX;
        iss >> pointY;

        points.push_back({pointX, pointY});
    }

    return points;
}

std::vector<meshkernel::Edge> ReadEdgeFile(std::string const& filePath)
{
    std::vector<meshkernel::Edge> edges;

    // read edge file
    std::string line;
    std::ifstream infile(filePath.c_str());

    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        meshkernel::UInt edgeStart;
        meshkernel::UInt edgeEnd;

        iss >> edgeStart;
        iss >> edgeEnd;

        edges.push_back({edgeStart, edgeEnd});
    }

    return edges;
}
