#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Exceptions.hpp"

std::tuple<int, double*, double*> meshkernel::ConvertFromNodesVector(const std::vector<Point>& nodes)
{

    double* xNodes = nullptr;
    double* yNodes = nullptr;

    try
    {
        if (nodes.empty())
        {
            return {0, nullptr, nullptr};
        }

        xNodes = new double[nodes.size()];
        yNodes = new double[nodes.size()];

        for (size_t i = 0; i < nodes.size(); ++i)
        {
            xNodes[i] = nodes[i].x;
            yNodes[i] = nodes[i].y;
        }

        return {static_cast<int>(nodes.size()), xNodes, yNodes};
    }
    catch (...)
    {
        if (xNodes != nullptr)
        {
            delete[] xNodes;
        }
        if (yNodes != nullptr)
        {
            delete[] yNodes;
        }

        throw;
    }
}
