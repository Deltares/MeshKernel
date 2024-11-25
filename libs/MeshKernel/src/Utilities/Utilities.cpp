#include "MeshKernel/Utilities/Utilities.hpp"

void meshkernel::Print(const std::vector<Point>& nodes, const std::vector<Edge>& edges, std::ostream& out)
{
    out << "nullId = " << constants::missing::uintValue << std::endl;
    out << "nullValue = " << constants::missing::doubleValue << std::endl;

    out << "nodex = zeros ( " << nodes.size() << ", 1);" << std::endl;
    out << "nodey = zeros ( " << nodes.size() << ", 1);" << std::endl;
    out << "edges = zeros ( " << edges.size() << ", 2);" << std::endl;

    for (UInt i = 0; i < nodes.size(); ++i)
    {
        out << "nodex (" << i + 1 << " ) = " << nodes[i].x << ";" << std::endl;
    }

    for (UInt i = 0; i < nodes.size(); ++i)
    {
        out << "nodey (" << i + 1 << " ) = " << nodes[i].y << ";" << std::endl;
    }

    out << "edges = zeros ( " << edges.size() << ", 2 );" << std::endl;

    for (UInt i = 0; i < edges.size(); ++i)
    {
        out << "edges ( " << i + 1 << ", 1 ) = " << edges[i].first + 1 << ";" << std::endl;
        out << "edges ( " << i + 1 << ", 2 ) = " << edges[i].second + 1 << ";" << std::endl;
    }
}
