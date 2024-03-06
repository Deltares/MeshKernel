#include "MeshKernel/OrthogonalizationAndSmoothingAction.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Mesh2D.hpp"

#include <algorithm>

std::unique_ptr<meshkernel::OrthogonalizationAndSmoothingAction> meshkernel::OrthogonalizationAndSmoothingAction::Create(Mesh2D& mesh, const std::vector<UInt>& nodeIndices)
{
    return std::make_unique<OrthogonalizationAndSmoothingAction>(mesh, nodeIndices);
}

meshkernel::OrthogonalizationAndSmoothingAction::OrthogonalizationAndSmoothingAction(Mesh2D& mesh, const std::vector<UInt>& nodeIndices)
    : BaseMeshUndoAction<OrthogonalizationAndSmoothingAction, Mesh2D>(mesh)
{

    if (nodeIndices.size() * (sizeof(UInt) + sizeof(Point)) > mesh.GetNumNodes() * sizeof(Point))
    {
        // save all nodes and only the nodes. No need to save the node indices
        // Based on approximation of total memory required for the action.
        m_nodes = mesh.Nodes();
    }
    else
    {
        // save only nodes to be moved.
        m_nodes.resize(nodeIndices.size());
        m_nodeIndices = nodeIndices;

        // Copy nodes to be saved from mesh.
        std::transform(nodeIndices.begin(), nodeIndices.end(), m_nodes.begin(), [&mesh = mesh](UInt pos)
                       { return mesh.Node(pos); });
    }
}

void meshkernel::OrthogonalizationAndSmoothingAction::Swap(std::vector<Point>& nodes)
{
    if (nodes.size() < m_nodes.size())
    {
        throw ConstraintError("Number of nodes passed is less than nodes stored. {} < {}",
                              nodes.size(), m_nodes.size());
    }

    if (m_nodeIndices.empty())
    {
        std::swap_ranges(m_nodes.begin(), m_nodes.end(), nodes.begin());
    }
    else
    {
        for (UInt i = 0; i < m_nodeIndices.size(); ++i)
        {
            std::swap(m_nodes[i], nodes[m_nodeIndices[i]]);
        }
    }
}

std::uint64_t meshkernel::OrthogonalizationAndSmoothingAction::MemorySize() const
{
    return sizeof(*this) + m_nodes.size() * sizeof(Point) + m_nodeIndices.size() * sizeof(UInt);
}

void meshkernel::OrthogonalizationAndSmoothingAction::Print(std::ostream& out) const
{
    out << "OrthogonalizationAndSmoothingAction: state " << to_string(State())
        << ", number of nodes: " << m_nodes.size()
        << std::endl;
}
