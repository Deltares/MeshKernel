#include "MeshKernel/AddNodeAction.hpp"
#include "MeshKernel/Mesh.hpp"

std::unique_ptr<meshkernel::AddNodeAction> meshkernel::AddNodeAction::create(Mesh& mesh, const UInt id, const Point& point)
{
    return std::make_unique<AddNodeAction>(mesh, id, point);
}

meshkernel::AddNodeAction::AddNodeAction(Mesh& mesh, const UInt id, const Point& p) : BaseMeshUndoAction<AddNodeAction, Mesh>(mesh), m_nodeId(id), m_node(p) {}

meshkernel::UInt meshkernel::AddNodeAction::NodeId() const
{
    return m_nodeId;
}

const meshkernel::Point& meshkernel::AddNodeAction::Node() const
{
    return m_node;
}
