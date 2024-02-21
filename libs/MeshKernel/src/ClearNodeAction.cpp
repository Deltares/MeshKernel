#include "MeshKernel/ClearNodeAction.hpp"
#include "MeshKernel/Mesh.hpp"

std::unique_ptr<meshkernel::ClearNodeAction> meshkernel::ClearNodeAction::Create(Mesh& mesh, const UInt id, const Point& node)
{
    return std::make_unique<ClearNodeAction>(mesh, id, node);
}

meshkernel::ClearNodeAction::ClearNodeAction(Mesh& mesh, const UInt id, const Point& node) : BaseMeshUndoAction<ClearNodeAction, Mesh>(mesh), m_nodeId(id), m_node(node) {}

meshkernel::UInt meshkernel::ClearNodeAction::NodeId() const
{
    return m_nodeId;
}

const meshkernel::Point& meshkernel::ClearNodeAction::Node() const
{
    return m_node;
}

void meshkernel::ClearNodeAction::Print(std::ostream& out) const
{
    out << "ClearNodeAction: state " << to_string(State()) << ", nodeId = " << m_nodeId << " node = {" << m_node.x << ", " << m_node.y << "}" << std::endl;
}
