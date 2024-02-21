#include "MeshKernel/DeleteEdgeAction.hpp"
#include "MeshKernel/Mesh.hpp"

std::unique_ptr<meshkernel::DeleteEdgeAction> meshkernel::DeleteEdgeAction::Create(Mesh& mesh, const UInt id, const UInt start, const UInt end)
{
    return std::make_unique<DeleteEdgeAction>(mesh, id, start, end);
}

meshkernel::DeleteEdgeAction::DeleteEdgeAction(Mesh& mesh, const UInt id, const UInt start, const UInt end) : BaseMeshUndoAction<DeleteEdgeAction, Mesh>(mesh), m_edgeId(id), m_edge(start, end) {}

meshkernel::UInt meshkernel::DeleteEdgeAction::EdgeId() const
{
    return m_edgeId;
}

const meshkernel::Edge& meshkernel::DeleteEdgeAction::GetEdge() const
{
    return m_edge;
}

void meshkernel::DeleteEdgeAction::Print(std::ostream& out) const
{
    out << "DeleteEdgeAction: state " << to_string(State()) << ", edgeId = " << m_edgeId << " edge = {" << m_edge.first << ", " << m_edge.second << "}" << std::endl;
}
