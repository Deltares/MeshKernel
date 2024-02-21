#include "MeshKernel/AddEdgeAction.hpp"
#include "MeshKernel/Mesh.hpp"


std::unique_ptr<meshkernel::AddEdgeAction> meshkernel::AddEdgeAction::Create(Mesh& mesh, const UInt id, const UInt start, const UInt end)
{
    return std::make_unique<AddEdgeAction>(mesh, id, start, end);
}

meshkernel::AddEdgeAction::AddEdgeAction(Mesh& mesh, const UInt id, const UInt start, const UInt end) : BaseMeshUndoAction<AddEdgeAction, Mesh>(mesh), m_edgeId(id), m_edge(start, end) {}

meshkernel::UInt meshkernel::AddEdgeAction::EdgeId() const
{
    return m_edgeId;
}

const meshkernel::Edge& meshkernel::AddEdgeAction::GetEdge() const
{
    return m_edge;
}
