#include "MeshKernel/UndoActions/ResetEdgeAction.hpp"
#include "MeshKernel/Mesh.hpp"

std::unique_ptr<meshkernel::ResetEdgeAction> meshkernel::ResetEdgeAction::Create(Mesh& mesh, const UInt id, const Edge& initial, const Edge& updated)
{
    return std::make_unique<ResetEdgeAction>(mesh, id, initial, updated);
}

meshkernel::ResetEdgeAction::ResetEdgeAction(Mesh& mesh, const UInt id, const Edge& initial, const Edge& updated) : BaseMeshUndoAction<ResetEdgeAction, Mesh>(mesh), m_edgeId(id), m_initialEdge(initial), m_updatedEdge(updated) {}

meshkernel::UInt meshkernel::ResetEdgeAction::EdgeId() const
{
    return m_edgeId;
}

const meshkernel::Edge& meshkernel::ResetEdgeAction::InitialEdge() const
{
    return m_initialEdge;
}

const meshkernel::Edge& meshkernel::ResetEdgeAction::UpdatedEdge() const
{
    return m_updatedEdge;
}

void meshkernel::ResetEdgeAction::Print(std::ostream& out) const
{
    out << "ResetEdgeAction: state " << to_string(State()) << ", edgeId = " << m_edgeId
        << " initial = {" << m_initialEdge.first << ", " << m_initialEdge.second << "}, "
        << " updated = {" << m_updatedEdge.first << ", " << m_updatedEdge.second << "}"
        << std::endl;
}
