#include "MeshKernel/UndoActions/DeleteEdgeAction.hpp"
#include "MeshKernel/Formatting.hpp"
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
    out << fmt_ns::vformat("DeleteEdgeAction: state {}, edgeId {}, edge {{{}, {}}}",
                           fmt_ns::make_format_args(to_string(GetState()), m_edgeId, m_edge.first, m_edge.second));
    out << std::endl;
}
