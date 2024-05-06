#include "MeshKernel/UndoActions/DeleteNodeAction.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Formatting.hpp"
#include "MeshKernel/Mesh.hpp"

#include <ranges>
#include <utility>

std::unique_ptr<meshkernel::DeleteNodeAction> meshkernel::DeleteNodeAction::Create(Mesh& mesh, const UInt id, const Point& node)
{
    return std::make_unique<DeleteNodeAction>(mesh, id, node);
}

meshkernel::DeleteNodeAction::DeleteNodeAction(Mesh& mesh, const UInt id, const Point& node) : m_mesh(mesh), m_nodeId(id), m_node(node) {}

void meshkernel::DeleteNodeAction::Add(std::unique_ptr<DeleteEdgeAction>&& action)
{
    if (action != nullptr)
    {
        if (action->GetState() == UndoAction::State::Restored)
        {
            throw ConstraintError("Cannot add an action in the {} state.", UndoAction::to_string(action->GetState()));
        }

        m_deletedEdges.emplace_back(std::move(action));
    }
}

meshkernel::UInt meshkernel::DeleteNodeAction::NodeId() const
{
    return m_nodeId;
}

const meshkernel::Point& meshkernel::DeleteNodeAction::Node() const
{
    return m_node;
}

void meshkernel::DeleteNodeAction::DoCommit()
{
    // This does seem to be exposing some information about the implementation
    // having to delete the edges first then the node, same for restore
    for (const std::unique_ptr<DeleteEdgeAction>& action : m_deletedEdges)
    {
        action->Commit();
    }

    m_mesh.CommitAction(*this);
}

void meshkernel::DeleteNodeAction::DoRestore()
{
    m_mesh.RestoreAction(*this);

    for (const std::unique_ptr<DeleteEdgeAction>& action : m_deletedEdges | std::views::reverse)
    {
        action->Restore();
    }
}

std::uint64_t meshkernel::DeleteNodeAction::MemorySize() const
{
    std::uint64_t size = sizeof(*this) + m_deletedEdges.capacity() * sizeof(std::unique_ptr<DeleteEdgeAction>);

    for (const std::unique_ptr<DeleteEdgeAction>& action : m_deletedEdges)
    {
        size += action->MemorySize();
    }

    return size;
}
