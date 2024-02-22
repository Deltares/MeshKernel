#include "MeshKernel/MoveNodeAction.hpp"
#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Mesh.hpp"

std::unique_ptr<meshkernel::MoveNodeAction> meshkernel::MoveNodeAction::Create(Mesh& mesh)
{
    return std::make_unique<MoveNodeAction>(mesh);
}

meshkernel::MoveNodeAction::MoveNodeAction(Mesh& mesh) : BaseMeshUndoAction<MoveNodeAction, Mesh>(mesh) {}

void meshkernel::MoveNodeAction::AddDisplacement(const UInt nodeId, const double xDisplacement, const double yDisplacement)
{
    if (nodeId != constants::missing::uintValue)
    {
        m_displacements.emplace_back(NodeDisplacement(Vector(xDisplacement, yDisplacement), nodeId));
    }
}

std::uint64_t meshkernel::MoveNodeAction::MemorySize() const
{
    std::uint64_t size = sizeof(*this) + m_displacements.size() * sizeof(NodeDisplacement);
    return size;
}

void meshkernel::MoveNodeAction::Print(std::ostream& out) const
{
    out << "MoveNodeAction: state " << to_string(State()) << std::endl;

    for (const NodeDisplacement& displacement : m_displacements)
    {
        out << "{ node => " << displacement.m_nodeId
            << ", displacement => { " << displacement.m_displacement.x() << ", "
            << displacement.m_displacement.y()
            << "}} ";
    }

    out << std::endl;
}
