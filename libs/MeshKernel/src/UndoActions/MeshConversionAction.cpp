#include "MeshKernel/UndoActions/MeshConversionAction.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Mesh.hpp"

#include <algorithm>

std::unique_ptr<meshkernel::MeshConversionAction> meshkernel::MeshConversionAction::Create(Mesh& mesh)
{
    return std::make_unique<MeshConversionAction>(mesh);
}

meshkernel::MeshConversionAction::MeshConversionAction(Mesh& mesh)
    : BaseMeshUndoAction<MeshConversionAction, Mesh>(mesh), m_nodes(mesh.Nodes()), m_projection(mesh.m_projection) {}

void meshkernel::MeshConversionAction::Swap(std::vector<Point>& nodes, Projection& projection)
{
    if (nodes.size() < m_nodes.size())
    {
        throw ConstraintError("Number of nodes passed is less than nodes stored. {} < {}",
                              nodes.size(), m_nodes.size());
    }

    std::swap_ranges(m_nodes.begin(), m_nodes.end(), nodes.begin());
    std::swap(m_projection, projection);
}

std::uint64_t meshkernel::MeshConversionAction::MemorySize() const
{
    return sizeof(*this) + m_nodes.capacity() * sizeof(Point);
}
