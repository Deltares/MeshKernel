#include "MeshKernel/OrthogonalizationAndSmoothingAction.hpp"
#include "MeshKernel/Mesh2D.hpp"


std::unique_ptr<meshkernel::OrthogonalizationAndSmoothingAction> meshkernel::OrthogonalizationAndSmoothingAction::Create(Mesh2D& mesh, const std::vector<Point>& nodes)
{
    return std::make_unique<OrthogonalizationAndSmoothingAction> (mesh, nodes);
}

meshkernel::OrthogonalizationAndSmoothingAction::OrthogonalizationAndSmoothingAction(Mesh2D& mesh, const std::vector<Point>& nodes)
    : BaseMeshUndoAction<OrthogonalizationAndSmoothingAction, Mesh2D> (mesh), m_nodes (nodes){}

void meshkernel::OrthogonalizationAndSmoothingAction::Swap(std::vector<Point>& nodes)
{
    std::swap_ranges(m_nodes.begin (), m_nodes.end (), nodes.begin ());
}

std::uint64_t meshkernel::OrthogonalizationAndSmoothingAction::MemorySize() const
{
    // Change with implementation on other computer.
    return 0;
}

void meshkernel::OrthogonalizationAndSmoothingAction::Print(std::ostream& out [[maybe_unused]]) const
{
    // Change with implementation on other computer.
}
