#include "MeshKernel/CurvilinearGrid/CurvilinearGridRefinementUndoAction.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"

std::unique_ptr<meshkernel::CurvilinearGridRefinementUndoAction> meshkernel::CurvilinearGridRefinementUndoAction::Create(CurvilinearGrid& grid)
{
    return std::make_unique<CurvilinearGridRefinementUndoAction>(grid);
}

meshkernel::CurvilinearGridRefinementUndoAction::CurvilinearGridRefinementUndoAction(CurvilinearGrid& grid)
    : BaseMeshUndoAction<CurvilinearGridRefinementUndoAction, CurvilinearGrid>(grid), m_nodes(grid.GetNodes()), m_lower(grid.StartOffset()), m_upper(grid.EndOffset()) {}

void meshkernel::CurvilinearGridRefinementUndoAction::Swap(lin_alg::Matrix<Point>& nodes, CurvilinearGridNodeIndices& lower, CurvilinearGridNodeIndices& upper)
{
    lin_alg::Matrix<Point> temp = nodes;
    nodes = m_nodes;
    m_nodes = temp;
    std::swap(m_lower, lower);
    std::swap(m_upper, upper);
    // m_nodes.swap(nodes);
}

std::uint64_t meshkernel::CurvilinearGridRefinementUndoAction::MemorySize() const
{
    std::uint64_t result = 0;
    result += sizeof(*this);
    result += static_cast<std::uint64_t>(m_nodes.rows() * m_nodes.cols() * sizeof(Point));
    return result;
}
