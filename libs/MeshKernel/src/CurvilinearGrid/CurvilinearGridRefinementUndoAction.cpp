#include "MeshKernel/CurvilinearGrid/CurvilinearGridRefinementUndoAction.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"

std::unique_ptr<meshkernel::CurvilinearGridRefinementUndoAction> meshkernel::CurvilinearGridRefinementUndoAction::Create(CurvilinearGrid& grid)
{
    return std::make_unique<CurvilinearGridRefinementUndoAction>(grid);
}

meshkernel::CurvilinearGridRefinementUndoAction::CurvilinearGridRefinementUndoAction(CurvilinearGrid& grid)
    : BaseMeshUndoAction<CurvilinearGridRefinementUndoAction, CurvilinearGrid>(grid), m_nodes(grid.GetNodes()), m_startOffset(grid.StartOffset()), m_endOffset(grid.EndOffset()) {}

void meshkernel::CurvilinearGridRefinementUndoAction::Swap(lin_alg::Matrix<Point>& nodes, CurvilinearGridNodeIndices& startOffset, CurvilinearGridNodeIndices& endOffset)
{
    std::swap(m_startOffset, startOffset);
    std::swap(m_endOffset, endOffset);
    m_nodes.swap(nodes);
}

std::uint64_t meshkernel::CurvilinearGridRefinementUndoAction::MemorySize() const
{
    std::uint64_t result = 0;
    result += sizeof(*this);
    result += static_cast<std::uint64_t>(m_nodes.rows() * m_nodes.cols() * sizeof(Point));
    return result;
}
