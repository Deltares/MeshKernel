#include "MeshKernel/CurvilinearGrid/CurvilinearGridRefinementUndoAction.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"

std::unique_ptr<meshkernel::CurvilinearGridRefinementUndoAction> meshkernel::CurvilinearGridRefinementUndoAction::Create(CurvilinearGrid& grid)
{
    return std::make_unique<CurvilinearGridRefinementUndoAction>(grid);
}

meshkernel::CurvilinearGridRefinementUndoAction::CurvilinearGridRefinementUndoAction(CurvilinearGrid& grid)
    : BaseMeshUndoAction<CurvilinearGridRefinementUndoAction, CurvilinearGrid>(grid), m_nodes(grid.GetNodes()), m_lower(grid.m_startOffset), m_upper(grid.m_endOffset) {}

void meshkernel::CurvilinearGridRefinementUndoAction::Swap(lin_alg::Matrix<Point>& nodes, CurvilinearGridNodeIndices& lower, CurvilinearGridNodeIndices& upper)
{
    lin_alg::Matrix<Point> temp = nodes;
    nodes = m_nodes;
    m_nodes = temp;
    std::swap(m_lower, lower);
    std::swap(m_upper, upper);
    // m_nodes.swap(nodes);
}
