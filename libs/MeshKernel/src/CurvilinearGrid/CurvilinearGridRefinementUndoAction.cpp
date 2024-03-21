#include "MeshKernel/CurvilinearGrid/CurvilinearGridRefinementUndoAction.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"

std::unique_ptr<meshkernel::CurvilinearGridRefinementUndoAction> meshkernel::CurvilinearGridRefinementUndoAction::Create(CurvilinearGrid& grid)
{
    return std::make_unique<CurvilinearGridRefinementUndoAction>(grid);
}

meshkernel::CurvilinearGridRefinementUndoAction::CurvilinearGridRefinementUndoAction(CurvilinearGrid& grid)
    : BaseMeshUndoAction<CurvilinearGridRefinementUndoAction, CurvilinearGrid>(grid), m_nodes(grid.GetNodes()) {}

void meshkernel::CurvilinearGridRefinementUndoAction::Swap(lin_alg::Matrix<Point>& nodes)
{
    lin_alg::Matrix<Point> temp = nodes;
    nodes = m_nodes;
    m_nodes = temp;
    // m_nodes.swap(nodes);
}
