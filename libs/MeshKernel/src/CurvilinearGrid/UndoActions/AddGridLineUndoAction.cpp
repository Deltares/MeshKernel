#include "MeshKernel/CurvilinearGrid/UndoActions/AddGridLineUndoAction.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"

std::unique_ptr<meshkernel::AddGridLineUndoAction> meshkernel::AddGridLineUndoAction::Create(CurvilinearGrid& grid,
                                                                                             const CurvilinearGridNodeIndices& startOffset,
                                                                                             const CurvilinearGridNodeIndices& endOffset)
{
    return std::make_unique<AddGridLineUndoAction>(grid, startOffset, endOffset);
}

meshkernel::AddGridLineUndoAction::AddGridLineUndoAction(CurvilinearGrid& grid,
                                                         const CurvilinearGridNodeIndices& startOffset,
                                                         const CurvilinearGridNodeIndices& endOffset)
    : BaseMeshUndoAction<AddGridLineUndoAction, CurvilinearGrid>(grid),
      m_startOffset(startOffset),
      m_endOffset(endOffset) {}
