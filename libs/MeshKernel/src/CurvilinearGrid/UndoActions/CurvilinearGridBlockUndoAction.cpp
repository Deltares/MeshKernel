#include "MeshKernel/CurvilinearGrid/UndoActions/CurvilinearGridBlockUndoAction.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"

std::unique_ptr<meshkernel::CurvilinearGridBlockUndoAction> meshkernel::CurvilinearGridBlockUndoAction::Create(CurvilinearGrid& grid,
                                                                                                               const CurvilinearGridNodeIndices& startOffset,
                                                                                                               const CurvilinearGridNodeIndices& endOffset)
{
    return std::make_unique<CurvilinearGridBlockUndoAction>(grid, startOffset, endOffset);
}

std::unique_ptr<meshkernel::CurvilinearGridBlockUndoAction> meshkernel::CurvilinearGridBlockUndoAction::Create(CurvilinearGrid& grid)
{
    return std::make_unique<CurvilinearGridBlockUndoAction>(grid,
                                                            CurvilinearGridNodeIndices(0, 0),
                                                            CurvilinearGridNodeIndices(grid.NumN() - 1, grid.NumM() - 1));
}

meshkernel::CurvilinearGridBlockUndoAction::CurvilinearGridBlockUndoAction(CurvilinearGrid& grid,
                                                                           const CurvilinearGridNodeIndices& startOffset,
                                                                           const CurvilinearGridNodeIndices& endOffset)
    : BaseMeshUndoAction<CurvilinearGridBlockUndoAction, CurvilinearGrid>(grid), m_block(startOffset, endOffset)
{
    m_block.CopyFrom(grid);
}

void meshkernel::CurvilinearGridBlockUndoAction::Swap(CurvilinearGrid& grid)
{
    m_block.Swap(grid);
}

std::uint64_t meshkernel::CurvilinearGridBlockUndoAction::MemorySize() const
{
    return static_cast<std::uint64_t>(sizeof(*this)) + m_block.MemorySize();
}
