#include "MeshKernel/CurvilinearGrid/CurvilinearGridBlockUndoAction.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"

std::unique_ptr<meshkernel::CurvilinearGridBlockUndoAction> meshkernel::CurvilinearGridBlockUndoAction::Create(CurvilinearGrid& grid,
                                                                                                               const CurvilinearGridNodeIndices& startOffset,
                                                                                                               const CurvilinearGridNodeIndices& endOffset)
{
    return std::make_unique<CurvilinearGridBlockUndoAction>(grid, startOffset, endOffset);
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

void meshkernel::CurvilinearGridBlockUndoAction::Print(std::ostream& out) const
{
    out << "CurvilinearGridBlockUndoAction: state " << to_string(State())
        << ", startOffset = {" << m_block.StartOffset().m_n << ", " << m_block.StartOffset().m_m << "}, "
        << "endOffset = {" << m_block.EndOffset().m_n << ", " << m_block.EndOffset().m_m << "}"
        << std::endl;
}
