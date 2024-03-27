#include "MeshKernel/CurvilinearGrid/CurvilinearGridBlockUndo.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"

std::unique_ptr<meshkernel::CurvilinearGridBlockUndo> meshkernel::CurvilinearGridBlockUndo::Create(CurvilinearGrid& grid,
                                                                                                   const CurvilinearGridNodeIndices& startOffset,
                                                                                                   const CurvilinearGridNodeIndices& endOffset)
{
    return std::make_unique<CurvilinearGridBlockUndo>(grid, startOffset, endOffset);
}

meshkernel::CurvilinearGridBlockUndo::CurvilinearGridBlockUndo(CurvilinearGrid& grid,
                                                               const CurvilinearGridNodeIndices& startOffset,
                                                               const CurvilinearGridNodeIndices& endOffset)
    : BaseMeshUndoAction<CurvilinearGridBlockUndo, CurvilinearGrid>(grid), m_block(startOffset, endOffset)
{
    m_block.CopyFrom(grid);
}

void meshkernel::CurvilinearGridBlockUndo::Swap(CurvilinearGrid& grid)
{
    m_block.Swap(grid);
}

std::uint64_t meshkernel::CurvilinearGridBlockUndo::MemorySize() const
{
    return static_cast<std::uint64_t>(sizeof(*this)) + m_block.MemorySize();
}

void meshkernel::CurvilinearGridBlockUndo::Print(std::ostream& out) const
{
    out << "CurvilinearGridBlockUndo: state " << to_string(State())
        << ", startOffset = {" << m_block.StartOffset().m_n << ", " << m_block.StartOffset().m_m << "}, "
        << "endOffset = {" << m_block.EndOffset().m_n << ", " << m_block.EndOffset().m_m << "}"
        << std::endl;
}
