#include "MeshKernel/CurvilinearGrid/AddGridLineUndoAction.hpp"
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

void meshkernel::AddGridLineUndoAction::Print(std::ostream& out) const
{
    out << "AddGridLineUndoAction: state " << to_string(State())
        << ", startOffset = {" << m_startOffset.m_n << ", " << m_startOffset.m_m << "}, "
        << "endOffset = {" << m_endOffset.m_n << ", " << m_endOffset.m_m << "}"
        << std::endl;
}
