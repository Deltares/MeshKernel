#include "MeshKernel/UndoAction.hpp"

void meshkernel::UndoAction::Commit()
{
    if (m_state == Restored)
    {
        DoCommit();
        m_state = Committed;
    }
}

void meshkernel::UndoAction::Restore()
{
    if (m_state == Committed)
    {
        DoRestore();
        m_state = Restored;
    }
}
