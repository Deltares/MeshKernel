#include "MeshKernel/UndoAction.hpp"

std::string meshkernel::UndoAction::to_string(const ActionState state)
{
    static std::string committedStr = "Committed";
    static std::string restoredStr = "Restored";
    static std::string unknownStr = "UNKNOWN";

    switch (state)
    {
    case Committed:
        return committedStr;
    case Restored:
        return restoredStr;
    default:
        return unknownStr;
    }
}

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

void meshkernel::UndoAction::Print(std::ostream& out [[maybe_unused]]) const
{
    // do nothing for the default
}

std::uint64_t meshkernel::UndoAction::MemorySize() const
{
    return sizeof(*this);
}
