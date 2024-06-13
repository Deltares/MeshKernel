#include "MeshKernel/UndoActions/UndoAction.hpp"

const std::string& meshkernel::UndoAction::to_string(const State state)
{
    static std::string committedStr = "Committed";
    static std::string restoredStr = "Restored";
    static std::string unknownStr = "UNKNOWN";

    switch (state)
    {
    case State::Committed:
        return committedStr;
    case State::Restored:
        return restoredStr;
    default:
        return unknownStr;
    }
}

meshkernel::UndoAction::State meshkernel::UndoAction::GetState() const
{
    return m_state;
}

void meshkernel::UndoAction::Commit()
{
    if (m_state == State::Restored)
    {
        DoCommit();
        m_state = State::Committed;
    }
}

void meshkernel::UndoAction::Restore()
{
    if (m_state == State::Committed)
    {
        DoRestore();
        m_state = State::Restored;
    }
}

std::uint64_t meshkernel::UndoAction::MemorySize() const
{
    return sizeof(*this);
}
