#include "MeshKernel/UndoActionStack.hpp"

void meshkernel::UndoActionStack::Add(UndoActionPtr&& action)
{
    if (action != nullptr)
    {
        m_committed.emplace_back(std::move(action));
        // TODO Do I need to do this?
        m_restored.reserve(m_committed.size());
    }
    else
    {
        // Logging message
    }
}

bool meshkernel::UndoActionStack::Undo()
{
    bool didUndo = false;

    if (m_committed.size() > 0)
    {
        m_committed.back()->Restore();
        m_restored.emplace_back(std::move(m_committed.back()));
        m_committed.pop_back();
        didUndo = true;
    }

    return didUndo;
}

bool meshkernel::UndoActionStack::Commit()
{
    bool didCommit = false;

    if (m_restored.size() > 0)
    {
        m_restored.back()->Commit();
        m_committed.emplace_back(std::move(m_restored.back()));
        m_restored.pop_back();
        didCommit = true;
    }

    return didCommit;
}
