#include "MeshKernel/UndoActionStack.hpp"

meshkernel::UndoActionStack::UndoActionStack()
{
    m_committed.reserve(DefaultReserveSize);
    m_restored.reserve(DefaultReserveSize);
}

void meshkernel::UndoActionStack::Add(UndoActionPtr&& action)
{
    if (action != nullptr)
    {
        m_committed.emplace_back(std::move(action));
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
        // Perform undo operation
        m_committed.back()->Restore();
        // Now move to restored stack
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
        // Perform commit (redo) operation
        m_restored.back()->Commit();
        // Now move to committed stack
        m_committed.emplace_back(std::move(m_restored.back()));
        m_restored.pop_back();
        didCommit = true;
    }

    return didCommit;
}