#include "MeshKernel/UndoActions/UndoActionStack.hpp"
#include "MeshKernel/Exceptions.hpp"

#include <algorithm>
#include <utility>

const meshkernel::UInt meshkernel::UndoActionStack::MaxUndoSize = 50;

void meshkernel::UndoActionStack::Add(UndoActionPtr&& action, const int actionId, const std::string& info)
{
    if (action != nullptr)
    {
        if (action->GetState() != UndoAction::State::Committed)
        {
            throw ConstraintError("Cannot add an action in the {} state.", UndoAction::to_string(action->GetState()));
        }

        m_committed.emplace_back(std::move(action), actionId, info);

        // Clear the restored actions for this actionId.
        // Adding a new undo action for an actionId, means that no action for this Id can be restored
        m_restored.remove_if([actionId](const UndoActionForMesh& undoAction)
                             { return undoAction.m_actionId == actionId; });

        if (m_committed.size() > MaxUndoSize)
        {
            // If the number of undo-actions is greater than the maximum, then remove the first item in the list.
            m_committed.pop_front();
        }
    }
    else
    {
        // Logging message
    }
}

void meshkernel::UndoActionStack::print() const
{

    std::cout << "stack info:" << std::endl;
    std::cout << "committed: " << m_committed.size() << std::endl;

    size_t count = 0;

    for (auto iter = m_committed.begin(); iter != m_committed.end(); ++iter)
    {
        std::cout << "undo " << count << "  " << iter->m_actionId << "  " << iter->m_info << std::endl;
        ++count;
    }

    std::cout << "restored: " << m_restored.size() << std::endl;

    count = 0;

    for (auto iter = m_restored.begin(); iter != m_restored.end(); ++iter)
    {
        std::cout << "redo " << count << "  " << iter->m_actionId << "  " << iter->m_info << std::endl;
        ++count;
    }
}

meshkernel::UInt meshkernel::UndoActionStack::Size() const
{
    return static_cast<UInt>(m_committed.size() + m_restored.size());
}

meshkernel::UInt meshkernel::UndoActionStack::CommittedSize(const int actionId) const
{
    if (actionId == constants::missing::intValue)
    {
        return static_cast<UInt>(m_committed.size());
    }
    else
    {
        return static_cast<UInt>(std::ranges::count_if(m_committed, [actionId](const UndoActionForMesh& undoAction)
                                                       { return undoAction.m_actionId == actionId; }));
    }
}

meshkernel::UInt meshkernel::UndoActionStack::RestoredSize(const int actionId) const
{
    if (actionId == constants::missing::intValue)
    {
        return static_cast<UInt>(m_restored.size());
    }
    else
    {
        return static_cast<UInt>(std::ranges::count_if(m_restored, [actionId](const UndoActionForMesh& undoAction)
                                                       { return undoAction.m_actionId == actionId; }));
    }
}

std::tuple<bool, int> meshkernel::UndoActionStack::Undo()
{
    bool didUndo = false;
    int actionId = constants::missing::intValue;

    if (!m_committed.empty())
    {
        // Perform undo operation
        m_committed.back().m_undoAction->Restore();
        actionId = m_committed.back().m_actionId;
        // Now move to restored stack
        m_restored.emplace_back(std::move(m_committed.back()));
        m_committed.pop_back();
        didUndo = true;
    }

    return {didUndo, actionId};
}

std::tuple<bool, int> meshkernel::UndoActionStack::Commit()
{
    bool didCommit = false;
    int actionId = constants::missing::intValue;

    if (!m_restored.empty())
    {
        // Perform commit (redo) operation
        m_restored.back().m_undoAction->Commit();
        actionId = m_restored.back().m_actionId;
        // Now move to committed stack
        m_committed.emplace_back(std::move(m_restored.back()));
        m_restored.pop_back();
        didCommit = true;
    }

    return {didCommit, actionId};
}

void meshkernel::UndoActionStack::Clear()
{
    m_committed.clear();
    m_restored.clear();
}

meshkernel::UInt meshkernel::UndoActionStack::Remove(const int actionId)
{
    UInt removedCount = 0;

    auto hasMatchingActionId = [actionId](const UndoActionForMesh& action)
    { return action.m_actionId == actionId; };

    removedCount = static_cast<UInt>(m_committed.remove_if(hasMatchingActionId));
    removedCount += static_cast<UInt>(m_restored.remove_if(hasMatchingActionId));

    return removedCount;
}

std::uint64_t meshkernel::UndoActionStack::MemorySize() const
{
    const std::uint64_t committedSize = std::accumulate(m_committed.begin(), m_committed.end(), static_cast<std::uint64_t>(0u),
                                                        [](const std::uint64_t partialSum, const UndoActionForMesh& action)
                                                        { return partialSum + action.m_undoAction->MemorySize(); });
    const std::uint64_t restoredSize = std::accumulate(m_restored.begin(), m_restored.end(), static_cast<std::uint64_t>(0u),
                                                       [](const std::uint64_t partialSum, const UndoActionForMesh& action)
                                                       { return partialSum + action.m_undoAction->MemorySize(); });

    return sizeof(*this) + committedSize + restoredSize;
}
