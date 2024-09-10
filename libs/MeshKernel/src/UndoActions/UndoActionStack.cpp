#include "MeshKernel/UndoActions/UndoActionStack.hpp"
#include "MeshKernel/Exceptions.hpp"

#include <algorithm>
#include <utility>

const meshkernel::UInt meshkernel::UndoActionStack::DefaultMaxUndoSize = 50;

meshkernel::UndoActionStack::UndoActionStack(const UInt maximumSize) : m_maxUndoSize(maximumSize) {}

void meshkernel::UndoActionStack::Add(UndoActionPtr&& action, const int actionId)
{
    if (m_maxUndoSize > 0 && action != nullptr)
    {
        if (action->GetState() != UndoAction::State::Committed)
        {
            throw ConstraintError("Cannot add an action in the {} state.", UndoAction::to_string(action->GetState()));
        }

        m_committed.emplace_back(std::move(action), actionId);

        // Clear the restored actions for this actionId.
        // Adding a new undo action for an actionId, means that no action for this Id can be restored
        m_restored.remove_if([actionId](const UndoActionForMesh& undoAction)
                             { return undoAction.m_actionId == actionId; });

        if (m_committed.size() > m_maxUndoSize)
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

void meshkernel::UndoActionStack::SetMaximumSize(const UInt maximumSize)
{
    if (maximumSize == 0)
    {
        m_committed.clear();
        m_restored.clear();
    }
    else if (maximumSize < m_committed.size())
    {
        UInt loopLimit = static_cast<UInt>(m_committed.size());

        for (UInt i = maximumSize; i < loopLimit; ++i)
        {
            m_committed.pop_front();
        }
    }

    m_maxUndoSize = maximumSize;
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

std::optional<int> meshkernel::UndoActionStack::Undo()
{

    if (!m_committed.empty())
    {
        // Perform undo operation
        m_committed.back().m_undoAction->Restore();
        int actionId = m_committed.back().m_actionId;
        // Now move to restored stack
        m_restored.emplace_back(std::move(m_committed.back()));
        m_committed.pop_back();
        return actionId;
    }

    return std::nullopt;
}

std::optional<int> meshkernel::UndoActionStack::Commit()
{

    if (!m_restored.empty())
    {
        // Perform commit (redo) operation
        m_restored.back().m_undoAction->Commit();
        int actionId = m_restored.back().m_actionId;
        // Now move to committed stack
        m_committed.emplace_back(std::move(m_restored.back()));
        m_restored.pop_back();

        if (m_committed.size() > m_maxUndoSize)
        {
            // If the number of undo-actions is greater than the maximum, then remove the first item in the list.
            m_committed.pop_front();
        }

        return actionId;
    }

    return std::nullopt;
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
