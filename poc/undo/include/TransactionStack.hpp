#ifndef TRANSACTION_STACK__HPP
#define TRANSACTION_STACK__HPP

#include <utility>
#include <vector>

#include "Transaction.hpp"

class TransactionStack
{
public:
    // When adding a new transaction, should probably delete all restored transactions
    // When adding new transactions, could check the size of the committed list and remove transactions more than some number ago
    // e.g. keep the undo list no longer than 10
    // Probably should remove this function if using unique_ptr
    void push_back(TransactionPtr& transaction)
    {
        committed.push_back(std::move(transaction));
        restored.reserve(committed.size());
    }

    // When adding a new transaction, should probably delete all restored transactions
    void emplace_back(TransactionPtr&& transaction)
    {
        committed.emplace_back(std::move(transaction));
        restored.reserve(committed.size());
    }

    bool undo()
    {
        bool didUndo = false;

#ifdef ADD_LOGGING
        std::cout << " % undo " << committed.size() << "  " << restored.size() << std::endl;
#endif

        if (committed.size() > 0)
        {
            committed.back()->restore();
            restored.emplace_back(std::move(committed.back()));
            committed.pop_back();
            didUndo = true;
        }

        return didUndo;
    }

    bool commit()
    {
        bool didCommit = false;

#ifdef ADD_LOGGING
        std::cout << " % restore " << committed.size() << "  " << restored.size() << std::endl;
#endif

        if (restored.size() > 0)
        {
            restored.back()->commit();
            committed.emplace_back(std::move(restored.back()));
            restored.pop_back();
            didCommit = true;
        }

        return didCommit;
    }

    void commitAll()
    {
        while (commit())
        {
            // Nothing else to do
        }
    }

    void restoreAll()
    {
        while (undo())
        {
            // Nothing else to do
        }
    }

private:
    std::vector<TransactionPtr> committed;
    std::vector<TransactionPtr> restored;
};

#endif // TRANSACTION_STACK__HPP
