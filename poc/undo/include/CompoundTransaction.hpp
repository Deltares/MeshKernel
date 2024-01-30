#ifndef COMPOUND_TRANSACTION__HPP
#define COMPOUND_TRANSACTION__HPP

#include "Transaction.hpp"
#include <utility>
#include <vector>

class CompoundTransaction : public Transaction
{
public:
    void push_back(TransactionPtr& transaction)
    {
        transactions.push_back(std::move(transaction));
    }

    void emplace_back(TransactionPtr&& transaction)
    {
        transactions.emplace_back(std::move(transaction));
    }

private:
    void doCommit()
    {
        for (size_t i = 0; i < transactions.size(); ++i)
        {
            transactions[i]->commit();
        }
    }

    void doRestore()
    {
        for (auto iter = transactions.rbegin(); iter != transactions.rend(); ++iter)
        {
            (*iter)->restore();
        }
    }

    std::vector<TransactionPtr> transactions;
};

#endif // COMPOUND_TRANSACTION__HPP
