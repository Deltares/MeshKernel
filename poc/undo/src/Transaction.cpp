#include "Transaction.hpp"

void Transaction::commit()
{
    if (state == Restored)
    {
        doCommit();
        state = Committed;
    }
}

void Transaction::restore()
{
    if (state == Committed)
    {
        doRestore();
        state = Restored;
    }
}

Transaction::State Transaction::getState() const
{
    return state;
}
