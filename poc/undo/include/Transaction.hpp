#ifndef TRANSACTION__HPP
#define TRANSACTION__HPP

#include <memory>
#include <vector>

// #define NULL_TRANSACTION 1

// All transactions must be created in the committed state
class Transaction
{
public:
    enum State
    {
        Committed,
        Restored
    };

    Transaction() = default;

    virtual ~Transaction() = default;

    // (re)Apply the changes
    void commit();

    // Undo the changes
    void restore();

    State getState() const;

private:
    virtual void doCommit() = 0;

    virtual void doRestore() = 0;

    // How best to set the initial state?
    // Here or pass as parameter in constructor
    State state = Committed;
};

using TransactionPtr = std::unique_ptr<Transaction>;

#endif // TRANSACTION__HPP
