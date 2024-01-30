#ifndef BASE_TRANSACTION__HPP
#define BASE_TRANSACTION__HPP

#include "Transaction.hpp"

template <typename DerivedTransction, class SimpleMesh>
class BaseMeshTransaction : public Transaction
{
public:
    BaseMeshTransaction(SimpleMesh& mesh) : mesh_(mesh) {}

protected:
    void doCommit();

    void doRestore();
private:
    SimpleMesh& mesh_;
};

template <typename DerivedTransction, class SimpleMesh>
void BaseMeshTransaction<DerivedTransction, SimpleMesh>::doCommit()
{
    mesh_.commit(static_cast<DerivedTransction&>(*this));
}

template <typename DerivedTransction, class SimpleMesh>
void BaseMeshTransaction<DerivedTransction, SimpleMesh>::doRestore()
{
    mesh_.restore(static_cast<DerivedTransction&>(*this));
}

#endif // BASE_TRANSACTION__HPP
