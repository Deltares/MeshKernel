#ifndef ADD_EDGE_TRANSACTION__HPP
#define ADD_EDGE_TRANSACTION__HPP

#include <memory>

#include "BaseMeshTransaction.hpp"
#include "MeshTypes.hpp"

class SimpleMesh;

class AddEdgeTransaction : public BaseMeshTransaction<AddEdgeTransaction, SimpleMesh>
{
public:
    static TransactionPtr create(SimpleMesh& mesh, const size_t id, const size_t start, const size_t end)
    {
        return std::make_unique<AddEdgeTransaction>(mesh, id, start, end);
    }

    AddEdgeTransaction(SimpleMesh& mesh, const size_t id, const size_t start, const size_t end) : BaseMeshTransaction<AddEdgeTransaction, SimpleMesh>(mesh), edgeId_(id), start_(start), end_(end) {}

    size_t edgeId() const
    {
        return edgeId_;
    }

    size_t start() const
    {
        return start_;
    }

    size_t end() const
    {
        return end_;
    }

private:
    size_t edgeId_;
    size_t start_;
    size_t end_;
};

#endif // ADD_EDGE_TRANSACTION__HPP
