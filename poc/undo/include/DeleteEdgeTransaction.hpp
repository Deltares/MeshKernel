#ifndef DELETE_EDGE_TRANSACTION__HPP
#define DELETE_EDGE_TRANSACTION__HPP

#include <memory>

#include "BaseMeshTransaction.hpp"
#include "MeshTypes.hpp"

class SimpleMesh;

class DeleteEdgeTransaction : public BaseMeshTransaction<DeleteEdgeTransaction, SimpleMesh>
{
public:
    static TransactionPtr create(SimpleMesh& mesh, const size_t id, const size_t start, const size_t end)
    {
        return std::make_unique<DeleteEdgeTransaction>(mesh, id, start, end);
    }

    DeleteEdgeTransaction(SimpleMesh& mesh, const size_t id, const size_t start, const size_t end) : BaseMeshTransaction<DeleteEdgeTransaction, SimpleMesh>(mesh), edgeId_(id), start_(start), end_(end) {}

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

#endif // DELETE_EDGE_TRANSACTION__HPP
