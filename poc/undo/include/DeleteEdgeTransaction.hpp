#ifndef DELETE_NODE_TRANSACTION__HPP
#define DELETE_NODE_TRANSACTION__HPP

#include <memory>

#include "BaseMeshTransaction.hpp"
#include "MeshTypes.hpp"

class SimpleMesh;

class DeleteEdgeTransaction : public BaseMeshTransaction<DeleteEdgeTransaction, SimpleMesh>
{
public:
    static TransactionPtr create(SimpleMesh& mesh, const int id, const int start, const int end)
    {
        return std::make_unique<DeleteEdgeTransaction>(mesh, id, start, end);
    }

    DeleteEdgeTransaction(SimpleMesh& mesh, const int id, const int start, const int end) : BaseMeshTransaction<DeleteEdgeTransaction, SimpleMesh>(mesh), edgeId_(id), start_(start), end_(end) {}

    int edgeId() const
    {
        return edgeId_;
    }

    int start() const
    {
        return start_;
    }

    int end() const
    {
        return end_;
    }

private:
    int edgeId_;
    int start_;
    int end_;
};

#endif // DELETE_NODE_TRANSACTION__HPP
