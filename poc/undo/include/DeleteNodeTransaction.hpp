#ifndef DELETE_NODE_TRANSACTION__HPP
#define DELETE_NODE_TRANSACTION__HPP

#include <memory>
#include <utility>
#include <vector>

#include "BaseMeshTransaction.hpp"
#include "MeshTypes.hpp"

class SimpleMesh;

class DeleteNodeTransaction : public BaseMeshTransaction<DeleteNodeTransaction, SimpleMesh>
{
public:
    static TransactionPtr create(SimpleMesh& mesh, const size_t id, const Point& point)
    {
        return std::make_unique<DeleteNodeTransaction>(mesh, id, point);
    }

    DeleteNodeTransaction(SimpleMesh& mesh, const size_t id, const Point& p) : BaseMeshTransaction<DeleteNodeTransaction, SimpleMesh>(mesh), nodeId_(id), point_(p) {}

    // Should we have a pointer to the base class here
    // or a pointer to the DeleteNodeTransaction?
    // If the latter, then perhaps better way of handling
    // edge deletion is available
    void emplace_back(TransactionPtr&& transaction)
    {
        deletedEdges.emplace_back(std::move(transaction));
    }

    size_t nodeId() const
    {
        return nodeId_;
    }

    const Point& point() const
    {
        return point_;
    }

private:
    void doCommit()
    {
        // Where is better to commit the edge deletions?
        // Here or in the mesh? Same for restoring edge deletions.
        for (auto& transaction : deletedEdges)
        {
            transaction->commit();
        }

        BaseMeshTransaction<DeleteNodeTransaction, SimpleMesh>::doCommit();
    }

    void doRestore()
    {
        BaseMeshTransaction<DeleteNodeTransaction, SimpleMesh>::doRestore();

        for (auto iter = deletedEdges.rbegin(); iter != deletedEdges.rend(); ++iter)
        {
            (*iter)->restore();
        }
    }

    // DO we need pre- and post-condition transactions?
    std::vector<TransactionPtr> deletedEdges;
    // std::vector<DeleteEdgeTransaction> deletedEdges;

    size_t nodeId_;
    Point point_;
};

#endif // DELETE_NODE_TRANSACTION__HPP