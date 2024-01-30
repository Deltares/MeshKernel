#ifndef ADD_EDGE_TRANSACTION__HPP
#define ADD_EDGE_TRANSACTION__HPP

#include <memory>

#include "BaseMeshTransaction.hpp"
#include "MeshTypes.hpp"

class SimpleMesh;

class AddEdgeTransaction : public BaseMeshTransaction<AddEdgeTransaction, SimpleMesh>
{
public:
    static TransactionPtr create(SimpleMesh& mesh, const int id, const int start, const int end)
    {
        return std::make_unique<AddEdgeTransaction>(mesh, id, start, end);
    }

    AddEdgeTransaction(SimpleMesh& mesh, const int id, const int start, const int end) : BaseMeshTransaction<AddEdgeTransaction, SimpleMesh>(mesh), edgeId_(id), start_(start), end_(end) {}

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

#endif // ADD_EDGE_TRANSACTION__HPP
