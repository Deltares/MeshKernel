#ifndef ADD_NODE_TRANSACTION__HPP
#define ADD_NODE_TRANSACTION__HPP

#include <memory>

#include "BaseMeshTransaction.hpp"
#include "MeshTypes.hpp"

class SimpleMesh;

class AddNodeTransaction : public BaseMeshTransaction<AddNodeTransaction, SimpleMesh>
{
public:
    static TransactionPtr create(SimpleMesh& mesh, const size_t id, const Point& point)
    {
        return std::make_unique<AddNodeTransaction>(mesh, id, point);
    }

    AddNodeTransaction(SimpleMesh& mesh, const size_t id, const Point& p) : BaseMeshTransaction<AddNodeTransaction, SimpleMesh>(mesh), nodeId_(id), node_(p) {}

    size_t nodeId() const
    {
        return nodeId_;
    }

    const Point& point() const
    {
        return node_;
    }

private:
    size_t nodeId_;
    Point node_;
};

#endif // ADD_NODE_TRANSACTION__HPP
