#ifndef SIMPLE_MESH__HPP
#define SIMPLE_MESH__HPP

#include <tuple>
#include <vector>

#include "AddEdgeTransaction.hpp"
#include "AddNodeTransaction.hpp"
#include "BaseMeshTransaction.hpp"
#include "DeleteEdgeTransaction.hpp"
#include "DeleteNodeTransaction.hpp"
#include "MeshTypes.hpp"
#include "Transaction.hpp"

class SimpleMesh
{
public:
    SimpleMesh(const std::vector<Point>& nodes,
               const std::vector<Edge>& edges);

    const Point& getNode(const int id) const;

    const Edge& getEdge(const int id) const;

    const std::vector<int>& getNodeConnectivity(const int id) const;

    // Could return pointer to AddNodeTransaction instead of the base class.
    std::tuple<int, TransactionPtr> addNode(const Point& p);

    void commit(AddNodeTransaction& transaction);

    void restore(AddNodeTransaction& transaction);

    void resetNode(const int id, const Point& p);

    TransactionPtr deleteNode(const int id);

    void commit(DeleteNodeTransaction& transaction);

    void restore(DeleteNodeTransaction& transaction);

    std::tuple<int, TransactionPtr> addEdge(const int start, const int end);

    void resetEdge(const int id, const int start, const int end);

    void commit(AddEdgeTransaction& transaction);

    void restore(AddEdgeTransaction& transaction);

    TransactionPtr deleteEdge(const int id);

    void commit(DeleteEdgeTransaction& transaction);

    void restore(DeleteEdgeTransaction& transaction);

    int getNumberOfNodes() const;

    int getNumberOfEdges() const;

    int getNumberOfValidNodes() const;

    int getNumberOfValidEdges() const;

    void administrate();

    void print() const;

private:
    std::vector<Point> nodes_;
    std::vector<Edge> edges_;
    std::vector<std::vector<int>> nodeEdges_;
};

//--------------------------------

#endif // SIMPLE_MESH__HPP
