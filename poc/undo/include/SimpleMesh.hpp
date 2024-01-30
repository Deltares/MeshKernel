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

    const Point& getNode(const size_t id) const;

    const Edge& getEdge(const size_t id) const;

    const std::vector<size_t>& getNodeConnectivity(const size_t id) const;

    std::tuple<size_t, std::unique_ptr<AddNodeTransaction>> addNode(const Point& p);

    void commit(AddNodeTransaction& transaction);

    void restore(AddNodeTransaction& transaction);

    void resetNode(const size_t id, const Point& p);

    std::unique_ptr<DeleteNodeTransaction> deleteNode(const size_t id);

    void commit(DeleteNodeTransaction& transaction);

    void restore(DeleteNodeTransaction& transaction);

    std::tuple<size_t, std::unique_ptr<AddEdgeTransaction>> addEdge(const size_t start, const size_t end);

    void resetEdge(const size_t id, const size_t start, const size_t end);

    void commit(AddEdgeTransaction& transaction);

    void restore(AddEdgeTransaction& transaction);

    // Return pointer to DeleteEdgeTransaction
    std::unique_ptr<DeleteEdgeTransaction> deleteEdge(const size_t id);

    void commit(DeleteEdgeTransaction& transaction);

    void restore(DeleteEdgeTransaction& transaction);

    size_t getNumberOfNodes() const;

    size_t getNumberOfEdges() const;

    size_t getNumberOfValidNodes() const;

    size_t getNumberOfValidEdges() const;

    void administrate();

    void print() const;

private:
    std::vector<Point> nodes_;
    std::vector<Edge> edges_;
    std::vector<std::vector<size_t>> nodeEdges_;
};

//--------------------------------

#endif // SIMPLE_MESH__HPP
