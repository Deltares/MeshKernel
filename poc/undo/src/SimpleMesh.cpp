#include "SimpleMesh.hpp"

#include <memory>
#include <iostream>

SimpleMesh::SimpleMesh(const std::vector<Point>& nodes,
                       const std::vector<Edge>& edges) : nodes_(nodes), edges_(edges)
{
    administrate();
}

size_t SimpleMesh::getNumberOfNodes() const
{
    return nodes_.size();
}

size_t SimpleMesh::getNumberOfEdges() const
{
    return edges_.size();
}

const Point& SimpleMesh::getNode(const size_t id) const
{
    return nodes_[id];
}

const Edge& SimpleMesh::getEdge(const size_t id) const
{
    return edges_[id];
}

const std::vector<size_t>& SimpleMesh::getNodeConnectivity(const size_t id) const
{
    return nodeEdges_[id];
}

std::tuple<size_t, TransactionPtr> SimpleMesh::addNode(const Point& p)
{
    nodes_.resize(nodes_.size() + 1);

    TransactionPtr transaction = std::make_unique<AddNodeTransaction>(*this, nodes_.size() - 1, p);

    resetNode(nodes_.size() - 1, p);
    return {nodes_.size() - 1, std::move(transaction)};
}

void SimpleMesh::resetNode(const size_t id, const Point& p)
{
#ifdef ADD_LOGGING
    std::cout << "% reset node " << id << "  " << p.x() << "  " << p.y() << std::endl;
#endif

    nodes_[id] = p;
}

void SimpleMesh::commit(AddNodeTransaction& transaction)
{
    resetNode(transaction.nodeId(), transaction.point());
}

void SimpleMesh::restore(AddNodeTransaction& transaction)
{
    resetNode(transaction.nodeId(), Point());
}

void SimpleMesh::commit(AddEdgeTransaction& transaction)
{
    resetEdge(transaction.edgeId(), transaction.start(), transaction.end());
}

void SimpleMesh::restore(AddEdgeTransaction& transaction)
{
    resetEdge(transaction.edgeId(), core::nullValueId, core::nullValueId);
}

TransactionPtr SimpleMesh::deleteNode(const size_t id)
{
    std::unique_ptr<DeleteNodeTransaction> transaction = std::make_unique<DeleteNodeTransaction>(*this, id, nodes_[id]);

    for (size_t i = 0; i < nodeEdges_[id].size(); ++i)
    {
        transaction->emplace_back(deleteEdge(nodeEdges_[id][i]));
    }

#ifdef ADD_LOGGING
    std::cout << "% deleting node " << id << "  " << nodes_[id].x() << "  " << nodes_[id].y() << std::endl;
#endif
    nodes_[id].setInvalid();
    return transaction;
}

void SimpleMesh::commit(DeleteNodeTransaction& transaction)
{
    resetNode(transaction.nodeId(), Point());
}

void SimpleMesh::restore(DeleteNodeTransaction& transaction)
{
    resetNode(transaction.nodeId(), transaction.point());
}

std::tuple<size_t, TransactionPtr> SimpleMesh::addEdge(const size_t start, const size_t end)
{
    size_t index = edges_.size();
    edges_.resize(edges_.size() + 1);

    TransactionPtr transaction = std::make_unique<AddEdgeTransaction>(*this, index, start, end);

    resetEdge(index, start, end);

    return {index, std::move(transaction)};
}

void SimpleMesh::resetEdge(const size_t id, const size_t start, const size_t end)
{
#ifdef ADD_LOGGING
    std::cout << "% reset edge " << id << "  " << start << "  " << end << std::endl;
#endif

    edges_[id].start() = start;
    edges_[id].end() = end;
}

TransactionPtr SimpleMesh::deleteEdge(const size_t id)
{
    std::unique_ptr<DeleteEdgeTransaction> transaction = std::make_unique<DeleteEdgeTransaction>(*this, id, edges_[id].start(), edges_[id].end());

#ifdef ADD_LOGGING
    std::cout << "% deleting edge " << id << "  " << edges_[id].start() << "  " << edges_[id].end() << std::endl;
#endif
    resetEdge(id, core::nullValueId, core::nullValueId);
    return transaction;
}

void SimpleMesh::commit(DeleteEdgeTransaction& transaction)
{
    resetEdge(transaction.edgeId(), core::nullValueId, core::nullValueId);
}

void SimpleMesh::restore(DeleteEdgeTransaction& transaction)
{
    resetEdge(transaction.edgeId(), transaction.start(), transaction.end());
}

size_t SimpleMesh::getNumberOfValidNodes() const
{
    size_t count = 0;

    for (size_t i = 0; i < nodes_.size(); ++i)
    {
        if (nodes_[i].isValid())
        {
            ++count;
        }
    }

    return count;
}

size_t SimpleMesh::getNumberOfValidEdges() const
{
    size_t count = 0;

    for (size_t i = 0; i < edges_.size(); ++i)
    {
        if (edges_[i].isValid())
        {
            ++count;
        }
    }

    return count;
}

void SimpleMesh::print() const
{
    std::cout << "nullId = " << core::nullValueId << ";" << std::endl;

    std::cout << "nodex = zeros (" << nodes_.size() << ", 1);" << std::endl;
    std::cout << "nodey = zeros (" << nodes_.size() << ", 1);" << std::endl;

    for (size_t i = 0; i < nodes_.size(); ++i)
    {
        std::cout << "nodex (" << i + 1 << ") = " << nodes_[i].x() << ";" << std::endl;
        std::cout << "nodey (" << i + 1 << ") = " << nodes_[i].y() << ";" << std::endl;
    }

    std::cout << "edges = zeros (" << edges_.size() << ", 2);" << std::endl;

    for (size_t i = 0, count = 1; i < edges_.size(); ++i)
    {
        if (edges_[i].isValid())
        {
            std::cout << "edges (" << count << ", 1) = " << edges_[i].start() + 1 << ";" << std::endl;
            std::cout << "edges (" << count << ", 2) = " << edges_[i].end() + 1 << ";" << std::endl;
            ++count;
        }
        else
        {
            std::cout << "edges (" << count << ", 1) = " << core::nullValueId << ";" << std::endl;
            std::cout << "edges (" << count << ", 2) = " << core::nullValueId << ";" << std::endl;
            ++count;
        }
    }
}

void SimpleMesh::administrate()
{
    nodeEdges_.clear();
    nodeEdges_.resize(nodes_.size());

    for (size_t i = 0; i < edges_.size(); ++i)
    {
        size_t start = edges_[i].start();
        size_t end = edges_[i].end();

        if (!edges_[i].isValid())
        {
            continue;
        }

        nodeEdges_[start].push_back(i);
        nodeEdges_[end].push_back(i);
    }
}
