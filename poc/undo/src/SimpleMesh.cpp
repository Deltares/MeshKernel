#include "SimpleMesh.hpp"

#include <memory>
#include <iostream>

SimpleMesh::SimpleMesh(const std::vector<Point>& nodes,
                       const std::vector<Edge>& edges) : nodes_(nodes), edges_(edges)
{
    administrate();
}

int SimpleMesh::getNumberOfNodes() const
{
    return nodes_.size();
}

int SimpleMesh::getNumberOfEdges() const
{
    return edges_.size();
}

const Point& SimpleMesh::getNode(const int id) const
{
    return nodes_[id];
}

const Edge& SimpleMesh::getEdge(const int id) const
{
    return edges_[id];
}

const std::vector<int>& SimpleMesh::getNodeConnectivity(const int id) const
{
    return nodeEdges_[id];
}

std::tuple<int, TransactionPtr> SimpleMesh::addNode(const Point& p)
{
    nodes_.resize(nodes_.size() + 1);

#ifdef NULL_TRANSACTION
    TransactionPtr transaction;
#else
    TransactionPtr transaction = std::make_unique<AddNodeTransaction>(*this, nodes_.size() - 1, p);
#endif

    resetNode(nodes_.size() - 1, p);
    return {nodes_.size() - 1, std::move(transaction)};
}

void SimpleMesh::resetNode(const int id, const Point& p)
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

TransactionPtr SimpleMesh::deleteNode(const int id)
{
#ifdef NULL_TRANSACTION
    TransactionPtr transaction;
#else
    std::unique_ptr<DeleteNodeTransaction> transaction = std::make_unique<DeleteNodeTransaction>(*this, id, nodes_[id]);

    for (size_t i = 0; i < nodeEdges_[id].size(); ++i)
    {
        transaction->emplace_back(deleteEdge(nodeEdges_[id][i]));
    }
#endif

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

std::tuple<int, TransactionPtr> SimpleMesh::addEdge(const int start, const int end)
{
    int index = edges_.size();
    edges_.resize(edges_.size() + 1);

#ifdef NULL_TRANSACTION
    TransactionPtr transaction;
#else
    TransactionPtr transaction = std::make_unique<AddEdgeTransaction>(*this, index, start, end);
#endif

    resetEdge(index, start, end);

    return {index, std::move(transaction)};
}

void SimpleMesh::resetEdge(const int id, const int start, const int end)
{
#ifdef ADD_LOGGING
    std::cout << "% reset edge " << id << "  " << start << "  " << end << std::endl;
#endif

    edges_[id].start() = start;
    edges_[id].end() = end;
}

TransactionPtr SimpleMesh::deleteEdge(const int id)
{
#ifdef NULL_TRANSACTION
    TransactionPtr transaction;
#else
    std::unique_ptr<DeleteEdgeTransaction> transaction = std::make_unique<DeleteEdgeTransaction>(*this, id, edges_[id].start(), edges_[id].end());
#endif

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

int SimpleMesh::getNumberOfValidNodes() const
{
    int count = 0;

    for (size_t i = 0; i < nodes_.size(); ++i)
    {
        if (nodes_[i].isValid())
        {
            ++count;
        }
    }

    return count;
}

int SimpleMesh::getNumberOfValidEdges() const
{
    int count = 0;

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
        int start = edges_[i].start();
        int end = edges_[i].end();

        if (!edges_[i].isValid())
        {
            continue;
        }

        nodeEdges_[start].push_back(i);
        nodeEdges_[end].push_back(i);
    }
}
