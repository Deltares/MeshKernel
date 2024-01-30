#include <iostream>
#include <vector>
#include <utility>

#include "CompoundTransaction.hpp"
#include "SimpleMesh.hpp"
#include "Transaction.hpp"
#include "TransactionStack.hpp"

std::vector<Point> getPoints(const size_t count)
{
    size_t pointCount = count * count;

    std::vector<Point> points(pointCount);

    double originX = 0.0;
    double originY = 0.0;

    double deltaX = 1.0;
    double deltaY = 1.0;

    double y = originY;
    size_t pos = 0;

    for (size_t i = 1; i <= count; ++i)
    {
        double x = originX;

        for (size_t j = 1; j <= count; ++j)
        {
            points[pos] = Point(x, y);
            x += deltaX;
            ++pos;
        }

        y += deltaY;
    }

    return points;
}

std::vector<Edge> getEdges(const size_t count)
{
    std::vector<Edge> edges;

    for (size_t i = 0; i < count; ++i)
    {
        for (size_t j = 0; j < count; ++j)
        {
            size_t start;
            size_t end;

            // Horizontal
            if (j < count - 1)
            {
                start = i * count + j;
                end = i * count + j + 1;
                edges.push_back(Edge(start, end));
            }

            // Vertical
            if (i < count - 1)
            {
                start = i * count + j;
                end = (i + 1) * count + j;
                edges.push_back(Edge(start, end));
            }
        }
    }

    return edges;
}

void addCentreNodes(SimpleMesh& mesh, TransactionStack& transactionStack, size_t count)
{

    std::unique_ptr<CompoundTransaction> transactions = std::make_unique<CompoundTransaction>();

    double deltaX = 1.0;
    double deltaY = 1.0;

    double originX = 0.0 + 0.5 * deltaX;
    double originY = 0.0 + 0.5 * deltaY;

    double y = originY;

    for (size_t i = 0; i < count - 1; ++i)
    {
        double x = originX;

        for (size_t j = 0; j < count - 1; ++j)
        {
            size_t node1 = i * count + j;
            size_t node2 = i * count + j + 1;
            size_t node3 = (i + 1) * count + j + 1;
            size_t node4 = (i + 1) * count + j;

            auto [newNodeId, transaction] = mesh.addNode(Point(x, y));
#ifndef NULL_TRANSACTION
            transactions->emplace_back(std::move(transaction));
#endif
            auto [edgeId1, transaction1] = mesh.addEdge(node1, newNodeId);

#ifndef NULL_TRANSACTION
            transactions->emplace_back(std::move(transaction1));
#endif

            auto [edgeId2, transaction2] = mesh.addEdge(node2, newNodeId);

#ifndef NULL_TRANSACTION
            transactions->emplace_back(std::move(transaction2));
#endif

            auto [edgeId3, transaction3] = mesh.addEdge(node3, newNodeId);

#ifndef NULL_TRANSACTION
            transactions->emplace_back(std::move(transaction3));
#endif

            auto [edgeId4, transaction4] = mesh.addEdge(node4, newNodeId);
#ifndef NULL_TRANSACTION

            transactions->emplace_back(std::move(transaction4));
#endif

            x += deltaX;
        }
        y += deltaY;
    }

    mesh.administrate();
    transactionStack.emplace_back(std::move(transactions));
}

void addNodeDeletion(SimpleMesh& mesh, TransactionStack& transactionStack, size_t count)
{
    size_t nodeToDelete = count * count;

    std::unique_ptr<CompoundTransaction> transactions = std::make_unique<CompoundTransaction>();

    for (size_t i = 0; i < count - 1; ++i)
    {
        for (size_t j = 0; j < count - 1; ++j)
        {
            TransactionPtr transaction = mesh.deleteNode(nodeToDelete);
#ifndef NULL_TRANSACTION
            transactions->emplace_back(std::move(transaction));
#endif
            ++nodeToDelete;
        }
    }

    mesh.administrate();
    transactionStack.emplace_back(std::move(transactions));
}

int main()
{
    size_t count = 10; // 3164;
    std::vector<Point> points(getPoints(count));
    std::vector<Edge> edges(getEdges(count));

    SimpleMesh mesh(points, edges);

    TransactionStack transactions;

    addCentreNodes(mesh, transactions, count);
    addNodeDeletion(mesh, transactions, count);

    // Undo node deletion
    transactions.undo();

    // // Undo centre node and edge addition
    // transactions.undo ();

    mesh.administrate();
    mesh.print();
    return 0;
}
