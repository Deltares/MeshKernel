#include <iostream>
#include <vector>
#include <utility>

#include "CompoundTransaction.hpp"
#include "SimpleMesh.hpp"
#include "Transaction.hpp"
#include "TransactionStack.hpp"

std::vector<Point> getPoints(const int count)
{
    int pointCount = count * count;

    std::vector<Point> points(pointCount);

    double originX = 0.0;
    double originY = 0.0;

    double deltaX = 1.0;
    double deltaY = 1.0;

    double y = originY;
    int pos = 0;

    for (int i = 1; i <= count; ++i)
    {
        double x = originX;

        for (int j = 1; j <= count; ++j)
        {
            points[pos] = Point(x, y);
            x += deltaX;
            ++pos;
        }

        y += deltaY;
    }

    return points;
}

std::vector<Edge> getEdges(const int count)
{
    std::vector<Edge> edges;

    for (int i = 0; i < count; ++i)
    {
        for (int j = 0; j < count; ++j)
        {
            int start;
            int end;

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

void addCentreNodes(SimpleMesh& mesh, TransactionStack& transactionStack, int count)
{

    std::unique_ptr<CompoundTransaction> transactions = std::make_unique<CompoundTransaction>();

    double deltaX = 1.0;
    double deltaY = 1.0;

    double originX = 0.0 + 0.5 * deltaX;
    double originY = 0.0 + 0.5 * deltaY;

    double y = originY;

    for (int i = 0; i < count - 1; ++i)
    {
        double x = originX;

        for (int j = 0; j < count - 1; ++j)
        {
            int node1 = i * count + j;
            int node2 = i * count + j + 1;
            int node3 = (i + 1) * count + j + 1;
            int node4 = (i + 1) * count + j;

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

void addNodeDeletion(SimpleMesh& mesh, TransactionStack& transactionStack, int count)
{
    int nodeToDelete = count * count;

    std::unique_ptr<CompoundTransaction> transactions = std::make_unique<CompoundTransaction>();

    for (int i = 0; i < count - 1; ++i)
    {
        for (int j = 0; j < count - 1; ++j)
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
    int count = 10; // 3164;
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
