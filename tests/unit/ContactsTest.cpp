
#include <gtest/gtest.h>

#include "TestUtils/MakeMeshes.hpp"
#include <MeshKernel/Contacts.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh1D.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Polygons.hpp>

TEST(Contacts, ComputeSingleConnections1dOutside2dMesh)
{
    // Create 1d mesh
    std::vector<meshkernel::Point> nodes{
        {-16.1886410000000, 0.89018900000000},
        {-16.1464995876014, 9.78201442138723},
        {-16.1043581752028, 18.6738398427745},
        {-16.0622167628042, 27.5656652641617},
        {-15.7539488236928, 36.1966603330179},
        {-6.86476658679268, 36.4175095626911},
        {2.02441565010741, 36.6383587923643},
        {10.9135970000000, 36.8592080000000}};
    std::vector<meshkernel::Edge> edges{{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}};
    const auto mesh1d = std::make_shared<meshkernel::Mesh1D>(edges, nodes, meshkernel::Projection::cartesian);

    // Create 2d mesh
    const auto mesh2d = MakeRectangularMeshForTesting(4, 4, 10, meshkernel::Projection::cartesian, {0.0, 0.0});

    // Create contacts
    std::vector<bool> onedNodeMask(nodes.size(), true);
    meshkernel::Contacts contacts(mesh1d, mesh2d, onedNodeMask);

    // Set the polygon where to generate the contacts
    std::vector<meshkernel::Point> polygonPoints{{-30, -20}, {40, -20}, {40, 50}, {-40, 50}, {-30, -20}};
    meshkernel::Polygons polygon(polygonPoints, meshkernel::Projection::cartesian);

    // Execute
    contacts.ComputeSingleConnections(polygon);

    //Assert
    ASSERT_EQ(4, contacts.m_mesh1dIndices.size());
    ASSERT_EQ(4, contacts.m_mesh2dIndices.size());

    ASSERT_EQ(1, contacts.m_mesh1dIndices[0]);
    ASSERT_EQ(2, contacts.m_mesh1dIndices[1]);
    ASSERT_EQ(3, contacts.m_mesh1dIndices[2]);
    ASSERT_EQ(6, contacts.m_mesh1dIndices[3]);

    ASSERT_EQ(0, contacts.m_mesh2dIndices[0]);
    ASSERT_EQ(1, contacts.m_mesh2dIndices[1]);
    ASSERT_EQ(2, contacts.m_mesh2dIndices[2]);
    ASSERT_EQ(2, contacts.m_mesh2dIndices[3]);
}

TEST(Contacts, ComputeSingleConnections1dMeshInside2dMesh)
{

    // Create 1d mesh
    std::vector<meshkernel::Point> nodes{
        {1.73493900000000, -7.6626510000000},
        {2.35659313023165, 1.67281447902331},
        {5.38347452702839, 10.3513746546384},
        {14.2980910429074, 12.4797224193970},
        {22.9324017677239, 15.3007317677239},
        {25.3723169493137, 24.1623588554512},
        {25.8072280000000, 33.5111870000000}};
    std::vector<meshkernel::Edge> edges{{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}};
    const auto mesh1d = std::make_shared<meshkernel::Mesh1D>(edges, nodes, meshkernel::Projection::cartesian);

    // Create 2d mesh
    const auto mesh2d = MakeRectangularMeshForTesting(4, 4, 10, meshkernel::Projection::cartesian, {0.0, 0.0});

    // Create contacts
    std::vector<bool> onedNodeMask(nodes.size(), true);
    meshkernel::Contacts contacts(mesh1d, mesh2d, onedNodeMask);

    // Set the polygon where to generate the contacts
    std::vector<meshkernel::Point> polygonPoints{{-30, -20}, {40, -20}, {40, 50}, {-40, 50}, {-30, -20}};
    meshkernel::Polygons polygon(polygonPoints, meshkernel::Projection::cartesian);

    // Execute
    contacts.ComputeSingleConnections(polygon);

    //Assert
    ASSERT_EQ(5, contacts.m_mesh1dIndices.size());
    ASSERT_EQ(5, contacts.m_mesh2dIndices.size());

    ASSERT_EQ(1, contacts.m_mesh1dIndices[0]);
    ASSERT_EQ(2, contacts.m_mesh1dIndices[1]);
    ASSERT_EQ(3, contacts.m_mesh1dIndices[2]);
    ASSERT_EQ(4, contacts.m_mesh1dIndices[3]);
    ASSERT_EQ(5, contacts.m_mesh1dIndices[4]);

    ASSERT_EQ(0, contacts.m_mesh2dIndices[0]);
    ASSERT_EQ(1, contacts.m_mesh2dIndices[1]);
    ASSERT_EQ(4, contacts.m_mesh2dIndices[2]);
    ASSERT_EQ(7, contacts.m_mesh2dIndices[3]);
    ASSERT_EQ(8, contacts.m_mesh2dIndices[4]);
}

TEST(Contacts, ComputeMultipleConnections1dMeshInside2dMesh)
{

    // Create 1d mesh
    std::vector<meshkernel::Point> nodes{
        {1.73493900000000, -7.6626510000000},
        {2.35659313023165, 1.67281447902331},
        {5.38347452702839, 10.3513746546384},
        {14.2980910429074, 12.4797224193970},
        {22.9324017677239, 15.3007317677239},
        {25.3723169493137, 24.1623588554512},
        {25.8072280000000, 33.5111870000000}};
    std::vector<meshkernel::Edge> edges{{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}};
    const auto mesh1d = std::make_shared<meshkernel::Mesh1D>(edges, nodes, meshkernel::Projection::cartesian);

    // Create 2d mesh
    const auto mesh2d = MakeRectangularMeshForTesting(4, 4, 10, meshkernel::Projection::cartesian, {0.0, 0.0});

    // Create contacts
    std::vector<bool> onedNodeMask(nodes.size(), true);
    meshkernel::Contacts contacts(mesh1d, mesh2d, onedNodeMask);

    // Set the polygon where to generate the contacts
    std::vector<meshkernel::Point> polygonPoints{{-30, -20}, {40, -20}, {40, 50}, {-40, 50}, {-30, -20}};
    meshkernel::Polygons polygon(polygonPoints, meshkernel::Projection::cartesian);

    // Execute
    contacts.ComputeMultipleConnections();

    //Assert
    ASSERT_EQ(5, contacts.m_mesh1dIndices.size());
    ASSERT_EQ(5, contacts.m_mesh2dIndices.size());

    ASSERT_EQ(1, contacts.m_mesh1dIndices[0]);
    ASSERT_EQ(2, contacts.m_mesh1dIndices[1]);
    ASSERT_EQ(3, contacts.m_mesh1dIndices[2]);
    ASSERT_EQ(4, contacts.m_mesh1dIndices[3]);
    ASSERT_EQ(5, contacts.m_mesh1dIndices[4]);

    ASSERT_EQ(0, contacts.m_mesh2dIndices[0]);
    ASSERT_EQ(1, contacts.m_mesh2dIndices[1]);
    ASSERT_EQ(4, contacts.m_mesh2dIndices[2]);
    ASSERT_EQ(7, contacts.m_mesh2dIndices[3]);
    ASSERT_EQ(8, contacts.m_mesh2dIndices[4]);
}