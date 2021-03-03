#include <gtest/gtest.h>

#include "MeshKernel/CurvilinearGrid.hpp"
#include "MeshKernel/CurvilinearGridCreateUniform.hpp"
#include "MeshKernelApi/MakeMeshParameters.hpp"

#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Polygons.hpp>

TEST(CurvilinearGrid, MakeCurvilinearInInPolygon)
{
    //1 Setup
    std::vector<meshkernel::Point> polygonNodes{{302.002502, 472.130371},
                                                {144.501526, 253.128174},
                                                {368.752930, 112.876755},
                                                {707.755005, 358.879242},
                                                {301.252502, 471.380371},
                                                {302.002502, 472.130371}};

    const auto polygons = std::make_shared<meshkernel::Polygons>(polygonNodes, meshkernel::Projection::cartesian);

    meshkernelapi::MakeMeshParameters makeMeshParameters;
    makeMeshParameters.GridType = 0;
    makeMeshParameters.GridAngle = 0.0;
    makeMeshParameters.OriginXCoordinate = 0.0;
    makeMeshParameters.OriginYCoordinate = 0.0;
    makeMeshParameters.OriginZCoordinate = 0.0;
    makeMeshParameters.NumberOfColumns = 3;
    makeMeshParameters.NumberOfRows = 3;
    makeMeshParameters.XGridBlockSize = 100.0;
    makeMeshParameters.YGridBlockSize = 100.0;

    // 2 Execution
    meshkernel::CurvilinearGridCreateUniform curvilinearGridCreateUniform(makeMeshParameters, polygons);
    const auto [nodes, edges, gridIndices] = curvilinearGridCreateUniform.Compute().ConvertCurvilinearToNodesAndEdges();
    ASSERT_EQ(43, edges.size());
    ASSERT_EQ(27, nodes.size());
}

TEST(CurvilinearGrid, MakeCurvilinearInPolygonSpherical)
{
    //1 Setup
    std::vector<meshkernel::Point> polygonNodes{{302.002502, 472.130371},
                                                {144.501526, 253.128174},
                                                {368.752930, 112.876755},
                                                {707.755005, 358.879242},
                                                {301.252502, 471.380371},
                                                {302.002502, 472.130371}};

    const auto polygons = std::make_shared<meshkernel::Polygons>(polygonNodes, meshkernel::Projection::spherical);

    meshkernelapi::MakeMeshParameters makeMeshParameters;
    makeMeshParameters.GridType = 0;
    makeMeshParameters.GridAngle = 0.0;
    makeMeshParameters.OriginXCoordinate = 0.0;
    makeMeshParameters.OriginYCoordinate = 0.0;
    makeMeshParameters.OriginZCoordinate = 0.0;
    makeMeshParameters.NumberOfColumns = 3;
    makeMeshParameters.NumberOfRows = 3;
    makeMeshParameters.XGridBlockSize = 5000000.0; //resolution in meters (when using spherical coordinates distances are usually much larger)
    makeMeshParameters.YGridBlockSize = 5000000.0;

    // 2 Execution: function not producing grid points (points gets transformed in meters, therfore everything is outside)
    meshkernel::CurvilinearGridCreateUniform curvilinearGridCreateUniform(makeMeshParameters, polygons);
    const auto [nodes, edges, gridIndices] = curvilinearGridCreateUniform.Compute().ConvertCurvilinearToNodesAndEdges();
    ASSERT_EQ(0, nodes.size());
    ASSERT_EQ(0, edges.size());
}

TEST(CurvilinearGrid, MakeCurvilinearInEmptyPolygonSpherical)
{
    //1 Setup
    std::vector<meshkernel::Point> polygonNodes;
    const auto polygons = std::make_shared<meshkernel::Polygons>(polygonNodes, meshkernel::Projection::spherical);

    meshkernelapi::MakeMeshParameters makeMeshParameters;
    makeMeshParameters.GridType = 0;
    makeMeshParameters.GridAngle = 0.0;
    makeMeshParameters.OriginXCoordinate = 0.0;
    makeMeshParameters.OriginYCoordinate = 0.0;
    makeMeshParameters.OriginZCoordinate = 0.0;
    makeMeshParameters.NumberOfColumns = 3;
    makeMeshParameters.NumberOfRows = 3;
    makeMeshParameters.XGridBlockSize = 5000000.0; //resolution in meters (when using spherical coordinates distances are usually much larger)
    makeMeshParameters.YGridBlockSize = 5000000.0;

    // 2 Execution
    meshkernel::CurvilinearGridCreateUniform curvilinearGridCreateUniform(makeMeshParameters, polygons);
    const auto [nodes, edges, gridIndices] = curvilinearGridCreateUniform.Compute().ConvertCurvilinearToNodesAndEdges();

    meshkernel::Mesh2D mesh(edges, nodes, meshkernel::Projection::spherical);

    // 3 Assert
    ASSERT_EQ(24, mesh.GetNumEdges());
    ASSERT_EQ(16, mesh.GetNumNodes());

    // x coordinate
    ASSERT_EQ(0.0, mesh.m_nodes[0].x);
    ASSERT_EQ(makeMeshParameters.XGridBlockSize, mesh.m_nodes[1].x);
    ASSERT_EQ(makeMeshParameters.XGridBlockSize * 2, mesh.m_nodes[2].x);
    ASSERT_EQ(makeMeshParameters.XGridBlockSize * 3, mesh.m_nodes[3].x);

    ASSERT_EQ(0.0, mesh.m_nodes[4].x);
    ASSERT_EQ(makeMeshParameters.XGridBlockSize, mesh.m_nodes[5].x);
    ASSERT_EQ(makeMeshParameters.XGridBlockSize * 2, mesh.m_nodes[6].x);
    ASSERT_EQ(makeMeshParameters.XGridBlockSize * 3, mesh.m_nodes[7].x);

    ASSERT_EQ(0.0, mesh.m_nodes[8].x);
    ASSERT_EQ(makeMeshParameters.XGridBlockSize, mesh.m_nodes[9].x);
    ASSERT_EQ(makeMeshParameters.XGridBlockSize * 2, mesh.m_nodes[10].x);
    ASSERT_EQ(makeMeshParameters.XGridBlockSize * 3, mesh.m_nodes[11].x);

    // y coordinate
    ASSERT_EQ(0.0, mesh.m_nodes[0].y);
    ASSERT_EQ(0.0, mesh.m_nodes[1].y);
    ASSERT_EQ(0.0, mesh.m_nodes[2].y);
    ASSERT_EQ(0.0, mesh.m_nodes[3].y);

    ASSERT_EQ(makeMeshParameters.XGridBlockSize, mesh.m_nodes[4].y);
    ASSERT_EQ(makeMeshParameters.XGridBlockSize, mesh.m_nodes[5].y);
    ASSERT_EQ(makeMeshParameters.XGridBlockSize, mesh.m_nodes[6].y);
    ASSERT_EQ(makeMeshParameters.XGridBlockSize, mesh.m_nodes[7].y);

    ASSERT_EQ(3830222.2156113400, mesh.m_nodes[8].y);
    ASSERT_EQ(3830222.2156113400, mesh.m_nodes[9].y);
    ASSERT_EQ(3830222.2156113400, mesh.m_nodes[10].y);
    ASSERT_EQ(3830222.2156113400, mesh.m_nodes[11].y);
}