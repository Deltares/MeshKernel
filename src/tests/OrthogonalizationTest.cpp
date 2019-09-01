#include "GridGeomTest.hpp"

TEST(OrthogonalizationTest, TestAspectRatio)
{
    using Mesh = Mesh<cartesianPoint>;

    //One gets the edges
    std::vector<cartesianPoint> nodes;
    nodes.push_back(cartesianPoint{ 0.0,0.0 });
    nodes.push_back(cartesianPoint{ 0.0,10.0 });
    nodes.push_back(cartesianPoint{ 10.0,0.0 });
    nodes.push_back(cartesianPoint{ 10.0,10.0 });

    std::vector<Edge> edges;
    // Local edges
    edges.push_back({ 0, 2 });
    edges.push_back({ 1, 3 });
    edges.push_back({ 0, 1 });
    edges.push_back({ 2, 3 });

    // now build node-edge mapping
    Mesh mesh(edges, nodes);

    // compute the aspect ratio
    mesh.aspectRatio();

    EXPECT_EQ(1, mesh.m_aspectRatio[0]);
    EXPECT_EQ(1, mesh.m_aspectRatio[1]);
    EXPECT_EQ(1, mesh.m_aspectRatio[2]);
    EXPECT_EQ(1, mesh.m_aspectRatio[3]);
}


