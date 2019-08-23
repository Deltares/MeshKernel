#include "GridGeomTest.hpp"
#include <chrono>

TEST(TestMesh, OneQuad) 
{
    using Mesh = Mesh<CartesianPoint>;

    //One gets the edges
    std::vector<CartesianPoint> nodes;
    nodes.push_back(CartesianPoint{ 0.0,0.0 });
    nodes.push_back(CartesianPoint{ 0.0,10.0 });
    nodes.push_back(CartesianPoint{ 10.0,0.0 });
    nodes.push_back(CartesianPoint{ 10.0,10.0 });

    std::vector<Mesh::Edge> edges;
    // Local edges
    edges.push_back({ 0, 2 });
    edges.push_back({ 1, 3 });
    edges.push_back({ 0, 1 });
    edges.push_back({ 2, 3 });

    // now build node-edge mapping
    std::vector<std::vector<size_t>> edgeNode(nodes.size(), std::vector<size_t>(8, 0));
    std::vector<size_t> numEdgesPeNode(nodes.size(), 0);
    // now build node-edge mapping
    Mesh mesh(edges, nodes, false);

	//check values
    EXPECT_EQ(1, mesh.getFaces().size());
    auto circumcenters = mesh.getFacesCircumcenters();
	EXPECT_DOUBLE_EQ(5.0, circumcenters[0].x);
    EXPECT_DOUBLE_EQ(5.0, circumcenters[0].y);
}


TEST(TestMesh, MillionQuads)
{
    int n = 1001; //x
    int m = 1001; //y

    std::cout << "start adding edges " << std::endl;
    auto start(std::chrono::steady_clock::now());

    // inject the point type (cartesian or spherical)
    using Mesh = Mesh<CartesianPoint>;

    std::vector<std::vector<int>> indexesValues(n, std::vector<int>(m));
    std::vector<CartesianPoint> nodes(n * m);
    size_t nodeIndex = 0;
    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            indexesValues[i][j] = i + j * n;
            nodes[nodeIndex] = { (double)i, (double)j };
            nodeIndex++;
        }
    }

    std::vector<Mesh::Edge> edges((n - 1) * m + (m - 1) * n);
    size_t edgeIndex = 0;
    for (int j = 0; j < m; ++j)
    {
        for (int i = 0; i < n - 1; ++i)
        {
            edges[edgeIndex] = { indexesValues[i][j], indexesValues[i + 1][j] };
            edgeIndex++;
        }
    }

    for (int j = 0; j < m - 1; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            edges[edgeIndex] = { indexesValues[i][j + 1], indexesValues[i][j] };
            edgeIndex++;
        }
    }
    auto end(std::chrono::steady_clock::now());
    std::cout << "Elapsed time " << std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count() << " s " << std::endl;

    std::cout << "start finding cells " << std::endl;
    start = std::chrono::steady_clock::now();
    // now build node-edge mapping
    Mesh mesh(edges, nodes, false);

    end = std::chrono::steady_clock::now();
    std::cout << "Elapsed time " << std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count() << " s " << std::endl;

    // the number of found faces is
    auto faces = mesh.getFaces();
    std::cout << "Number of found cells " << faces.size() << std::endl;
    std::cout << "First face " << faces[0][0] << " " << faces[0][1] << " " << faces[0][2] << " " << faces[0][3] << std::endl;
    std::cout << "Second face " << faces[1][0] << " " << faces[1][1] << " " << faces[1][2] << " " << faces[1][3] << std::endl;
}


