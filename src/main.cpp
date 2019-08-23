#define _USE_MATH_DEFINES
#include <iostream>
#include <chrono>
#include "Mesh.cpp"
#include "CartesianPoint.cpp"

int main()
{
    // test 1
    //{
    //    //One gets the edges
    //    std::vector<Point> nodes;
    //    nodes.push_back(Point{ 0.0,0.0 });
    //    nodes.push_back(Point{ 0.0,10.0 });
    //    nodes.push_back(Point{ 10.0,0.0 });
    //    nodes.push_back(Point{ 10.0,10.0 });

    //    std::vector<Edge> edges;
    //    // Local edges
    //    edges.push_back({ 0, 2 });
    //    edges.push_back({ 1, 3 });
    //    edges.push_back({ 0, 1 });
    //    edges.push_back({ 2, 3 });

    //    // now build node-edge mapping
    //    std::vector<std::vector<size_t>> edgeNode(nodes.size(), std::vector<size_t>(8, 0));
    //    std::vector<size_t> numEdgesPeNode(nodes.size(), 0);
    //    BuildNodeEdge(edges, edgeNode, numEdgesPeNode);
    //    SortEdges(nodes, edges, numEdgesPeNode, edgeNode);

    //    // now find cells
    //    std::vector<std::vector<size_t>> faces;
    //    size_t numFaces = 0;
    //    std::vector<std::vector<size_t>> faceIndexsesForEachEdge(2, std::vector<size_t>(edges.size()));
    //    findFaces(nodes, edges, numEdgesPeNode, edgeNode, 4, faces, numFaces, faceIndexsesForEachEdge);


    //    // the number of found faces is
    //    std::cout << "Number of found cells " << faces.size();
    //}


    // test 2
    {

        int n = 11; //x
        int m = 11; //y

        std::cout << "start adding edges " << std::endl;
        auto start(std::chrono::steady_clock::now());

        // inject the point type (cartesian or spherical)
        using Mesh = Mesh<CartesianPoint>;

        std::vector<std::vector<int>> indexesValues(n, std::vector<int>(m));
        std::vector<CartesianPoint> nodes(n*m);
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
        Mesh mesh(edges,nodes, false);

        end = std::chrono::steady_clock::now();
        std::cout << "Elapsed time " << std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count() << " s " << std::endl;

        // the number of found faces is
        auto faces = mesh.getFaces();
        std::cout << "Number of found cells " << faces.size() << std::endl;
        std::cout << "First face " <<  faces[0][0] << " " << faces[0][1] << " " << faces[0][2] << " " << faces[0][3] << std::endl;
        std::cout << "Second face " << faces[1][0] << " " << faces[1][1] << " " << faces[1][2] << " " << faces[1][3] << std::endl;


        ////-------------------------------------------------------------------------------------------------------------------------------------//
        //// constructing surface mesh ds
        //std::cout << "Costructing cgal surface mesh data structure ... " << std::endl;
        //start = std::chrono::steady_clock::now();
        //surface_mesh mesh;
        //// add nodes
        //std::vector<surface_mesh_vertex_descriptor> surface_mesh_nodes_descriptors(nodes.size());
        //std::vector<surface_mesh_point> surface_mesh_nodes(nodes.size());
        //for (int n = 0; n < nodes.size(); n++)
        //{
        //    surface_mesh_nodes[n] = { nodes[n].first, nodes[n].second };
        //    surface_mesh_nodes_descriptors[n] = mesh.add_vertex(surface_mesh_nodes[n]);
        //}
        //// add faces
        //std::vector<surface_mesh_face_descriptor> surface_mesh_faces(faces.size());
        //std::vector<std::vector<surface_mesh_vertex_descriptor>> surface_mesh_face_node(faces.size());
        //for (int f = 0; f < faces.size(); f++)
        //{
        //    for (int n = 0; n < faces[f].size(); n++)
        //    {
        //        surface_mesh_face_node[f].push_back(surface_mesh_nodes_descriptors[faces[f][n]]);
        //    }
        //    surface_mesh_faces[f] = mesh.add_face(surface_mesh_face_node[f]);
        //}

        //end = std::chrono::steady_clock::now();
        //std::cout << "surface mesh : Number of faces " << mesh.num_faces() <<", number of edges "<< mesh.num_edges() << std::endl;
        //std::cout << "Elapsed time " << std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count() << " s " << std::endl;

        ////-------------------------------------------------------------------------------------------------------------------------------------//
        //// compute triangles circumcenters
        //auto v = surface_mesh_faces[0];
        //CGAL::circumcenter(surface_mesh_nodes[0], surface_mesh_nodes[1], surface_mesh_nodes[2]);

    }

    //getchar();
    return 0;
}