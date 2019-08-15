// Surface mesh part
#define _USE_MATH_DEFINES
#include <vector>
#include <cmath>
#include <chrono>
#include <numeric>
#include <algorithm>
#include <iostream>

typedef std::pair<size_t, size_t> Edge;
typedef std::pair<double, double> Point;

// Cgal surface mesh
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Simple_cartesian<double> surface_mesh_kernel;
typedef CGAL::Surface_mesh<surface_mesh_kernel::Point_2> surface_mesh;
typedef surface_mesh_kernel::Point_2 surface_mesh_point;
typedef surface_mesh::Vertex_index surface_mesh_vertex_descriptor;
typedef surface_mesh::Face_index surface_mesh_face_descriptor;


//Set node admin
void buildNodeEdgeMapping(std::vector<Edge>& edges, std::vector<std::vector<size_t>>& edgeNode, std::vector<size_t>& numEdgesPeNode)
{
    // assume no duplicated linkscma
    for (size_t e = 0; e < edges.size(); e++)
    {
        const size_t firstNode = edges[e].first;
        const size_t secondNode = edges[e].second;
        edgeNode[firstNode][numEdgesPeNode[firstNode]] = e;
        edgeNode[secondNode][numEdgesPeNode[secondNode]] = e;
        numEdgesPeNode[firstNode]++;
        numEdgesPeNode[secondNode]++;
    }

    // resize
    for (int node = 0; node < edgeNode.size(); node++)
    {
        edgeNode[node].resize(numEdgesPeNode[node]);
    }
}


void sortEdgesNodeCounterClockWise(const std::vector<Point>& nodes,
    const std::vector<Edge>& edges,
    const std::vector<size_t>& numEdgesPeNode,
    std::vector<std::vector<size_t>>& edgeNode)
{
    const double minimumDelta = 1e-14;
    for (size_t node = 0; node < nodes.size(); node++)
    {
        double phi0 = 0.0;
        double phi;
        std::vector<double> edgesAngles(numEdgesPeNode[node], 0.0);
        for (size_t edgeIndex = 0; edgeIndex < numEdgesPeNode[node]; edgeIndex++)
        {


            auto firstNode = edges[edgeNode[node][edgeIndex]].first;
            auto secondNode = edges[edgeNode[node][edgeIndex]].second;
            if (secondNode == node)
            {
                secondNode = firstNode;
                firstNode = node;
            }

            double deltaX = nodes[secondNode].first - nodes[firstNode].first;
            double deltaY = nodes[secondNode].second - nodes[firstNode].second;

            if (abs(deltaX) < minimumDelta && abs(deltaY) < minimumDelta)
            {
                if (deltaY < 0.0)
                {
                    phi = -M_PI / 2.0;
                }
                else
                {
                    phi = M_PI / 2.0;
                }
            }
            else
            {
                phi = atan2(deltaY, deltaX);
            }


            if (edgeIndex == 0)
            {
                phi0 = phi;
            }

            edgesAngles[edgeIndex] = phi - phi0;
            if (edgesAngles[edgeIndex] < 0.0)
            {
                edgesAngles[edgeIndex] = edgesAngles[edgeIndex] + 2.0 * M_PI;
            }
        }

        // Performing sorting
        std::vector<size_t> indexes(numEdgesPeNode[node]);
        std::vector<size_t> edgeNodeCopy{ edgeNode[node] };
        iota(indexes.begin(), indexes.end(), 0);
        sort(indexes.begin(), indexes.end(), [&edgesAngles](size_t i1, size_t i2) {return edgesAngles[i1] < edgesAngles[i2]; });

        for (size_t edgeIndex = 0; edgeIndex < numEdgesPeNode[node]; edgeIndex++)
        {
            edgeNode[node][edgeIndex] = edgeNodeCopy[indexes[edgeIndex]];
        }
    }
}

void findFaces(const std::vector<Point>& nodes,
    const std::vector<Edge>& edges,
    const std::vector<size_t>& numEdgesPeNode,
    const std::vector<std::vector<size_t>>& edgeNode,
    const int& numEdges,
    std::vector<std::vector<size_t>>& faces,
    size_t& numFaces,
    std::vector<std::vector<size_t>>& faceIndexsesForEachEdge)
{

    std::vector<size_t> numFacesForEachEdge(edges.size(), 0);
    std::vector<size_t> foundEdges(numEdges);
    std::vector<size_t> foundNodes(numEdges);

    for (size_t node = 0; node < nodes.size(); node++)
    {
        for (size_t firstEdgeLocalIndex = 0; firstEdgeLocalIndex < numEdgesPeNode[node]; firstEdgeLocalIndex++)
        {
            size_t indexFoundNodes = 0;
            size_t indexFoundEdges = 0;

            size_t currentEdge = edgeNode[node][firstEdgeLocalIndex];
            size_t currentNode = node;
            foundEdges[indexFoundEdges] = currentEdge;
            foundNodes[indexFoundNodes] = currentNode;
            int numFoundEdges = 1;

            if (numFacesForEachEdge[currentEdge] >= 2)
            {
                continue;
            }

            while (numFoundEdges < numEdges)
            {

                // the new node index
                if (numFacesForEachEdge[currentEdge] >= 2)
                {
                    break;
                }

                currentNode = edges[currentEdge].first + edges[currentEdge].second - currentNode;
                indexFoundNodes++;
                foundNodes[indexFoundNodes] = currentNode;

                int edgeIndex = 0;
                for (size_t localEdgeIndex = 0; localEdgeIndex < numEdgesPeNode[currentNode]; localEdgeIndex++)
                {
                    if (edgeNode[currentNode][localEdgeIndex] == currentEdge)
                    {
                        edgeIndex = localEdgeIndex;
                        break;
                    }
                }

                edgeIndex = edgeIndex - 1;
                if (edgeIndex < 0)
                {
                    edgeIndex = edgeIndex + numEdgesPeNode[currentNode];
                }
                if (edgeIndex > numEdgesPeNode[currentNode] - 1)
                {
                    edgeIndex = edgeIndex - numEdgesPeNode[currentNode];
                }
                currentEdge = edgeNode[currentNode][edgeIndex];
                indexFoundEdges++;
                foundEdges[indexFoundEdges] = currentEdge;

                numFoundEdges++;
            }

            // now check if the last node coincides
            if (numFacesForEachEdge[currentEdge] >= 2)
            {
                continue;
            }
            currentNode = edges[currentEdge].first + edges[currentEdge].second - currentNode;
            //indexFoundNodes++;
            //foundNodes[indexFoundNodes] = currentNode;


            if (currentNode == foundNodes[0])
            {
                // a cell has been found
                bool isFaceAlreadyFound = false;
                for (size_t localEdgeIndex = 0; localEdgeIndex <= indexFoundEdges; localEdgeIndex++)
                {
                    if (numFacesForEachEdge[foundEdges[localEdgeIndex]] >= 2)
                    {
                        isFaceAlreadyFound = true;
                        break;
                    }
                }
                if (isFaceAlreadyFound)
                {
                    continue;
                }

                bool allEdgesHaveAFace = true;
                for (size_t localEdgeIndex = 0; localEdgeIndex <= indexFoundEdges; localEdgeIndex++)
                {
                    if (numFacesForEachEdge[foundEdges[localEdgeIndex]] < 1)
                    {
                        allEdgesHaveAFace = false;
                        break;
                    }
                }

                bool isAnAlreadyFoundBoundaryFace = true;
                for (size_t localEdgeIndex = 0; localEdgeIndex <= indexFoundEdges - 1; localEdgeIndex++)
                {
                    if (faceIndexsesForEachEdge[0][foundEdges[localEdgeIndex]] != faceIndexsesForEachEdge[0][foundEdges[localEdgeIndex + 1]])
                    {
                        isAnAlreadyFoundBoundaryFace = false;
                        break;
                    }
                }

                if (allEdgesHaveAFace && isAnAlreadyFoundBoundaryFace)
                {
                    continue;
                }

                // increase numFacesForEachEdge 
                numFaces += 1;
                for (size_t localEdgeIndex = 0; localEdgeIndex <= indexFoundEdges; localEdgeIndex++)
                {
                    numFacesForEachEdge[foundEdges[localEdgeIndex]] += 1;
                    int numFace = numFacesForEachEdge[foundEdges[localEdgeIndex]];
                    faceIndexsesForEachEdge[numFace - 1][foundEdges[localEdgeIndex]] = numFaces;
                }

                // store the result
                faces.push_back(foundNodes);
            }
        }
    }
}

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
    //    buildNodeEdgeMapping(edges, edgeNode, numEdgesPeNode);
    //    sortEdgesNodeCounterClockWise(nodes, edges, numEdgesPeNode, edgeNode);

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

        int n = 1001; //x
        int m = 1001; //y

        std::cout << "start adding edges " << std::endl;
        auto start(std::chrono::steady_clock::now());

        std::vector<std::vector<int>> indexesValues(n, std::vector<int>(m));
        std::vector<Point> nodes(n*m);
        size_t nodeIndex = 0;
        for (int j = 0; j < m; ++j)
        {
            for (int i = 0; i < n; ++i)
            {
                indexesValues[i][j] = i + j * n;
                nodes[nodeIndex] = { i, j };
                nodeIndex++;
            }
        }


        std::vector<Edge> edges((n - 1) * m + (m - 1) * n);
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
        std::vector<std::vector<size_t>> edgeNode(nodes.size(), std::vector<size_t>(8, 0));
        std::vector<size_t> numEdgesPeNode(nodes.size(), 0);
        buildNodeEdgeMapping(edges, edgeNode, numEdgesPeNode);
        sortEdgesNodeCounterClockWise(nodes, edges, numEdgesPeNode, edgeNode);

        // now find cells
        std::vector<std::vector<size_t>> faces;
        size_t numFaces = 0;
        std::vector<std::vector<size_t>> faceIndexsesForEachEdge(2, std::vector<size_t>(edges.size()));
        findFaces(nodes, edges, numEdgesPeNode, edgeNode, 3, faces, numFaces, faceIndexsesForEachEdge);
        findFaces(nodes, edges, numEdgesPeNode, edgeNode, 4, faces, numFaces, faceIndexsesForEachEdge);
        findFaces(nodes, edges, numEdgesPeNode, edgeNode, 5, faces, numFaces, faceIndexsesForEachEdge);
        findFaces(nodes, edges, numEdgesPeNode, edgeNode, 6, faces, numFaces, faceIndexsesForEachEdge);

        end = std::chrono::steady_clock::now();
        std::cout << "Elapsed time " << std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count() << " s " << std::endl;

        // the number of found faces is
        std::cout << "Number of found cells " << faces.size() << std::endl;
        std::cout << "First face " << faces[0][0] << " " << faces[0][1] << " " << faces[0][2] << " " << faces[0][3] << std::endl;
        std::cout << "Second face " << faces[1][0] << " " << faces[1][1] << " " << faces[1][2] << " " << faces[1][3] << std::endl;


        //-------------------------------------------------------------------------------------------------------------------------------------//
        // constructing surface mesh ds
        std::cout << "Costructing cgal surface mesh data structure ... " << std::endl;
        start = std::chrono::steady_clock::now();
        surface_mesh mesh;
        // add nodes
        std::vector<surface_mesh_vertex_descriptor> surface_mesh_nodes_descriptors(nodes.size());
        std::vector<surface_mesh_point> surface_mesh_nodes(nodes.size());
        for (int n = 0; n < nodes.size(); n++)
        {
            surface_mesh_nodes[n] = { nodes[n].first, nodes[n].second };
            surface_mesh_nodes_descriptors[n] = mesh.add_vertex(surface_mesh_nodes[n]);
        }
        // add faces
        std::vector<surface_mesh_face_descriptor> surface_mesh_faces(faces.size());
        std::vector<std::vector<surface_mesh_vertex_descriptor>> surface_mesh_face_node(faces.size());
        for (int f = 0; f < faces.size(); f++)
        {
            for (int n = 0; n < faces[f].size(); n++)
            {
                surface_mesh_face_node[f].push_back(surface_mesh_nodes_descriptors[faces[f][n]]);
            }
            surface_mesh_faces[f] = mesh.add_face(surface_mesh_face_node[f]);
        }

        end = std::chrono::steady_clock::now();
        std::cout << "surface mesh : Number of faces " << mesh.num_faces() <<", number of edges "<< mesh.num_edges() << std::endl;
        std::cout << "Elapsed time " << std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count() << " s " << std::endl;

        //-------------------------------------------------------------------------------------------------------------------------------------//
        // compute triangles circumcenters
        auto v = surface_mesh_faces[0];
        CGAL::circumcenter(surface_mesh_nodes[0], surface_mesh_nodes[1], surface_mesh_nodes[2]);

    }

    //getchar();
    return 0;
}
