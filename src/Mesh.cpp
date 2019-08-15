// Surface mesh part
#define _USE_MATH_DEFINES
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>

// cGAL surface mesh
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

class Mesh
{
public:

    typedef CGAL::Simple_cartesian<double> surface_mesh_kernel;
    typedef CGAL::Surface_mesh<surface_mesh_kernel::Point_2> surface_mesh;
    typedef surface_mesh_kernel::Point_2 surface_mesh_point;
    typedef surface_mesh::Vertex_index surface_mesh_vertex_descriptor;
    typedef surface_mesh::Face_index surface_mesh_face_descriptor;

    typedef std::pair<size_t, size_t> Edge;
    typedef std::pair<double, double> Point;

    Mesh(const std::vector<Edge>& edges, const std::vector<Point>& nodes) :
        _edges(edges), _nodes(nodes)
    {
        //allocate
        _numEdgesPeNode.resize(_nodes.size());
        _nodeEdges.resize(_nodes.size(), std::vector<size_t>(maximumNumberOfEdgesPerNode, 0));
        _numFacesForEachEdge.resize(_edges.size());
        _faceIndex = 0;
        _faceIndexsesForEachEdge.resize(edges.size(), std::vector<size_t>(2));


        //run administration and find the faces
        NodeAdministration();
        SortEdgesInCounterClockWiseOrder();
        findFaces(3);
        findFaces(4);
        findFaces(5);
        findFaces(6);
    };

    // Set node admin
    void NodeAdministration()
    {
        // assume no duplicated linkscma
        for (size_t e = 0; e < _edges.size(); e++)
        {
            const size_t firstNode = _edges[e].first;
            const size_t secondNode = _edges[e].second;
            _nodeEdges[firstNode][_numEdgesPeNode[firstNode]] = e;
            _nodeEdges[secondNode][_numEdgesPeNode[secondNode]] = e;
            _numEdgesPeNode[firstNode]++;
            _numEdgesPeNode[secondNode]++;
        }

        // resize
        for (int node = 0; node < _nodeEdges.size(); node++)
        {
            _nodeEdges[node].resize(_numEdgesPeNode[node]);
        }
    }

    void SortEdgesInCounterClockWiseOrder()
    {
        
        for (size_t node = 0; node < _nodes.size(); node++)
        {
            double phi0 = 0.0;
            double phi;
            std::vector<double> edgesAngles(_numEdgesPeNode[node], 0.0);
            for (size_t edgeIndex = 0; edgeIndex < _numEdgesPeNode[node]; edgeIndex++)
            {

                auto firstNode = _edges[_nodeEdges[node][edgeIndex]].first;
                auto secondNode = _edges[_nodeEdges[node][edgeIndex]].second;
                if (secondNode == node)
                {
                    secondNode = firstNode;
                    firstNode = node;
                }

                double deltaX = _nodes[secondNode].first - _nodes[firstNode].first;
                double deltaY = _nodes[secondNode].second - _nodes[firstNode].second;

                if (abs(deltaX) < minimumDeltaCoordinate && abs(deltaY) < minimumDeltaCoordinate)
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
            std::vector<size_t> indexes(_numEdgesPeNode[node]);
            std::vector<size_t> edgeNodeCopy{ _nodeEdges[node] };
            iota(indexes.begin(), indexes.end(), 0);
            sort(indexes.begin(), indexes.end(), [&edgesAngles](size_t i1, size_t i2) {return edgesAngles[i1] < edgesAngles[i2]; });

            for (size_t edgeIndex = 0; edgeIndex < _numEdgesPeNode[node]; edgeIndex++)
            {
                _nodeEdges[node][edgeIndex] = edgeNodeCopy[indexes[edgeIndex]];
            }
        }
    }

    // find cells
    void findFaces(const int& numEdges)
    {
        
        std::vector<size_t> foundEdges(numEdges);
        std::vector<size_t> foundNodes(numEdges);

        for (size_t node = 0; node < _nodes.size(); node++)
        {
            for (size_t firstEdgeLocalIndex = 0; firstEdgeLocalIndex < _numEdgesPeNode[node]; firstEdgeLocalIndex++)
            {
                size_t indexFoundNodes = 0;
                size_t indexFoundEdges = 0;

                size_t currentEdge = _nodeEdges[node][firstEdgeLocalIndex];
                size_t currentNode = node;
                foundEdges[indexFoundEdges] = currentEdge;
                foundNodes[indexFoundNodes] = currentNode;
                int numFoundEdges = 1;

                if (_numFacesForEachEdge[currentEdge] >= 2)
                {
                    continue;
                }

                while (numFoundEdges < numEdges)
                {

                    // the new node index
                    if (_numFacesForEachEdge[currentEdge] >= 2)
                    {
                        break;
                    }

                    currentNode = _edges[currentEdge].first + _edges[currentEdge].second - currentNode;
                    indexFoundNodes++;
                    foundNodes[indexFoundNodes] = currentNode;

                    int edgeIndex = 0;
                    for (size_t localEdgeIndex = 0; localEdgeIndex < _numEdgesPeNode[currentNode]; localEdgeIndex++)
                    {
                        if (_nodeEdges[currentNode][localEdgeIndex] == currentEdge)
                        {
                            edgeIndex = localEdgeIndex;
                            break;
                        }
                    }

                    edgeIndex = edgeIndex - 1;
                    if (edgeIndex < 0)
                    {
                        edgeIndex = edgeIndex + _numEdgesPeNode[currentNode];
                    }
                    if (edgeIndex > _numEdgesPeNode[currentNode] - 1)
                    {
                        edgeIndex = edgeIndex - _numEdgesPeNode[currentNode];
                    }
                    currentEdge = _nodeEdges[currentNode][edgeIndex];
                    indexFoundEdges++;
                    foundEdges[indexFoundEdges] = currentEdge;

                    numFoundEdges++;
                }

                // now check if the last node coincides
                if (_numFacesForEachEdge[currentEdge] >= 2)
                {
                    continue;
                }
                currentNode = _edges[currentEdge].first + _edges[currentEdge].second - currentNode;
                //indexFoundNodes++;
                //foundNodes[indexFoundNodes] = currentNode;


                if (currentNode == foundNodes[0])
                {
                    // a cell has been found
                    bool isFaceAlreadyFound = false;
                    for (size_t localEdgeIndex = 0; localEdgeIndex <= indexFoundEdges; localEdgeIndex++)
                    {
                        if (_numFacesForEachEdge[foundEdges[localEdgeIndex]] >= 2)
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
                        if (_numFacesForEachEdge[foundEdges[localEdgeIndex]] < 1)
                        {
                            allEdgesHaveAFace = false;
                            break;
                        }
                    }

                    bool isAnAlreadyFoundBoundaryFace = true;
                    for (size_t localEdgeIndex = 0; localEdgeIndex <= indexFoundEdges - 1; localEdgeIndex++)
                    {
                        if (_faceIndexsesForEachEdge[foundEdges[localEdgeIndex]][0] != _faceIndexsesForEachEdge[foundEdges[localEdgeIndex + 1]][0])
                        {
                            isAnAlreadyFoundBoundaryFace = false;
                            break;
                        }
                    }

                    if (allEdgesHaveAFace && isAnAlreadyFoundBoundaryFace)
                    {
                        continue;
                    }

                    // increase _numFacesForEachEdge 
                    _faceIndex += 1;
                    for (size_t localEdgeIndex = 0; localEdgeIndex <= indexFoundEdges; localEdgeIndex++)
                    {
                        _numFacesForEachEdge[foundEdges[localEdgeIndex]] += 1;
                        const int numFace = _numFacesForEachEdge[foundEdges[localEdgeIndex]];
                        _faceIndexsesForEachEdge[foundEdges[localEdgeIndex]][numFace - 1] = _faceIndex;
                    }

                    // store the result
                    _faces.push_back(foundNodes);
                }
            }
        }
    }

    
    const std::vector<std::vector<size_t>>& getFaces() const
    {
        return _faces;
    }

private:
    std::vector<Edge> _edges;                    //KN
    std::vector<Point> _nodes;                   //KN
    std::vector<std::vector<size_t>> _nodeEdges; //NOD
    std::vector<size_t> _numEdgesPeNode;         //NMK
    std::vector<std::vector<size_t>> _faces;     //netcell
    std::vector<size_t> _numFacesForEachEdge;    //LNN
    size_t _faceIndex;                           //NUMP
    std::vector<std::vector<size_t>> _faceIndexsesForEachEdge; //LNE

    //constants
    const double minimumDeltaCoordinate = 1e-14;
    const size_t maximumNumberOfEdgesPerNode = 8;
};
