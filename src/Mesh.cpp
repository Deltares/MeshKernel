// Surface mesh part
// TODO: Use CGAL kernel points to easy switch between cartesian and spherical basic computations (e.g. triangle circumcenters or distances)
// Spherical coordinate kernel https://doc.cgal.org/latest/Circular_kernel_3/index.html (ja spheric)

#define _USE_MATH_DEFINES
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <CGAL/Kernel/global_functions.h>

#include"AlgebraicOperations.cpp"

template<typename Point>
class Mesh
{
public:

    typedef std::pair<size_t, size_t>        Edge;

    Mesh(const std::vector<Edge>& edges, const std::vector<Point>& nodes, const bool isSpheric) :
        _edges(edges), _nodes(nodes), _isSpheric(isSpheric)
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
        
        // find mesh circumcenters
        //faceCircumcenters(1.0);
    };

    const std::vector<std::vector<size_t>>& getFaces() const
    {
        return _facesNodes;
    }

    const std::vector<Point>& getFacesCircumcenters() const
    {
        return _facesCircumcenters;
    }
    
private:

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

                double deltaX = Point::getDx(_nodes[secondNode], _nodes[firstNode]);
                double deltaY = Point::getDy(_nodes[secondNode], _nodes[firstNode]);
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
                    _facesNodes.push_back(foundNodes);
                    _facesEdges.push_back(foundEdges);
                }
            }
        }
    }
    
    void facePolygon(const int faceIndex, std::vector<Point>& localPolygon, std::vector<size_t>& localEdges, std::vector<size_t>& localEdgeNumber)
    {
        localPolygon.resize(_facesNodes[faceIndex].size());
        localEdgeNumber.resize(_facesNodes[faceIndex].size());
        localEdges.resize(_facesNodes[faceIndex].size());
        for(int n=0;n<_facesNodes[faceIndex].size();n++)
        {
            localPolygon[n] = _nodes[_facesNodes[faceIndex][n]];
            localEdges[n] = _numFacesForEachEdge[_facesEdges[faceIndex][n]];
            localEdgeNumber[n] = n;
        }
    }

    void faceCircumcenters(const double & weightCircumCenter)
    {
        _facesCircumcenters.resize(_facesNodes.size());
        std::vector<Point> localPolygon;
        std::vector<size_t> localEdges;
        std::vector<size_t> localEdgeNumber;
        for (int f = 0; f < _facesNodes.size(); f++)
        {
            // for triangles, for now assume cartesian kernel
            if(_facesNodes[f].size()==3)
            {
                //_nodes[_facesNodes[f][0]] 
                //_nodes[_facesNodes[f][1]]
                //_nodes[_facesNodes[f][2]]
                //auto circumCenter = CGAL::circumcenter(_nodes[_facesNodes[f][0]], _nodes[_facesNodes[f][1]], _nodes[_facesNodes[f][2]]);
                //_facesCircumcenters[f] = { circumCenter.x(),circumCenter.y() };
                _facesCircumcenters[f] = { 0.0,0.0 };
            }
            else
            {
                double xCenter = 0.0;
                double yCenter = 0.0;

                facePolygon(f, localPolygon, localEdges, localEdgeNumber);

                int numberOfInteriorEdges = 0;
                for (int n = 0; n < localPolygon.size(); n++)
                {
                    xCenter += localPolygon[n].x;
                    yCenter += localPolygon[n].y;
                    if (localEdges[n]==2)
                    {
                        numberOfInteriorEdges += 1;
                    }
                }

                // average centers
                xCenter = xCenter / localPolygon.size();
                yCenter = yCenter / localPolygon.size();
                
                if( numberOfInteriorEdges == 0 )
                {
                    _facesCircumcenters[f].x=  xCenter;
                    _facesCircumcenters[f].y = yCenter;

                }
                else
                {
                    Point estimatedCircumCenter{ xCenter, yCenter };
                    const double eps = 1e-3;
                    
                    std::vector<Point> middlePoints(localPolygon.size());
                    std::vector<Point> normals(localPolygon.size());
                    for (int n = 0; n < localPolygon.size(); n++)
                    {
                        int nextNode = n + 1;
                        if (nextNode == localPolygon.size()) nextNode = 0;
                        middlePoints[n] = { 0.5 * (localPolygon[n].x + localPolygon[nextNode].x),  0.5 * (localPolygon[n].y + localPolygon[nextNode].y) };
                        normals[n] = Point::normalVector(localPolygon[n], localPolygon[nextNode], middlePoints[n]);
                    }

                    for (int iter = 0; iter < 100; iter++)
                    {
                        Point previousCircumCenter = estimatedCircumCenter;
                        // compute quantities

                        for (int n = 0; n < localPolygon.size(); n++)
                        {
                            if (localEdges[n] == 2 || localPolygon.size() == 3)
                            {
                                int nextNode =  n + 1; 
                                if (nextNode == localPolygon.size()) nextNode = 0;
                                double dx = Point::getDx(middlePoints[n], estimatedCircumCenter);
                                double dy = Point::getDy(middlePoints[n], estimatedCircumCenter);
                                double increment = - 0.1 * AlgebraicOperations::dotProduct(dx, dy, normals[n].x, normals[n].y);
                                Point::add(estimatedCircumCenter, normals[n], increment);
                            }
                        }
                        if (iter > 0 && 
                            abs(estimatedCircumCenter.x - previousCircumCenter.x) < eps && 
                            abs(estimatedCircumCenter.y - previousCircumCenter.y) < eps)
                        {
                            _facesCircumcenters[f] = estimatedCircumCenter;
                            break;
                        }
                    }
                }
            }
        }

        if (weightCircumCenter <= 1.0 && weightCircumCenter >=0.0)
        {
            //double localWeightCircumCenter = 1.0;
            //if( localPolygon.size()>3 )
            //{
            //    localWeightCircumCenter = weightCircumCenter;
            //}
            // to finish
        }

    }

    std::vector<Edge> _edges;                                  //KN
    std::vector<Point> _nodes;                                 //KN
    std::vector<std::vector<size_t>> _nodeEdges;               //NOD
    std::vector<size_t> _numEdgesPeNode;                       //NMK
    std::vector<std::vector<size_t>> _facesNodes;              //netcell%Nod
    std::vector<std::vector<size_t>> _facesEdges;              //netcell%lin
    std::vector<Point>   _facesCircumcenters;                  //xw
    std::vector<size_t> _numFacesForEachEdge;                  //LNN
    size_t _faceIndex;                                         //NUMP
    std::vector<std::vector<size_t>> _faceIndexsesForEachEdge; //LNE
    bool _isSpheric;                                           //jsferic      

    //constants
    const double minimumDeltaCoordinate = 1e-14;
    const size_t maximumNumberOfEdgesPerNode = 8;
    const double dcenterinside = 1.0;                          //Force cell center inside cell with factor dcenterinside, 1: confined by cell edges, 0 : at mass center

};
