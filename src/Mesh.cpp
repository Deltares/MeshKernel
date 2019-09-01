// Surface mesh part
// TODO: compute cell area
// Spherical coordinate kernel https://doc.cgal.org/latest/Circular_kernel_3/index.html (ja spheric)

#ifndef MESH_CPP
#define MESH_CPP

#define _USE_MATH_DEFINES
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include "Operations.cpp"

using namespace GridGeom;

template<typename Point>
class Mesh
{
public:

    Mesh(const std::vector<Edge>& edges, const std::vector<Point>& nodes) :
        m_edges(edges), m_nodes(nodes)
    {
        //allocate
        m_nodesNumEdges.resize(m_nodes.size());
        m_nodesEdges.resize(m_nodes.size(), std::vector<size_t>(maximumNumberOfEdgesPerNode, 0));
        m_edgesNumFaces.resize(m_edges.size());
        m_numFaces = 0;
        m_edgesFaces.resize(edges.size(), std::vector<size_t>(2));


        //run administration and find the faces
        NodeAdministration();
        SortEdgesInCounterClockWiseOrder();
        findFaces(3);
        findFaces(4);
        findFaces(5);
        findFaces(6);
        
        // find mesh circumcenters
        faceCircumcenters(1.0);

        //compute faces areas and centers of mass
        facesAreasAndCentersOfMass();
    };

    const std::vector<std::vector<size_t>>& getNodesEdges() const
    {
        return m_nodesEdges;
    }

    const  std::vector<size_t>& getNodesNumEdges() const
    {
        return m_nodesNumEdges;
    }

    const std::vector<std::vector<size_t>>& getFacesNodes() const
    {
        return m_facesNodes;
    }

    const std::vector<std::vector<size_t>>& getFacesEdges() const
    {
        return m_facesEdges;
    }

    const std::vector<Point>& getFacesCircumcenters() const
    {
        return m_facesCircumcenters;
    }

    const std::vector<size_t>& getEdgesNumFaces() const
    {
        return m_edgesNumFaces;
    }

    const std::vector<std::vector<size_t>>& getEdgesFaces() const
    {
        return m_edgesFaces;
    }

    std::vector<double> m_aspectRatio;

    bool aspectRatio()
    {
        std::vector<std::vector<double>> averageEdgesLength(m_edges.size(), std::vector<double>(2, doubleMissingValue));
        std::vector<double> averageFlowEdgesLength(m_edges.size(), doubleMissingValue);
        std::vector<bool> curvilinearGridIndicator(m_nodes.size(), true);
        std::vector<double> edgesLength(m_edges.size(), 0.0);
        m_aspectRatio.resize(m_edges.size(), 0.0);

        for (size_t e = 0; e < m_edges.size(); e++)
        {
            size_t first = m_edges[e].first;
            size_t second = m_edges[e].second;

            if (first == second) continue;
            double edgeLength = Operations<Point>::distance(m_nodes[first], m_nodes[second]);
            edgesLength[e] = edgeLength;

            Point leftCenter;
            Point rightCenter;
            if (m_edgesNumFaces[e] > 0)
            {
                leftCenter = m_facesCircumcenters[m_edgesFaces[e][0]];
            }
            else
            {
                leftCenter = m_nodes[first];
            }

            //find right cell center, if it exists
            if (m_edgesNumFaces[e] == 2)
            {
                rightCenter = m_facesCircumcenters[m_edgesFaces[e][1]];
            }
            else
            {
                //otherwise, make ghost node by imposing boundary condition
                double dinry = Operations<Point>::innerProductTwoSegments(m_nodes[first], m_nodes[second], m_nodes[first], leftCenter);
                dinry = dinry / std::max(edgeLength * edgeLength, minimumEdgeLength);

                double x0_bc = (1.0 - dinry) * m_nodes[first].x + dinry * m_nodes[second].x;
                double y0_bc = (1.0 - dinry) * m_nodes[first].y + dinry * m_nodes[second].y;
                rightCenter.x = 2.0 * x0_bc - leftCenter.x;
                rightCenter.y = 2.0 * y0_bc - leftCenter.y;
            }

            averageFlowEdgesLength[e] = Operations<Point>::distance(leftCenter, rightCenter);
        }

        // Compute normal length
        for (int f = 0; f < m_facesNodes.size(); f++)
        {
            size_t numberOfFaceNodes = m_facesNodes[f].size();
            if (numberOfFaceNodes < 3) continue;

            for (int n = 0; n < numberOfFaceNodes; n++)
            {
                if (numberOfFaceNodes != 4) curvilinearGridIndicator[m_facesNodes[f][n]] = false;
                size_t edgeIndex = m_facesEdges[f][n];

                if (m_edgesNumFaces[edgeIndex] < 1) continue;

                //get the other links in the right numbering
                //TODO: ask why only 3 are requested, why an average lenght stored in averageEdgesLength is needed?
                //int kkm1 = n - 1; if (kkm1 < 0) kkm1 = kkm1 + numberOfFaceNodes;
                //int kkp1 = n + 1; if (kkp1 >= numberOfFaceNodes) kkp1 = kkp1 - numberOfFaceNodes;
                //
                //size_t klinkm1 = m_facesEdges[f][kkm1];
                //size_t klinkp1 = m_facesEdges[f][kkp1];
                //

                double edgeLength = edgesLength[edgeIndex];
                if (edgeLength != 0.0)
                {
                    m_aspectRatio[edgeIndex] = averageFlowEdgesLength[edgeIndex] / edgeLength;
                }

                //quads
                if (numberOfFaceNodes == 4)
                {
                    int kkp2 = n + 2; if (kkp2 >= numberOfFaceNodes) kkp2 = kkp2 - numberOfFaceNodes;
                    size_t klinkp2 = m_facesEdges[f][kkp2];
                    edgeLength = 0.5 * (edgesLength[edgeIndex] + edgesLength[klinkp2]);
                }

                if (averageEdgesLength[edgeIndex][0] == doubleMissingValue)
                {
                    averageEdgesLength[edgeIndex][0] = edgeLength;
                }
                else
                {
                    averageEdgesLength[edgeIndex][1] = edgeLength;
                }
            }
        }

        if (curvilinearToOrthogonalRatio == 1.0)
            return true;

        for (size_t e = 0; e < m_edges.size(); e++)
        {
            size_t first = m_edges[e].first;
            size_t second = m_edges[e].second;

            if (first < 0 || second < 0) continue;
            if (m_edgesNumFaces[e] < 1) continue;
            // Consider only quads
            if (!curvilinearGridIndicator[first] || !curvilinearGridIndicator[second]) continue;

            if (m_edgesNumFaces[e] == 1)
            {
                if (averageEdgesLength[e][0] != 0.0 && averageEdgesLength[e][0] != doubleMissingValue)
                {
                    m_aspectRatio[e] = averageFlowEdgesLength[e] / averageEdgesLength[e][0];
                }
            }
            else
            {
                if (averageEdgesLength[e][0] != 0.0 && averageEdgesLength[e][1] != 0.0 &&
                    averageEdgesLength[e][0] != doubleMissingValue && averageEdgesLength[e][1] != doubleMissingValue)
                {
                    m_aspectRatio[e] = curvilinearToOrthogonalRatio * m_aspectRatio[e] +
                        (1.0 - curvilinearToOrthogonalRatio) * averageFlowEdgesLength[e] / (0.5 * (averageEdgesLength[e][0] + averageEdgesLength[e][1]));
                }
            }
        }
        return true;
    }

private:

    // Set node admin
    void NodeAdministration()
    {
        // assume no duplicated linkscma
        for (size_t e = 0; e < m_edges.size(); e++)
        {
            const size_t firstNode = m_edges[e].first;
            const size_t secondNode = m_edges[e].second;
            m_nodesEdges[firstNode][m_nodesNumEdges[firstNode]] = e;
            m_nodesEdges[secondNode][m_nodesNumEdges[secondNode]] = e;
            m_nodesNumEdges[firstNode]++;
            m_nodesNumEdges[secondNode]++;
        }

        // resize
        for (int node = 0; node < m_nodesEdges.size(); node++)
        {
            m_nodesEdges[node].resize(m_nodesNumEdges[node]);
        }
    }

    void SortEdgesInCounterClockWiseOrder()
    {
        std::vector<double> edgesAngles(maximumNumberOfEdgesPerNode, 0.0);
        for (size_t node = 0; node < m_nodes.size(); node++)
        {
            double phi0 = 0.0;
            double phi;
            std::fill(edgesAngles.begin(), edgesAngles.end(), 0.0);
            for (size_t edgeIndex = 0; edgeIndex < m_nodesNumEdges[node]; edgeIndex++)
            {

                auto firstNode = m_edges[m_nodesEdges[node][edgeIndex]].first;
                auto secondNode = m_edges[m_nodesEdges[node][edgeIndex]].second;
                if (secondNode == node)
                {
                    secondNode = firstNode;
                    firstNode = node;
                }

                double deltaX = Operations<Point>::getDx(m_nodes[secondNode], m_nodes[firstNode]);
                double deltaY = Operations<Point>::getDy(m_nodes[secondNode], m_nodes[firstNode]);
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
            std::vector<size_t> indexes(m_nodesNumEdges[node]);
            std::vector<size_t> edgeNodeCopy{ m_nodesEdges[node] };
            iota(indexes.begin(), indexes.end(), 0);
            sort(indexes.begin(), indexes.end(), [&edgesAngles](size_t i1, size_t i2) {return edgesAngles[i1] < edgesAngles[i2]; });

            for (size_t edgeIndex = 0; edgeIndex < m_nodesNumEdges[node]; edgeIndex++)
            {
                m_nodesEdges[node][edgeIndex] = edgeNodeCopy[indexes[edgeIndex]];
            }
        }
    }

    // find cells
    void findFaces(const int& numEdges)
    {

        std::vector<size_t> foundEdges(numEdges);
        std::vector<size_t> foundNodes(numEdges); // close the

        for (size_t node = 0; node < m_nodes.size(); node++)
        {
            for (size_t firstEdgeLocalIndex = 0; firstEdgeLocalIndex < m_nodesNumEdges[node]; firstEdgeLocalIndex++)
            {
                size_t indexFoundNodes = 0;
                size_t indexFoundEdges = 0;

                size_t currentEdge = m_nodesEdges[node][firstEdgeLocalIndex];
                size_t currentNode = node;
                foundEdges[indexFoundEdges] = currentEdge;
                foundNodes[indexFoundNodes] = currentNode;
                int numFoundEdges = 1;

                if (m_edgesNumFaces[currentEdge] >= 2)
                {
                    continue;
                }

                while (numFoundEdges < numEdges)
                {

                    // the new node index
                    if (m_edgesNumFaces[currentEdge] >= 2)
                    {
                        break;
                    }

                    currentNode = m_edges[currentEdge].first + m_edges[currentEdge].second - currentNode;
                    indexFoundNodes++;
                    foundNodes[indexFoundNodes] = currentNode;

                    int edgeIndex = 0;
                    for (size_t localEdgeIndex = 0; localEdgeIndex < m_nodesNumEdges[currentNode]; localEdgeIndex++)
                    {
                        if (m_nodesEdges[currentNode][localEdgeIndex] == currentEdge)
                        {
                            edgeIndex = localEdgeIndex;
                            break;
                        }
                    }

                    edgeIndex = edgeIndex - 1;
                    if (edgeIndex < 0)
                    {
                        edgeIndex = edgeIndex + m_nodesNumEdges[currentNode];
                    }
                    if (edgeIndex > m_nodesNumEdges[currentNode] - 1)
                    {
                        edgeIndex = edgeIndex - m_nodesNumEdges[currentNode];
                    }
                    currentEdge = m_nodesEdges[currentNode][edgeIndex];
                    indexFoundEdges++;
                    foundEdges[indexFoundEdges] = currentEdge;

                    numFoundEdges++;
                }

                // now check if the last node coincides
                if (m_edgesNumFaces[currentEdge] >= 2)
                {
                    continue;
                }
                currentNode = m_edges[currentEdge].first + m_edges[currentEdge].second - currentNode;
                //indexFoundNodes++;
                //foundNodes[indexFoundNodes] = currentNode;

                if (currentNode == foundNodes[0])
                {
                    // a cell has been found
                    bool isFaceAlreadyFound = false;
                    for (size_t localEdgeIndex = 0; localEdgeIndex <= indexFoundEdges; localEdgeIndex++)
                    {
                        if (m_edgesNumFaces[foundEdges[localEdgeIndex]] >= 2)
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
                        if (m_edgesNumFaces[foundEdges[localEdgeIndex]] < 1)
                        {
                            allEdgesHaveAFace = false;
                            break;
                        }
                    }

                    bool isAnAlreadyFoundBoundaryFace = true;
                    for (size_t localEdgeIndex = 0; localEdgeIndex <= indexFoundEdges - 1; localEdgeIndex++)
                    {
                        if (m_edgesFaces[foundEdges[localEdgeIndex]][0] != m_edgesFaces[foundEdges[localEdgeIndex + 1]][0])
                        {
                            isAnAlreadyFoundBoundaryFace = false;
                            break;
                        }
                    }

                    if (allEdgesHaveAFace && isAnAlreadyFoundBoundaryFace)
                    {
                        continue;
                    }

                    // increase m_edgesNumFaces 
                    m_numFaces += 1;
                    for (size_t localEdgeIndex = 0; localEdgeIndex <= indexFoundEdges; localEdgeIndex++)
                    {
                        m_edgesNumFaces[foundEdges[localEdgeIndex]] += 1;
                        const int numFace = m_edgesNumFaces[foundEdges[localEdgeIndex]];
                        m_edgesFaces[foundEdges[localEdgeIndex]][numFace - 1] = m_numFaces - 1;
                    }

                    // store the result
                    m_facesNodes.push_back(foundNodes);
                    m_facesEdges.push_back(foundEdges);
                }
            }
        }
    }
    
    void faceCircumcenters(const double& weightCircumCenter)
    {
        m_facesCircumcenters.resize(m_facesNodes.size());
        std::vector<Point> middlePoints(maximumNumberOfNodesPerFace);
        std::vector<Point> normals(maximumNumberOfNodesPerFace);
        std::vector<Point> localFace(maximumNumberOfNodesPerFace + 1); //closed face
        Point centerOfMass;
        const int maximumNumberCircumcenterIterations = 100;
        for (int f = 0; f < m_facesNodes.size(); f++)
        {
            // for triangles, for now assume cartesian kernel
            size_t numberOfFaceNodes = m_facesNodes[f].size();
            if (numberOfFaceNodes == 3)
            {
                //m_nodes[m_facesNodes[f][0]] 
                //m_nodes[m_facesNodes[f][1]]
                //m_nodes[m_facesNodes[f][2]]
                //auto circumCenter = CGAL::circumcenter(m_nodes[m_facesNodes[f][0]], m_nodes[m_facesNodes[f][1]], m_nodes[m_facesNodes[f][2]]);
                //m_facesCircumcenters[f] = { circumCenter.x(),circumCenter.y() };
                m_facesCircumcenters[f] = { 0.0,0.0 };
            }
            else
            {
                double xCenter = 0.0;
                double yCenter = 0.0;
                size_t numberOfInteriorEdges = 0;
                for (int n = 0; n < numberOfFaceNodes; n++)
                {
                    localFace[n] = m_nodes[m_facesNodes[f][n]];
                    xCenter += m_nodes[m_facesNodes[f][n]].x;
                    yCenter += m_nodes[m_facesNodes[f][n]].y;
                    if (m_edgesNumFaces[m_facesEdges[f][n]] == 2)
                    {
                        numberOfInteriorEdges += 1;
                    }
                }
                localFace[numberOfFaceNodes] = localFace[0]; 
                // average centers
                centerOfMass = { xCenter / numberOfFaceNodes, yCenter / numberOfFaceNodes };

                if (numberOfInteriorEdges == 0)
                {
                    m_facesCircumcenters[f] = centerOfMass;
                }
                else
                {
                    Point estimatedCircumCenter = centerOfMass;
                    const double eps = 1e-3;

                    for (int n = 0; n < numberOfFaceNodes; n++)
                    {
                        int nextNode = n + 1;
                        if (nextNode == numberOfFaceNodes) nextNode = 0;
                        middlePoints[n].x = 0.5 * (localFace[n].x + localFace[nextNode].x);
                        middlePoints[n].y = 0.5 * (localFace[n].y + localFace[nextNode].y);
                        Operations<Point>::normalVector(localFace[n], localFace[nextNode], middlePoints[n], normals[n]);
                    }

                    Point previousCircumCenter = estimatedCircumCenter;
                    for (int iter = 0; iter < maximumNumberCircumcenterIterations; iter++)
                    {
                        previousCircumCenter = estimatedCircumCenter;
                        for (int n = 0; n < numberOfFaceNodes; n++)
                        {
                            if (m_edgesNumFaces[m_facesEdges[f][n]] == 2 || numberOfFaceNodes == 3)
                            {
                                int nextNode = n + 1;
                                if (nextNode == numberOfFaceNodes) nextNode = 0;
                                double dx = Operations<Point>::getDx(middlePoints[n], estimatedCircumCenter);
                                double dy = Operations<Point>::getDy(middlePoints[n], estimatedCircumCenter);
                                double increment = -0.1 * dotProduct(dx, dy, normals[n].x, normals[n].y);
                                Operations<Point>::add(estimatedCircumCenter, normals[n], increment);
                            }
                        }
                        if (iter > 0 &&
                            abs(estimatedCircumCenter.x - previousCircumCenter.x) < eps &&
                            abs(estimatedCircumCenter.y - previousCircumCenter.y) < eps)
                        {
                            m_facesCircumcenters[f] = estimatedCircumCenter;
                            break;
                        }
                    }
                }
                if (weightCircumCenter <= 1.0 && weightCircumCenter >= 0.0)
                {
                    double localWeightCircumCenter = 1.0;
                    if (numberOfFaceNodes > 3)
                    {
                        localWeightCircumCenter = weightCircumCenter;
                    }

                    for (int n = 0; n < numberOfFaceNodes; n++)
                    {
                        localFace[n].x = localWeightCircumCenter * localFace[n].x + (1.0 - localWeightCircumCenter) * centerOfMass.x;
                        localFace[n].y = localWeightCircumCenter * localFace[n].y + (1.0 - localWeightCircumCenter) * centerOfMass.y;
                    }

                    bool isCircumcenterInside = pointInPolygon(m_facesCircumcenters[f], localFace, numberOfFaceNodes + 1);

                    if(!isCircumcenterInside)
                    {
                        for (int n = 0; n < numberOfFaceNodes; n++)
                        {
                            int nextNode = n + 1;
                            if (nextNode == numberOfFaceNodes) nextNode = 0;
                            Point intersection;
                            bool isLineCrossing = lineCrossing(centerOfMass, m_facesCircumcenters[f], localFace[n], localFace[nextNode], intersection);
                            if(isLineCrossing)
                            {
                                m_facesCircumcenters[f] = intersection;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    void facesAreasAndCentersOfMass()
    {
        // polygon coordinates 
        m_faceArea.resize(m_facesNodes.size());
        m_facesMasscenters.resize(m_facesNodes.size());
        std::vector<Point> localFace(maximumNumberOfNodesPerFace + 1);
        for (int f = 0; f < m_facesNodes.size(); f++)
        {
            size_t numberOfFaceNodes = m_facesNodes[f].size();
            
            for (int n = 0; n < numberOfFaceNodes; n++)
            {
                localFace[n] = m_nodes[m_facesNodes[f][n]];
            }
            localFace[numberOfFaceNodes] = localFace[0];
            double area = 0.0;
            Point centerOfMass;
            faceAreaAndCenterOfMass(localFace, numberOfFaceNodes, area, centerOfMass);
            m_faceArea[f] = area;
            m_facesMasscenters[f] = centerOfMass;
        }
    }


    bool orthogonalize()
    {
        size_t maxNumNeighbours = std::max_element(m_nodesNumEdges);
        maxNumNeighbours += 1;

        std::vector< std::vector<size_t>> nodesNodes(m_nodes.size(), std::vector<size_t>(maxNumNeighbours, 0));
        //for each node, determine the neighbours
        for (size_t n = 0; n < m_nodes.size(); n++)
        {
            for (size_t nn = 0; nn < m_nodesNumEdges[n]; nn++)
            {
                Edge edge = m_edges[m_nodesEdges[n][nn]];
                size_t neighbour = edge.first == n ? edge.second : edge.first;
                nodesNodes[n] = neighbour;
            }
        }

        //Orthogonalization loop
        //size_t maxOuterIter = 2;
        //size_t maxInnerIter = 25;
        //std::vector<Point> nodesBackUp(m_nodes.size());
        //std::vector < std::vector < double>> weights(m_nodes.size(), std::vector < double>(maxNumNeighbours, 0.0));
        //std::vector < std::vector < double>> rightHandSide(m_nodes.size(), std::vector < double>(maxNumNeighbours, 0.0));
        //for (size_t outerIter = 0; outerIter < maxOuterIter; outerIter++)
        //{
        //    nodesBackUp = m_nodes;

        //    for (size_t innerIter = 0; innerIter < maxInnerIter; innerIter++)
        //    {

        //    }
        //}
    }



    const double m_dcenterinside = 1.0;

    std::vector<Edge> m_edges;                                  //KN
    std::vector<Point> m_nodes;                                 //KN
    std::vector<std::vector<size_t>> m_nodesEdges;              //NOD
    std::vector<size_t> m_nodesNumEdges;                        //NMK

    //edges
    std::vector<size_t> m_edgesNumFaces;                        //LNN
    std::vector<std::vector<size_t>> m_edgesFaces;              //LNE

    // faces
    std::vector<std::vector<size_t>> m_facesNodes;              //netcell%Nod, the nodes composing the faces, in ccw order
    std::vector<std::vector<size_t>> m_facesEdges;              //netcell%lin
    std::vector<Point>   m_facesCircumcenters;                  //xw
    std::vector<Point>   m_facesMasscenters;                    //the faces canters of mass


    size_t m_numFaces;                                          //NUMP

    std::vector<double> m_faceArea;                             //Face area



};

#endif
