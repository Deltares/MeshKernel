// TODO: duplicated face detection issue for large number of edges
#ifndef MESH_CPP
#define MESH_CPP

#define _USE_MATH_DEFINES
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>

#include "Mesh.hpp"
#include "Constants.cpp"
#include "Operations.cpp"

template <typename Point>
Mesh<Point>::Mesh(const std::vector<GridGeom::Edge>& edges, const std::vector<Point>& nodes) :
    m_edges(edges), m_nodes(nodes)
{
    //allocate
    m_nodesNumEdges.resize(m_nodes.size());
    m_nodesEdges.resize(m_nodes.size(), std::vector<size_t>(GridGeom::maximumNumberOfEdgesPerNode, 0));
    m_edgesNumFaces.resize(m_edges.size());
    m_numFaces = 0;
    m_edgesFaces.resize(edges.size(), std::vector<int>(2, -1));


    //run administration and find the faces
    NodeAdministration();
    SortEdgesInCounterClockWiseOrder();

    for(int numEdgesPerFace =3; numEdgesPerFace<= 4; numEdgesPerFace++)
    {
        findFaces(numEdgesPerFace);
    }
 
    // find mesh circumcenters
    faceCircumcenters(1.0);

    //compute faces areas and centers of mass
    facesAreasAndCentersOfMass();
};

template <typename Point>
void Mesh<Point>::NodeAdministration()
{
    // assume no duplicated links
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
};

template <typename Point>
void Mesh<Point>::SortEdgesInCounterClockWiseOrder()
{
    std::vector<double> edgesAngles(GridGeom::maximumNumberOfEdgesPerNode, 0.0);
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

            double deltaX = GridGeom::Operations<Point>::getDx(m_nodes[secondNode], m_nodes[firstNode]);
            double deltaY = GridGeom::Operations<Point>::getDy(m_nodes[secondNode], m_nodes[firstNode]);
            if (abs(deltaX) < GridGeom::minimumDeltaCoordinate && abs(deltaY) < GridGeom::minimumDeltaCoordinate)
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
template <typename Point>
void Mesh<Point>::findFaces(const int& numEdges)
{

    std::vector<size_t> foundEdges(numEdges);
    std::vector<size_t> foundNodes(numEdges); 

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

template <typename Point>
void Mesh<Point>::faceCircumcenters(const double& weightCircumCenter)
{
    m_facesCircumcenters.resize(m_facesNodes.size());
    std::vector<Point> middlePoints(GridGeom::maximumNumberOfNodesPerFace);
    std::vector<Point> normals(GridGeom::maximumNumberOfNodesPerFace);
    std::vector<Point> localFace(GridGeom::maximumNumberOfNodesPerFace + 1); //closed face
    Point centerOfMass;
    const int maximumNumberCircumcenterIterations = 100;
    for (int f = 0; f < m_facesNodes.size(); f++)
    {
        // for triangles, for now assume cartesian kernel
        size_t numberOfFaceNodes = m_facesNodes[f].size();
        if (numberOfFaceNodes == 3)
        {
            GridGeom::Operations<Point>::circumcenterOfTriangle(m_nodes[m_facesNodes[f][0]], m_nodes[m_facesNodes[f][1]], m_nodes[m_facesNodes[f][2]], m_facesCircumcenters[f]);
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
                    GridGeom::Operations<Point>::normalVector(localFace[n], localFace[nextNode], middlePoints[n], normals[n]);
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
                            double dx = GridGeom::Operations<Point>::getDx(middlePoints[n], estimatedCircumCenter);
                            double dy = GridGeom::Operations<Point>::getDy(middlePoints[n], estimatedCircumCenter);
                            double increment = -0.1 * GridGeom::dotProduct(dx, dy, normals[n].x, normals[n].y);
                            GridGeom::Operations<Point>::add(estimatedCircumCenter, normals[n], increment);
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

                if (!isCircumcenterInside)
                {
                    for (int n = 0; n < numberOfFaceNodes; n++)
                    {
                        int nextNode = n + 1;
                        if (nextNode == numberOfFaceNodes) nextNode = 0;
                        Point intersection;
                        bool isLineCrossing = lineCrossing(centerOfMass, m_facesCircumcenters[f], localFace[n], localFace[nextNode], intersection);
                        if (isLineCrossing)
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

template <typename Point>
void Mesh<Point>::facesAreasAndCentersOfMass()
{
    // polygon coordinates 
    m_faceArea.resize(m_facesNodes.size());
    m_facesMasscenters.resize(m_facesNodes.size());
    std::vector<Point> localFace(GridGeom::maximumNumberOfNodesPerFace + 1);
    Point centerOfMass;
    double area = 0.0;
    for (int f = 0; f < m_facesNodes.size(); f++)
    {
        size_t numberOfFaceNodes = m_facesNodes[f].size();

        for (int n = 0; n < numberOfFaceNodes; n++)
        {
            localFace[n] = m_nodes[m_facesNodes[f][n]];
        }
        localFace[numberOfFaceNodes] = localFace[0];
        area = 0.0;
        faceAreaAndCenterOfMass(localFace, numberOfFaceNodes, area, centerOfMass);
        m_faceArea[f] = area;
        m_facesMasscenters[f].x = centerOfMass.x;
        m_facesMasscenters[f].y = centerOfMass.y;
    }
}


#endif
