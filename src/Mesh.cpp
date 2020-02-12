// TODO: duplicated face detection issue for large number of edges
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>

#include "Mesh.hpp"
#include "Constants.cpp"
#include "Operations.cpp"

bool GridGeom::Mesh::Set(const std::vector<Edge>& edges, const std::vector<Point>& nodes, Projections projection)
{
     // copy edges and nodes
    m_edges = edges;
    m_nodes = nodes;
    m_projection = projection;

    return true;
};

bool GridGeom::Mesh::Administrate() 
{
    m_nodesEdges.resize(m_nodes.size(),
        std::vector<std::size_t>(maximumNumberOfEdgesPerNode, 0));
    m_nodesNumEdges.resize(m_nodes.size());
    m_edgesNumFaces.resize(m_edges.size());

    m_numFaces = 0;
    m_facesNodes.resize(0);
    m_facesEdges.resize(0);
    m_facesCircumcenters.resize(0);
    m_facesMassCenters.resize(0);
    m_faceArea.resize(0);

    m_edgesFaces.resize(m_edges.size(), std::vector<int>(2, -1));

    if (m_edges.size() == 0 || m_nodes.size() == 0)
    {
        return true;
    }

    // run administration and find the faces    
    NodeAdministration();
    SortEdgesInCounterClockWiseOrder();

    // find faces
    FindFaces();

    // find mesh circumcenters
    FaceCircumcenters(1.0);

    // compute faces areas and centers of mass
    FacesAreasAndMassCenters();

    // classify node types
    ClassifyNodes();

    // return value
    return true;


}

GridGeom::Mesh::Mesh(const CurvilinearGrid& curvilinearGrid) 
{
    if (curvilinearGrid.m_grid.size() == 0) 
    {
        return;
    }

    m_nodes.resize(curvilinearGrid.m_grid.size()*curvilinearGrid.m_grid[0].size());
    m_edges.resize(curvilinearGrid.m_grid.size() * (curvilinearGrid.m_grid[0].size() - 1) + (curvilinearGrid.m_grid.size() - 1) * curvilinearGrid.m_grid[0].size());
    std::vector<std::vector<int>> indexses(curvilinearGrid.m_grid.size(), std::vector<int>(curvilinearGrid.m_grid[0].size()));

    int ind = 0;
    for (int m = 0; m < curvilinearGrid.m_grid.size(); m++)
    {
        for (int n = 0; n < curvilinearGrid.m_grid[0].size(); n++)
        {
            m_nodes[ind] = curvilinearGrid.m_grid[m][n];
            indexses[m][n] = ind;
            ind++;
        }
    }

    ind = 0;
    for (int m = 0; m < curvilinearGrid.m_grid.size(); m++)
    {
        for (int n = 0; n < curvilinearGrid.m_grid[0].size() - 1; n++)
        {
            m_edges[ind].first = indexses[m][n];
            m_edges[ind].second = indexses[m][n + 1];
            ind++;
        }
    }

    for (int n = 0; n < curvilinearGrid.m_grid[0].size(); n++)
    {
        for (int m = 0; m < curvilinearGrid.m_grid.size()-1; m++)
        {
            m_edges[ind].first = indexses[m][n];
            m_edges[ind].second = indexses[m+1][n];
            ind++;
        }
    }

    Administrate();
}


bool GridGeom::Mesh::SetFlatCopies()
{
    // Used for internal state
    if (m_nodes.size() >0)
    {
        m_nodex.resize(m_nodes.size());
        m_nodey.resize(m_nodes.size());
        m_nodez.resize(m_nodes.size(), 0.0);
        for (int n = 0; n < m_nodex.size(); n++)
        {
            m_nodex[n] = m_nodes[n].x;
            m_nodey[n] = m_nodes[n].y;
        }
        m_edgeNodes.resize(m_edges.size() * 2);
        int ei = 0;
        for (int e = 0; e < m_edges.size(); e++)
        {
            m_edgeNodes[ei] = m_edges[e].first;
            ei++;
            m_edgeNodes[ei] = m_edges[e].second;
            ei++;
        }    
    }
    else
    {
        m_nodex.resize(1);
        m_nodey.resize(1);
        m_nodez.resize(1);
        m_edgeNodes.resize(1);    
    }

    return true;
}

bool GridGeom::Mesh::DeleteFlatCopies()
{
    //Used for internal state
    m_edges.resize(0);
    m_nodes.resize(0);
    m_nodesEdges.resize(0);
    m_nodesNumEdges.resize(0);
    m_edgesNumFaces.resize(0);
    m_edgesFaces.resize(0);
    m_facesNodes.resize(0);
    m_facesEdges.resize(0);
    m_facesCircumcenters.resize(0);
    m_facesMassCenters.resize(0);
    m_faceArea.resize(0);
    SetFlatCopies();

    return true;
}

void GridGeom::Mesh::NodeAdministration()
{
    // assume no duplicated links
    for (std::size_t e = 0; e < m_edges.size(); e++)
    {
        const std::size_t firstNode = m_edges[e].first;
        const std::size_t secondNode = m_edges[e].second;
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


void GridGeom::Mesh::SortEdgesInCounterClockWiseOrder()
{
    std::vector<double> edgesAngles(GridGeom::maximumNumberOfEdgesPerNode, 0.0);
    for (std::size_t node = 0; node < m_nodes.size(); node++)
    {
        double phi0 = 0.0;
        double phi;
        std::fill(edgesAngles.begin(), edgesAngles.end(), 0.0);
        for (std::size_t edgeIndex = 0; edgeIndex < m_nodesNumEdges[node]; edgeIndex++)
        {

            auto firstNode = m_edges[m_nodesEdges[node][edgeIndex]].first;
            auto secondNode = m_edges[m_nodesEdges[node][edgeIndex]].second;
            if (secondNode == node)
            {
                secondNode = firstNode;
                firstNode = node;
            }

            double deltaX = GetDx(m_nodes[secondNode], m_nodes[firstNode], m_projection);
            double deltaY = GetDy(m_nodes[secondNode], m_nodes[firstNode], m_projection);
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
        std::vector<std::size_t> indexes(m_nodesNumEdges[node]);
        std::vector<std::size_t> edgeNodeCopy{ m_nodesEdges[node] };
        iota(indexes.begin(), indexes.end(), 0);
        sort(indexes.begin(), indexes.end(), [&edgesAngles](std::size_t i1, std::size_t i2) {return edgesAngles[i1] < edgesAngles[i2]; });

        for (std::size_t edgeIndex = 0; edgeIndex < m_nodesNumEdges[node]; edgeIndex++)
        {
            m_nodesEdges[node][edgeIndex] = edgeNodeCopy[indexes[edgeIndex]];
        }
    }
}

// find cells
void GridGeom::Mesh::FindFaces(const int& numEdges)
{

    std::vector<std::size_t> foundEdges(numEdges);
    std::vector<std::size_t> foundNodes(numEdges); 

    for (std::size_t node = 0; node < m_nodes.size(); node++)
    {
        for (std::size_t firstEdgeLocalIndex = 0; firstEdgeLocalIndex < m_nodesNumEdges[node]; firstEdgeLocalIndex++)
        {
            std::size_t indexFoundNodes = 0;
            std::size_t indexFoundEdges = 0;

            std::size_t currentEdge = m_nodesEdges[node][firstEdgeLocalIndex];
            std::size_t currentNode = node;
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
                for (std::size_t localEdgeIndex = 0; localEdgeIndex < m_nodesNumEdges[currentNode]; localEdgeIndex++)
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
                for (std::size_t localEdgeIndex = 0; localEdgeIndex <= indexFoundEdges; localEdgeIndex++)
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
                for (std::size_t localEdgeIndex = 0; localEdgeIndex <= indexFoundEdges; localEdgeIndex++)
                {
                    if (m_edgesNumFaces[foundEdges[localEdgeIndex]] < 1)
                    {
                        allEdgesHaveAFace = false;
                        break;
                    }
                }

                bool isAnAlreadyFoundBoundaryFace = true;
                for (std::size_t localEdgeIndex = 0; localEdgeIndex <= indexFoundEdges - 1; localEdgeIndex++)
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
                for (std::size_t localEdgeIndex = 0; localEdgeIndex <= indexFoundEdges; localEdgeIndex++)
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

void GridGeom::Mesh::FindFaces()
{
    for (int numEdgesPerFace = 3; numEdgesPerFace <= maximumNumberOfEdgesPerFace; numEdgesPerFace++)
    {
        FindFaces(numEdgesPerFace);
    }
}


void GridGeom::Mesh::FaceCircumcenters(const double& weightCircumCenter)
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
        std::size_t numberOfFaceNodes = m_facesNodes[f].size();
        double xCenter = 0.0;
        double yCenter = 0.0;
        std::size_t numberOfInteriorEdges = 0;
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
        centerOfMass.x = xCenter / numberOfFaceNodes;
        centerOfMass.y = yCenter / numberOfFaceNodes;
        localFace[numberOfFaceNodes] = localFace[0];
        if (numberOfFaceNodes == 3)
        {
            circumcenterOfTriangle(m_nodes[m_facesNodes[f][0]], m_nodes[m_facesNodes[f][1]], m_nodes[m_facesNodes[f][2]], m_facesCircumcenters[f], m_projection);
        }
        else
        {
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
                    NormalVector(localFace[n], localFace[nextNode], middlePoints[n], normals[n], m_projection);
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
                            double dx = GetDx(middlePoints[n], estimatedCircumCenter, m_projection);
                            double dy = GetDy(middlePoints[n], estimatedCircumCenter, m_projection);
                            double increment = -0.1 * DotProduct(dx, dy, normals[n].x, normals[n].y);
                            Add(estimatedCircumCenter, normals[n], increment, m_projection);
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
            localFace[numberOfFaceNodes] = localFace[0];

            bool isCircumcenterInside = IsPointInPolygon(m_facesCircumcenters[f], localFace, numberOfFaceNodes);

            if (!isCircumcenterInside)
            {
                for (int n = 0; n < numberOfFaceNodes; n++)
                {
                    int nextNode = n + 1;
                    if (nextNode == numberOfFaceNodes) nextNode = 0;
                    Point intersection;
                    double crossProduct;
                    double firstRatio;
                    double secondRatio;
                    bool areLineCrossing = AreLinesCrossing(centerOfMass, m_facesCircumcenters[f], localFace[n], localFace[nextNode], false, intersection, crossProduct, firstRatio, secondRatio, m_projection);
                    if (areLineCrossing)
                    {
                        m_facesCircumcenters[f] = intersection;
                        break;
                    }
                }
            }
        }

    }
}

void GridGeom::Mesh::FacesAreasAndMassCenters()
{
    // polygon coordinates 
    m_faceArea.resize(m_facesNodes.size());
    m_facesMassCenters.resize(m_facesNodes.size());
    std::vector<Point> localFace(maximumNumberOfNodesPerFace + 1);
    Point centerOfMass;
    double area = 0.0;
    for (int f = 0; f < m_facesNodes.size(); f++)
    {
        std::size_t numberOfFaceNodes = m_facesNodes[f].size();

        for (int n = 0; n < numberOfFaceNodes; n++)
        {
            localFace[n] = m_nodes[m_facesNodes[f][n]];
        }
        localFace[numberOfFaceNodes] = localFace[0];
        area = 0.0;
        faceAreaAndCenterOfMass(localFace, numberOfFaceNodes, area, centerOfMass, m_projection);
        m_faceArea[f] = area;
        m_facesMassCenters[f].x = centerOfMass.x;
        m_facesMassCenters[f].y = centerOfMass.y;
    }
}

bool GridGeom::Mesh::ClassifyNodes()
{
    m_nodesTypes.resize(m_nodes.size(), 0);

    for (int e = 0; e < m_edges.size(); e++)
    {
        std::size_t first = m_edges[e].first;
        std::size_t second = m_edges[e].second;

        if (m_edgesNumFaces[e] == 0)
        {
            m_nodesTypes[first] = -1;
            m_nodesTypes[second] = -1;
        }
        else if (m_edgesNumFaces[e] == 1)
        {
            m_nodesTypes[first] += 1;
            m_nodesTypes[second] += 1;
        }
    }

    for (int n = 0; n < m_nodes.size(); n++)
    {
        if (m_nodesTypes[n] == 1 || m_nodesTypes[n] == 2)
        {
            if (m_nodesNumEdges[n] == 2)
            {
                //corner point
                m_nodesTypes[n] = 3;
            }
            else {}
        }
        else if (m_nodesTypes[n] > 2)
        {
            // corner point
            m_nodesTypes[n] = 3;
        }
        else if (m_nodesTypes[n] != -1)
        {
            //internal node
            m_nodesTypes[n] = 1;
        }

        if (m_nodesNumEdges[n] < 2)
        {
            //hanging node
            m_nodesTypes[n] = -1;
        }
    }
    return true;
}
