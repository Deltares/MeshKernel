// TODO: duplicated face detection issue for large number of edges
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>

#include "Mesh.hpp"
#include "Constants.cpp"
#include "Operations.cpp"
#include "Polygons.hpp"
#include <stack>

bool GridGeom::Mesh::Set(const std::vector<Edge>& edges, const std::vector<Point>& nodes, Projections projection)
{
    // copy edges and nodes
    m_edges = edges;
    m_nodes = nodes;
    m_projection = projection;

    Administrate();

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
//gridtonet
GridGeom::Mesh::Mesh(const CurvilinearGrid& curvilinearGrid, Projections projection)
{
    m_projection = projection;

    if (curvilinearGrid.m_grid.size() == 0)
    {
        return;
    }

    m_nodes.resize(curvilinearGrid.m_grid.size()*curvilinearGrid.m_grid[0].size());
    m_edges.resize(curvilinearGrid.m_grid.size() * (curvilinearGrid.m_grid[0].size() - 1) + (curvilinearGrid.m_grid.size() - 1) * curvilinearGrid.m_grid[0].size());
    std::vector<std::vector<int>> indexses(curvilinearGrid.m_grid.size(), std::vector<int>(curvilinearGrid.m_grid[0].size(), intMissingValue));

    int ind = 0;
    for (int m = 0; m < curvilinearGrid.m_grid.size(); m++)
    {
        for (int n = 0; n < curvilinearGrid.m_grid[0].size(); n++)
        {
            if (curvilinearGrid.m_grid[m][n].IsValid())
            {
                m_nodes[ind] = curvilinearGrid.m_grid[m][n];
                indexses[m][n] = ind;
                ind++;
            }
        }
    }
    m_nodes.resize(ind);

    ind = 0;
    for (int m = 0; m < curvilinearGrid.m_grid.size() - 1; m++)
    {
        for (int n = 0; n < curvilinearGrid.m_grid[0].size(); n++)
        {
            if (indexses[m][n] != intMissingValue && indexses[m + 1][n] != intMissingValue)
            {
                m_edges[ind].first = indexses[m][n];
                m_edges[ind].second = indexses[m + 1][n];
                ind++;
            }
        }
    }

    for (int m = 0; m < curvilinearGrid.m_grid.size(); m++)
    {
        for (int n = 0; n < curvilinearGrid.m_grid[0].size() - 1; n++)
        {
            if (indexses[m][n] != intMissingValue && indexses[m][n + 1] != intMissingValue)
            {
                m_edges[ind].first = indexses[m][n];
                m_edges[ind].second = indexses[m][n + 1];
                ind++;
            }
        }
    }
    m_edges.resize(ind);

    Administrate();
}

GridGeom::Mesh::Mesh(std::vector<Point>& inputNodes, const GridGeom::Polygons& polygons, Projections projection)
{
    m_projection = projection;
    std::vector<double> xLocalPolygon(inputNodes.size());
    std::vector<double> yLocalPolygon(inputNodes.size());
    for (int i = 0; i < inputNodes.size(); ++i)
    {
        xLocalPolygon[i] = inputNodes[i].x;
        yLocalPolygon[i] = inputNodes[i].y;
    }

    int numtri = -1;
    int jatri = 3;
    int numPointsIn = inputNodes.size();
    int numPointsOut = 0;
    int numberOfTriangles = numPointsIn * 6 + 10;
    double averageTriangleArea = 0.0;
    int numedge = 0;
    std::vector<int> faceNodesFlat;
    std::vector<int> edgeNodesFlat;
    std::vector<int> faceEdgesFlat;
    std::vector<double> xNodesFlat;
    std::vector<double> yNodesFlat;
    // if the number of estimated triangles is not sufficent, tricall must be repeated
    while (numtri < 0)
    {
        numtri = numberOfTriangles;
        faceNodesFlat.resize(numberOfTriangles * 3);
        edgeNodesFlat.resize(numberOfTriangles * 2);
        faceEdgesFlat.resize(numberOfTriangles * 3);
        xNodesFlat.resize(numberOfTriangles * 3, doubleMissingValue);
        yNodesFlat.resize(numberOfTriangles * 3, doubleMissingValue);
        Triangulation(&jatri,
            &xLocalPolygon[0],
            &yLocalPolygon[0],
            &numPointsIn,
            &faceNodesFlat[0],   // INDX
            &numtri,
            &edgeNodesFlat[0], // EDGEINDX
            &numedge,
            &faceEdgesFlat[0], // TRIEDGE
            &xNodesFlat[0],
            &yNodesFlat[0],
            &numPointsOut,
            &averageTriangleArea);
        if (numberOfTriangles)
        {
            numberOfTriangles = -numtri;
        }
    }

    // create face nodes
    std::vector<std::vector<int>> faceNodes(numtri, std::vector<int>(3, -1));
    std::vector<std::vector<int>> faceEdges(numtri, std::vector<int>(3, -1));
    int index = 0;
    for (int i = 0; i < numtri; ++i)
    {
        faceNodes[i][0] = faceNodesFlat[index] - 1;
        faceEdges[i][0] = faceEdgesFlat[index] - 1;
        index++;
        faceNodes[i][1] = faceNodesFlat[index] - 1;
        faceEdges[i][1] = faceEdgesFlat[index] - 1;
        index++;
        faceNodes[i][2] = faceNodesFlat[index] - 1;
        faceEdges[i][2] = faceEdgesFlat[index] - 1;
        index++;
    }

    // create edges
    std::vector<std::vector<int>> edgeNodes(numedge, std::vector<int>(2, 0));
    index = 0;
    for (int i = 0; i < numedge; ++i)
    {
        edgeNodes[i][0] = edgeNodesFlat[index] - 1;
        index++;
        edgeNodes[i][1] = edgeNodesFlat[index] - 1;
        index++;
    }


    // for each triangle we have to check
    // 1. validity of its internal angles
    // 2. is inside the polygon
    // if so we mark the edges and we add them to kn table
    std::vector<bool> edgeNodesFlag(numedge, false);
    for (int i = 0; i < numtri; ++i)
    {
        bool goodTriangle = CheckTriangle(faceNodes[i], inputNodes);

        if (!goodTriangle)
        {
            continue;
        }
        Point approximateCenter = (inputNodes[faceNodes[i][0]] + inputNodes[faceNodes[i][1]] + inputNodes[faceNodes[i][2]]) * oneThird;

        bool isTriangleInPolygon = IsPointInPolygon(approximateCenter, polygons.m_nodes, polygons.m_numNodes - 1);
        if (!isTriangleInPolygon)
        {
            continue;
        }

        // mark all edges of this triangle as good ones
        for (int j = 0; j < 3; ++j)
        {
            edgeNodesFlag[faceEdges[i][j]] = true;
        }
    }

    // now add all points and all valid edges
    m_nodes = inputNodes;
    int validEdges = 0;
    for (int i = 0; i < numedge; ++i)
    {
        if (!edgeNodesFlag[i])
            continue;
        validEdges++;
    }

    m_edges.resize(validEdges);
    validEdges = 0;
    for (int i = 0; i < numedge; ++i)
    {
        if (!edgeNodesFlag[i])
            continue;

        m_edges[validEdges].first = std::abs(edgeNodes[i][0]);
        m_edges[validEdges].second = edgeNodes[i][1];
        validEdges++;
    }

    Administrate();

}

bool GridGeom::Mesh::CheckTriangle(const std::vector<int>& faceNodes, const std::vector<Point>& nodes)
{
    double phiMin = 1e3;
    double phiMax = 0.0;
    static std::vector<std::vector<int>> nodePermutations
    {
        {2,0,1}, {0,1,2}, {1,2,0}
    };

    for (int i = 0; i < faceNodes.size(); ++i)
    {
        Point x0 = nodes[faceNodes[nodePermutations[i][0]]];
        Point x1 = nodes[faceNodes[nodePermutations[i][1]]];
        Point x2 = nodes[faceNodes[nodePermutations[i][2]]];

        double cosphi = NormalizedInnerProductTwoSegments(x1, x0, x1, x2, m_projection);
        double phi = std::acos(std::min(std::max(cosphi, -1.0), 1.0)) * raddeg_hp;
        phiMin = std::min(phiMin, phi);
        phiMax = std::max(phiMax, phi);
        if (phi < m_triangleMinimumAngle || phi > m_triangleMaximumAngle)
        {
            return false;
        }
    }
    return true;
}




bool GridGeom::Mesh::SetFlatCopies()
{
    // Used for internal state
    if (m_nodes.size() > 0)
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
    // you cold use std::sort + std::unique instead
    for (std::size_t e = 0; e < m_edges.size(); e++)
    {
        const std::size_t firstNode = m_edges[e].first;
        const std::size_t secondNode = m_edges[e].second;

        // Search for previously connected edges
        bool alreadyAddedEdge = false;
        for (int i = 0; i < m_nodesNumEdges[firstNode]; ++i)
        {
            auto currentEdge = m_edges[m_nodesEdges[firstNode][i]];
            if (currentEdge.first == secondNode || currentEdge.second == secondNode)
            {
                alreadyAddedEdge = true;
                break;
            }
        }
        if (!alreadyAddedEdge)
        {
            m_nodesEdges[firstNode][m_nodesNumEdges[firstNode]] = e;
            m_nodesNumEdges[firstNode]++;
        }

        // Search for previously connected edges
        alreadyAddedEdge = false;
        for (int i = 0; i < m_nodesNumEdges[secondNode]; ++i)
        {
            auto currentEdge = m_edges[m_nodesEdges[secondNode][i]];
            if (currentEdge.first == firstNode || currentEdge.second == firstNode)
            {
                alreadyAddedEdge = true;
                break;
            }
        }
        if (!alreadyAddedEdge)
        {
            m_nodesEdges[secondNode][m_nodesNumEdges[secondNode]] = e;
            m_nodesNumEdges[secondNode]++;
        }
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

bool GridGeom::Mesh::FindFacesRecursive(
    int startingNode, 
    int node, 
    int index, 
    int previusEdge,
    std::vector<size_t>& edges,
    std::vector<size_t>& nodes,
    std::vector<size_t>& sortedEdgesFaces,
    std::vector<size_t>& sortedNodes)
{
    if (index >= edges.size())
    {
        return false;
    }

    if (m_edgesNumFaces[previusEdge] >= 2)
    {
        return false;
    }

    edges[index] = previusEdge;
    nodes[index] = node;
    const int otherNode = m_edges[previusEdge].first + m_edges[previusEdge].second - node;

    // enclosure found
    if (otherNode == startingNode && index == edges.size() - 1)
    {
        // all nodes must be unique
        sortedNodes = nodes;
        std::sort(sortedNodes.begin(), sortedNodes.end());
        for (int n = 0; n < sortedNodes.size() - 1; n++)
        {
            if (sortedNodes[n + 1] == sortedNodes[n])
            {
                return false;
            }
        }
        // we need to add a face when at least one edge has no faces
        bool oneEdgeHasNoFace = false;
        for (int ee = 0; ee < edges.size(); ee++)
        {
            if (m_edgesNumFaces[edges[ee]] == 0)
            {
                oneEdgeHasNoFace = true;
                break;
            }
        }

        if (!oneEdgeHasNoFace)
        {
            // is an internal face only if all edges have a different face
            for (int ee = 0; ee < edges.size(); ee++)
            {
                sortedEdgesFaces[ee] = m_edgesFaces[edges[ee]][0];
            }
            std::sort(sortedEdgesFaces.begin(), sortedEdgesFaces.end());
            for (int n = 0; n < sortedEdgesFaces.size() - 1; n++)
            {
                if (sortedEdgesFaces[n + 1] == sortedEdgesFaces[n])
                {
                    return false;
                }
            }
        }

        // increase m_edgesNumFaces 
        m_numFaces += 1;
        for (int ee = 0; ee < edges.size(); ee++)
        {
            m_edgesNumFaces[edges[ee]] += 1;
            const int numFace = m_edgesNumFaces[edges[ee]];
            m_edgesFaces[edges[ee]][numFace - 1] = m_numFaces - 1;
        }

        // store the result
        m_facesNodes.push_back(nodes);
        m_facesEdges.push_back(edges);
        return true;
    }
    else
    {
        int edgeIndexOtherNode = 0;
        for (int e = 0; e < m_nodesNumEdges[otherNode]; e++)
        {
            if (m_nodesEdges[otherNode][e] == previusEdge)
            {
                edgeIndexOtherNode = e;
                break;
            }
        }

        edgeIndexOtherNode = edgeIndexOtherNode - 1;
        if (edgeIndexOtherNode < 0)
        {
            edgeIndexOtherNode = edgeIndexOtherNode + m_nodesNumEdges[otherNode];
        }
        if (edgeIndexOtherNode > m_nodesNumEdges[otherNode] - 1)
        {
            edgeIndexOtherNode = edgeIndexOtherNode - m_nodesNumEdges[otherNode];
        }

        const int edge = m_nodesEdges[otherNode][edgeIndexOtherNode];
        FindFacesRecursive(startingNode, otherNode, index + 1, edge, edges, nodes, sortedEdgesFaces, sortedNodes);
    }

}

void GridGeom::Mesh::FindFaces()
{
    for (int numEdgesPerFace = 3; numEdgesPerFace <= maximumNumberOfEdgesPerFace; numEdgesPerFace++)
    {
        std::vector<size_t> edges(numEdgesPerFace);
        std::vector<size_t> nodes(numEdgesPerFace);
        std::vector<size_t> sortedEdgesFaces(numEdgesPerFace);
        std::vector<size_t> sortedNodes(numEdgesPerFace);
        for (int n = 0; n < m_nodes.size(); n++)
        {
            for (int e = 0; e < m_nodesNumEdges[n]; e++)
            {
                FindFacesRecursive(n, n, 0, m_nodesEdges[n][e], edges, nodes, sortedEdgesFaces, sortedNodes);
            }
        }
    }
}

// Find cells, needs to be corrected, for example look at sub_findelemcontours in IrregularGridClass.f90
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


bool GridGeom::Mesh::MakeMesh(const GridGeomApi::MakeGridParametersNative& makeGridParametersNative, const Polygons& polygons)
{
    CurvilinearGrid CurvilinearGrid;
    m_projection = polygons.m_projection;
    if (makeGridParametersNative.GridType == 0)
    {
        // regular grid
        int numM = makeGridParametersNative.NumberOfRows + 1;
        int numN = makeGridParametersNative.NumberOfColumns + 1;
        double XGridBlockSize = makeGridParametersNative.XGridBlockSize;
        double YGridBlockSize = makeGridParametersNative.YGridBlockSize;
        double cosineAngle = std::cos(makeGridParametersNative.GridAngle*degrad_hp);
        double sinAngle = std::sin(makeGridParametersNative.GridAngle*degrad_hp);
        double OriginXCoordinate = makeGridParametersNative.OriginXCoordinate;
        double OriginYCoordinate = makeGridParametersNative.OriginYCoordinate;

        // in case a polygon is there, re-compute parameters
        if (polygons.m_numNodes >= 3)
        {
            Point referencePoint;
            // rectangular grid in polygon
            for (int i = 0; i < polygons.m_numNodes; ++i)
            {
                if (polygons.m_nodes[i].IsValid())
                {
                    referencePoint = polygons.m_nodes[i];
                    break;
                }
            }

            // get polygon min/max in rotated (xi,eta) coordinates
            double xmin = std::numeric_limits<double>::max();
            double xmax = -xmin;
            double etamin = std::numeric_limits<double>::max();
            double etamax = -etamin;
            for (int i = 0; i < polygons.m_numNodes; ++i)
            {
                if (polygons.m_nodes[i].IsValid())
                {
                    double dx = GetDx(referencePoint, polygons.m_nodes[i], m_projection);
                    double dy = GetDy(referencePoint, polygons.m_nodes[i], m_projection);
                    double xi = dx * cosineAngle + dy * sinAngle;
                    double eta = -dx * sinAngle + dy * cosineAngle;
                    xmin = std::min(xmin, xi);
                    xmax = std::max(xmax, xi);
                    etamin = std::min(etamin, eta);
                    etamax = std::max(etamax, eta);
                }
            }

            double xShift = xmin*cosineAngle - etamin * sinAngle;
            double yShift = xmin*sinAngle + etamin* cosineAngle;
            if (m_projection == Projections::spherical)
            {
                xShift = xShift / earth_radius *raddeg_hp;
                yShift = yShift / (earth_radius *std::cos(referencePoint.y*degrad_hp)) * raddeg_hp;
            }

            OriginXCoordinate = referencePoint.x + xShift;
            OriginYCoordinate = referencePoint.y + yShift;
            numN = std::ceil((etamax - etamin) / XGridBlockSize) + 1;
            numM = std::ceil((xmax - xmin) / YGridBlockSize) + 1;
        }


        CurvilinearGrid.IncreaseGrid(numN, numM);
        for (int n = 0; n < numN; ++n)
        {
            for (int m = 0; m < numM; ++m)
            {
                CurvilinearGrid.m_grid[n][m] = Point
                {
                    OriginXCoordinate + m*XGridBlockSize * cosineAngle - n * YGridBlockSize * sinAngle,
                    OriginYCoordinate + m*XGridBlockSize * sinAngle + n * YGridBlockSize * cosineAngle
                };
            }
        }

        // in case a polygon is there, remove nodes outside
        if (polygons.m_numNodes >= 3)
        {
            std::vector<std::vector<bool>> nodeBasedMask(numN, std::vector<bool>(numM, false));
            std::vector<std::vector<bool>> faceBasedMask(numN - 1, std::vector<bool>(numM - 1, false));
            // mark points inside a polygon
            for (int n = 0; n < numN; ++n)
            {
                for (int m = 0; m < numM; ++m)
                {
                    bool isInPolygon = IsPointInPolygon(CurvilinearGrid.m_grid[n][m], polygons.m_nodes, polygons.m_numNodes - 1);
                    if (isInPolygon)
                    {
                        nodeBasedMask[n][m] = true;
                    }
                }
            }

            // mark faces when at least one node is inside
            for (int n = 0; n < numN - 1; ++n)
            {
                for (int m = 0; m < numM - 1; ++m)
                {
                    if (nodeBasedMask[n][m] || nodeBasedMask[n + 1][m] || nodeBasedMask[n][m + 1] || nodeBasedMask[n + 1][m + 1])
                    {
                        faceBasedMask[n][m] = true;
                    }
                }
            }

            //mark nodes that are member of a cell inside the polygon(s)
            for (int n = 0; n < numN - 1; ++n)
            {
                for (int m = 0; m < numM - 1; ++m)
                {
                    if (faceBasedMask[n][m])
                    {
                        nodeBasedMask[n][m] = true;
                        nodeBasedMask[n + 1][m] = true;
                        nodeBasedMask[n][m + 1] = true;
                        nodeBasedMask[n + 1][m + 1] = true;
                    }
                }
            }

            // mark points inside a polygon
            for (int n = 0; n < numN; ++n)
            {
                for (int m = 0; m < numM; ++m)
                {
                    if (!nodeBasedMask[n][m])
                    {
                        CurvilinearGrid.m_grid[n][m].x = doubleMissingValue;
                        CurvilinearGrid.m_grid[n][m].y = doubleMissingValue;
                    }
                }
            }
        }
    }

    // Assign mesh
    *this = Mesh(CurvilinearGrid, m_projection);

    return true;
}