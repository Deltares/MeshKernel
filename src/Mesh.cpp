#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>

#include "Mesh.hpp"
#include "Constants.cpp"
#include "Operations.cpp"
#include "Polygons.hpp"
#include "SpatialTrees.hpp"
#include "CurvilinearGrid.hpp"
#include "Entities.hpp"
#include "MakeGridParametersNative.hpp"

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
    //Nodes, edges and faces indexses. Shift valid values up
    m_nodeMask.resize(m_nodes.size());
    std::fill(m_nodeMask.begin(), m_nodeMask.end(), 1);
    int validIndex = 0;
    for (int n = 0; n < m_nodes.size(); ++n)
    {
        if(m_nodes[n].IsValid())
        {
            m_nodeMask[n] = validIndex;
            validIndex++;
        }
    }

    auto endNodeVector = std::remove_if(m_nodes.begin(), m_nodes.end(), [](const Point& n) {return !n.IsValid(); });
    m_numNodes = endNodeVector - m_nodes.begin();

    auto endEdgeVector = std::remove_if(m_edges.begin(), m_edges.end(), [](const Edge& e) {return e.first < 0 || e.second < 0; });
    m_numEdges = endEdgeVector - m_edges.begin();

    if (m_numNodes == 0 || m_numEdges == 0)
    {
        return true;
    }

    for (int e = 0; e < m_numEdges; ++e)
    {
        m_edges[e].first = m_nodeMask[m_edges[e].first];
        m_edges[e].second = m_nodeMask[m_edges[e].second];
    }

    m_numFaces = 0;
    std::fill(m_nodeMask.begin(), m_nodeMask.end(), 1);

    ResizeVectorIfNeeded(m_nodes.size(), m_nodesEdges);
    std::fill(m_nodesEdges.begin(), m_nodesEdges.end(), std::vector<int>(maximumNumberOfEdgesPerNode, 0));

    ResizeVectorIfNeeded(m_nodes.size(), m_nodesNumEdges);
    std::fill(m_nodesNumEdges.begin(), m_nodesNumEdges.end(), 0);

    ResizeVectorIfNeeded(m_edges.size(), m_edgesNumFaces);
    std::fill(m_edgesNumFaces.begin(), m_edgesNumFaces.end(), 0);

    ResizeVectorIfNeeded(m_edges.size(), m_edgesFaces);
    std::fill(m_edgesFaces.begin(), m_edgesFaces.end(), std::vector<int>(2, -1));

    m_facesNodes.resize(0);
    m_facesEdges.resize(0);
    m_facesCircumcenters.resize(0);
    m_facesMassCenters.resize(0);
    m_faceArea.resize(0);

    // run administration and find the faces    
    NodeAdministration();

    SortEdgesInCounterClockWiseOrder();

    // find faces
    FindFaces();

    // find mesh circumcenters
    ComputeFaceCircumcentersMassCentersAreas();

    // classify node types
    ClassifyNodes();

    m_isAdministrationDone = true;

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
    if (GetNumNodes() > 0)
    {
        m_nodex.resize(GetNumNodes());
        m_nodey.resize(GetNumNodes());
        m_nodez.resize(GetNumNodes(), 0.0);
        for (int n = 0; n < GetNumNodes(); n++)
        {
            m_nodex[n] = m_nodes[n].x;
            m_nodey[n] = m_nodes[n].y;
        }
        m_edgeNodes.resize(GetNumEdges() * 2);
        int ei = 0;
        for (int e = 0; e < GetNumEdges(); e++)
        {
            m_edgeNodes[ei] = m_edges[e].first;
            ei++;
            m_edgeNodes[ei] = m_edges[e].second;
            ei++;
        }
    }
    else
    {
        DeleteFlatCopies();
    }

    return true;
}

bool GridGeom::Mesh::DeleteMesh()
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

    m_isAdministrationDone = false;

    return true;
}

bool GridGeom::Mesh::DeleteFlatCopies()
{
    //Used for internal state
    m_nodex.resize(1);
    m_nodey.resize(1);
    m_nodez.resize(1);
    m_edgeNodes.resize(1);

    return true;
}

void GridGeom::Mesh::NodeAdministration()
{
    // assume no duplicated links
    // you cold use std::sort + std::unique instead
    for (int e = 0; e < GetNumEdges(); e++)
    {
        const auto firstNode = m_edges[e].first;
        const auto secondNode = m_edges[e].second;

        if (firstNode < 0 || secondNode < 0)
        {
            continue;
        }

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
    for (auto n = 0; n < GetNumNodes(); n++)
    {
        m_nodesEdges[n].resize(m_nodesNumEdges[n]);
    }
};


void GridGeom::Mesh::SortEdgesInCounterClockWiseOrder()
{
    std::vector<double> edgesAngles(GridGeom::maximumNumberOfEdgesPerNode, 0.0);
    for (auto node = 0; node < GetNumNodes(); node++)
    {
        if (!m_nodes[node].IsValid())
        {
            continue;
        }

        double phi0 = 0.0;
        double phi;
        std::fill(edgesAngles.begin(), edgesAngles.end(), 0.0);
        for (auto edgeIndex = 0; edgeIndex < m_nodesNumEdges[node]; edgeIndex++)
        {

            auto firstNode = m_edges[m_nodesEdges[node][edgeIndex]].first;
            auto secondNode = m_edges[m_nodesEdges[node][edgeIndex]].second;
            if (firstNode < 0 || secondNode < 0)
            {
                continue;
            }

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
        std::vector<int> edgeNodeCopy{ m_nodesEdges[node] };
        iota(indexes.begin(), indexes.end(), 0);
        sort(indexes.begin(), indexes.end(), [&edgesAngles](std::size_t i1, std::size_t i2) {return edgesAngles[i1] < edgesAngles[i2]; });

        for (std::size_t edgeIndex = 0; edgeIndex < m_nodesNumEdges[node]; edgeIndex++)
        {
            m_nodesEdges[node][edgeIndex] = edgeNodeCopy[indexes[edgeIndex]];
        }
    }
}

// look at sub_findelemcontours in IrregularGridClass.f90 for a similar implementation
bool GridGeom::Mesh::FindFacesRecursive(
    int startingNode,
    int node,
    int index,
    int previusEdge,
    std::vector<int>& edges,
    std::vector<int>& nodes,
    std::vector<int>& sortedEdgesFaces,
    std::vector<int>& sortedNodes)
{
    if (index >= edges.size())
    {
        return false;
    }

    if (m_edgesNumFaces[previusEdge] >= 2)
    {
        return false;
    }

    if (m_edges[previusEdge].first < 0 || m_edges[previusEdge].second < 0)
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


    return true;

}

void GridGeom::Mesh::FindFaces()
{
    for (int numEdgesPerFace = 3; numEdgesPerFace <= maximumNumberOfEdgesPerFace; numEdgesPerFace++)
    {
        std::vector<int> edges(numEdgesPerFace);
        std::vector<int> nodes(numEdgesPerFace);
        std::vector<int> sortedEdgesFaces(numEdgesPerFace);
        std::vector<int> sortedNodes(numEdgesPerFace);
        for (int n = 0; n < GetNumNodes(); n++)
        {
            if (!m_nodes[n].IsValid())
            {
                continue;
            }

            for (int e = 0; e < m_nodesNumEdges[n]; e++)
            {
                FindFacesRecursive(n, n, 0, m_nodesEdges[n][e], edges, nodes, sortedEdgesFaces, sortedNodes);
            }
        }
    }

    m_numFacesNodes.resize(m_numFaces);
    for (int f = 0; f < m_numFaces; ++f)
    {
        m_numFacesNodes[f] = m_facesNodes[f].size();
    }
}

void GridGeom::Mesh::ComputeFaceCircumcentersMassCentersAreas()
{
    m_facesCircumcenters.resize(GetNumFaces());
    m_faceArea.resize(GetNumFaces());
    m_facesMassCenters.resize(GetNumFaces());
    Point centerOfMass;
    const int maximumNumberCircumcenterIterations = 100;
    std::vector<Point> middlePointsCache(maximumNumberOfNodesPerFace);
    std::vector<Point> normalsCache(maximumNumberOfNodesPerFace);
    std::vector<int> numEdgeFacesCache(maximumNumberOfEdgesPerFace);
    m_polygonNodesCache.resize(maximumNumberOfNodesPerFace + 1);
    for (int f = 0; f < GetNumFaces(); f++)
    {
        //need to account for spherical coordinates. Build a polygon around a face
        int numPolygonPoints;
        bool successful = FacePolygon(f, m_polygonNodesCache, numPolygonPoints);
        if (!successful)
        {
            return;
        }

        auto numberOfFaceNodes = GetNumFaceEdges(f);
        double area;
        Point centerOfMass;
        successful = FaceAreaAndCenterOfMass(m_polygonNodesCache, numberOfFaceNodes, m_projection, area, centerOfMass);
        if (!successful)
        {
            return;
        }

        m_faceArea[f] = area;
        m_facesMassCenters[f] = centerOfMass;

        int numberOfInteriorEdges = 0;
        for (int n = 0; n < numberOfFaceNodes; n++)
        {
            if (m_edgesNumFaces[m_facesEdges[f][n]] == 2)
            {
                numberOfInteriorEdges += 1;
            }
        }
        if (numberOfInteriorEdges == 0)
        {
            m_facesCircumcenters[f] = centerOfMass;
            continue;
        }

        for (int n = 0; n < numberOfFaceNodes; n++)
        {
            numEdgeFacesCache[n] = m_edgesNumFaces[m_facesEdges[f][n]];
        }

        successful = ComputePolygonCircumenter(m_polygonNodesCache,
            middlePointsCache,
            normalsCache,
            numberOfFaceNodes,
            numEdgeFacesCache,
            m_projection,
            weightCircumCenter,
            m_facesCircumcenters[f]);

        if (!successful)
        {
            return;
        }
    }
}

bool GridGeom::Mesh::ClassifyNodes()
{
    m_nodesTypes.resize(GetNumNodes(), 0);

    for (int e = 0; e < GetNumEdges(); e++)
    {
        const auto firstNode = m_edges[e].first;
        const auto secondNode = m_edges[e].second;

        if (firstNode < 0 || secondNode < 0)
        {
            continue;
        }

        if (m_edgesNumFaces[e] == 0)
        {
            m_nodesTypes[firstNode] = -1;
            m_nodesTypes[secondNode] = -1;
        }
        else if (m_edgesNumFaces[e] == 1)
        {
            m_nodesTypes[firstNode] += 1;
            m_nodesTypes[secondNode] += 1;
        }
    }

    for (int n = 0; n < GetNumNodes(); n++)
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

///MERGENODESINPOLYGON
bool GridGeom::Mesh::MergeNodesInPolygon(const Polygons& polygon)
{
    // first filter the nodes in polygon
    std::vector<Point> filteredPoints(GetNumNodes());
    int index = 0;
    for (int i = 0; i < GetNumNodes(); i++)
    {
        bool inPolygon = IsPointInPolygon(m_nodes[i], polygon.m_nodes, polygon.m_numNodes - 1);
        if (inPolygon)
        {
            filteredPoints[index] = m_nodes[i];
            index++;
        }
    }
    filteredPoints.resize(index);

    // build the R-Tree
    GridGeom::SpatialTrees::RTree rtree;
    rtree.BuildTree(filteredPoints, m_projection);

    // merge the closest nodes
    for (int i = 0; i < filteredPoints.size(); i++)
    {
        auto result = rtree.NearestNeighbours(filteredPoints[i], mergingDistance);

        if (result.size() > 1)
        {
            for (int j = 0; j < result.size(); j++)
            {
                if (result[j] != i)
                {
                    rtree.RemoveNode(result[j]);
                    MergeTwoNodes(i, result[j]);
                }
            }
        }
    }

    Administrate();

    return true;
}

///mergenodes
bool GridGeom::Mesh::MergeTwoNodes(int firstNodeIndex, int secondNodeIndex)
{
    int edgeIndex;
    FindEdge(firstNodeIndex, secondNodeIndex, edgeIndex);
    if (edgeIndex >= 0)
    {
        m_edges[edgeIndex].first = -1;
        m_edges[edgeIndex].second = -1;
        m_numEdges--;
    }

    // check if there is another edge starting at firstEdgeOtherNode and ending at secondNode
    for (auto n = 0; n < m_nodesNumEdges[firstNodeIndex]; n++)
    {
        auto firstEdgeIndex = m_nodesEdges[firstNodeIndex][n];
        auto firstEdge = m_edges[firstEdgeIndex];
        auto firstEdgeOtherNode = firstEdge.first + firstEdge.second - firstNodeIndex;
        if (firstEdgeOtherNode >= 0 && firstEdgeOtherNode != secondNodeIndex)
        {
            for (int nn = 0; nn < m_nodesNumEdges[firstEdgeOtherNode]; nn++)
            {
                auto secondEdgeIndex = m_nodesEdges[firstEdgeOtherNode][nn];
                auto secondEdge = m_edges[secondEdgeIndex];
                auto secondEdgeOtherNode = secondEdge.first + secondEdge.second - firstEdgeOtherNode;
                if (secondEdgeOtherNode == secondNodeIndex)
                {
                    m_edges[secondEdgeIndex].first = -1;
                    m_edges[secondEdgeIndex].second = -1;
                    m_numEdges--;
                }
            }
        }
    }

    // add all valid edges starting at secondNode
    std::vector<int> secondNodeEdges(maximumNumberOfEdgesPerNode);
    int numSecondNodeEdges = 0;
    for (auto n = 0; n < m_nodesNumEdges[secondNodeIndex]; n++)
    {
        auto edgeIndex = m_nodesEdges[secondNodeIndex][n];
        if (m_edges[edgeIndex].first >= 0)
        {
            secondNodeEdges[numSecondNodeEdges] = edgeIndex;
            numSecondNodeEdges++;
        }
    }

    // add all valid edges starting at firstNode
    for (auto n = 0; n < m_nodesNumEdges[firstNodeIndex]; n++)
    {
        auto edgeIndex = m_nodesEdges[firstNodeIndex][n];
        if (m_edges[edgeIndex].first >= 0)
        {
            secondNodeEdges[numSecondNodeEdges] = edgeIndex;
            if (m_edges[edgeIndex].first == firstNodeIndex)
            {
                m_edges[edgeIndex].first = secondNodeIndex;
            }
            if (m_edges[edgeIndex].second == firstNodeIndex)
            {
                m_edges[edgeIndex].second = secondNodeIndex;
            }
            numSecondNodeEdges++;
        }
    }

    // re-assign edges to second node
    m_nodesEdges[secondNodeIndex] = std::move(std::vector<int>(secondNodeEdges.begin(), secondNodeEdges.begin() + numSecondNodeEdges));
    m_nodesNumEdges[secondNodeIndex] = numSecondNodeEdges;

    // remove edges to first node
    //m_nodesEdges[firstNodeIndex] = std::move(std::vector<int>(0));
    //m_nodesNumEdges[firstNodeIndex] = 0;
    //m_nodes[firstNodeIndex] = { doubleMissingValue, doubleMissingValue };
    //m_numNodes--;
    DeleteNode(firstNodeIndex);

    m_isAdministrationDone = false;

    return true;
}

bool GridGeom::Mesh::IsAdministrationDone()  const
{
    return m_isAdministrationDone;
}

bool GridGeom::Mesh::ConnectNodes(int startNode, int endNode, int& newEdgeIndex)
{
    int edgeIndex;
    bool successful = FindEdge(startNode, endNode, edgeIndex);
    if (!successful)
    {
        return false;
    }
    if (edgeIndex >= 0)
    {
        return true;
    }

    // increment the edges container
    newEdgeIndex = GetNumEdges();
    ResizeVectorIfNeeded(newEdgeIndex + 1, m_edges, std::make_pair(intMissingValue,intMissingValue));
    m_edges[newEdgeIndex].first = startNode;
    m_edges[newEdgeIndex].second = endNode;
    m_numEdges++;

    // add the new edge to the nodes
    if (m_nodesNumEdges[startNode] + 1 > maximumNumberOfEdgesPerNode)
    {
        return false;
    }
    m_nodesNumEdges[startNode] = m_nodesNumEdges[startNode] + 1;
    m_nodesEdges[startNode].push_back(newEdgeIndex);


    if (m_nodesNumEdges[endNode] + 1 > maximumNumberOfEdgesPerNode)
    {
        return false;
    }
    m_nodesNumEdges[endNode] = m_nodesNumEdges[endNode] + 1;
    m_nodesEdges[endNode].push_back(newEdgeIndex);

    m_isAdministrationDone = false;

    return true;
}

bool GridGeom::Mesh::InsertNode(const Point& newPoint, int& newNodeIndex)
{
    int newSize = GetNumNodes() + 1;
    newNodeIndex = GetNumNodes();

    ResizeVectorIfNeeded(newSize, m_nodes);
    ResizeVectorIfNeeded(newSize, m_nodeMask);
    ResizeVectorIfNeeded(newSize, m_nodesNumEdges);
    ResizeVectorIfNeeded(newSize, m_nodesEdges);
    m_numNodes++;

    m_nodes[newNodeIndex] = newPoint;
    m_nodeMask[newNodeIndex] = newNodeIndex;
    m_nodesNumEdges[newNodeIndex] = 0;

    m_isAdministrationDone = false;

    return true;
}

bool GridGeom::Mesh::DeleteNode(int nodeIndex)
{
    //TODO
    for (int e = 0; e <  m_nodesNumEdges[nodeIndex]; e++)
    {
        auto edgeIndex = m_nodesEdges[nodeIndex][e];
        auto startNode = m_edges[edgeIndex].first;
        auto endNode = m_edges[edgeIndex].second;
        DeleteEdge(startNode, endNode);
    }
    m_nodesNumEdges[nodeIndex] = 0;
    m_nodes[nodeIndex] = { doubleMissingValue,doubleMissingValue };
    m_numNodes--;

    m_isAdministrationDone = false;

    return true;
}

//to do: check if accessing after remove if works
bool GridGeom::Mesh::DeleteEdge(int startNode, int endNode)
{
    if (startNode < 0 || endNode < 0)
    {
        return true; 
    }

    auto isNotValidIndex = [](const int& v) { return v < 0; };
    for (int e = 0; e <  m_nodesNumEdges[startNode]; e++)
    {
        auto firstEdgeIndex = m_nodesEdges[startNode][e];
        for (int ee = 0; ee < m_nodesNumEdges[endNode]; ee++)
        {
            auto secondEdgeIndex = m_nodesEdges[endNode][ee];
            if (firstEdgeIndex == secondEdgeIndex)
            {
                m_nodesEdges[startNode][e] = -1;
                m_nodesEdges[endNode][ee] = -1;

                std::remove_if(m_nodesEdges[startNode].begin(), m_nodesEdges[startNode].end(), isNotValidIndex);
                std::remove_if(m_nodesEdges[endNode].begin(), m_nodesEdges[endNode].end(), isNotValidIndex);
                m_nodesNumEdges[startNode]--;
                m_nodesNumEdges[endNode]--;
                m_numEdges--;
            }
        }
    }

    // remove startNode
    if (m_nodesNumEdges[startNode] == 0)
    {
        m_nodes[startNode] = { doubleMissingValue, doubleMissingValue };
        m_numNodes--;
    }

    // remove endNode
    if (m_nodesNumEdges[endNode] == 0)
    {
        m_nodes[endNode] = { doubleMissingValue, doubleMissingValue };
        m_numNodes--;
    }

    m_isAdministrationDone = false;

    return true;
}


bool GridGeom::Mesh::DeleteEdge(int edgeIndex)
{
    if(edgeIndex<0)
    {
        return true;
    }

    m_edges[edgeIndex].first = intMissingValue;
    m_edges[edgeIndex].second = intMissingValue;
    m_edgesNumFaces[edgeIndex] = 0;
    m_numEdges--;

    m_isAdministrationDone = false;

    return false;
}


bool GridGeom::Mesh::FacePolygon(int faceIndex, std::vector<Point>& polygonNodesCache, int& numPolygonPoints) const
{
    auto numFaceNodes = GetNumFaceEdges(faceIndex);
    if (polygonNodesCache.size() < numFaceNodes + 1)
    {
        polygonNodesCache.resize(numFaceNodes + 1);
    }

    for (int n = 0; n < numFaceNodes; n++)
    {
        polygonNodesCache[n] = m_nodes[m_facesNodes[faceIndex][n]];
    }
    polygonNodesCache[numFaceNodes] = polygonNodesCache[0];

    numPolygonPoints = numFaceNodes + 1;

    return true;
}


bool GridGeom::Mesh::FacePolygon(int faceIndex, std::vector<Point>& polygonNodesCache, std::vector<int>& localNodeIndexsesCache, std::vector<int>& edgeIndexsesCache) const
{
    auto numFaceNodes = GetNumFaceEdges(faceIndex);
    if (polygonNodesCache.size() < numFaceNodes + 1)
    {
        polygonNodesCache.resize(numFaceNodes + 1);
    }

    if (localNodeIndexsesCache.size() < numFaceNodes + 1)
    {
        localNodeIndexsesCache.resize(numFaceNodes + 1);
    }

    if (edgeIndexsesCache.size() < numFaceNodes + 1)
    {
        edgeIndexsesCache.resize(numFaceNodes + 1);
    }

    for (int n = 0; n < numFaceNodes; n++)
    {
        polygonNodesCache[n] = m_nodes[m_facesNodes[faceIndex][n]];
        localNodeIndexsesCache[n] = n;
        edgeIndexsesCache[n] = m_facesEdges[faceIndex][n];
    }
    polygonNodesCache[numFaceNodes] = polygonNodesCache[0];
    localNodeIndexsesCache[numFaceNodes] = 0;
    edgeIndexsesCache[numFaceNodes] = m_facesEdges[faceIndex][0];

    return true;
}


bool GridGeom::Mesh::SelectNodesInPolygon(const Polygons& polygon, int inside)
{

    int numberOfMeshVertices = 0;
    for (int i = 0; i < GetNumNodes(); ++i)
    {
        bool isInPolygon = IsPointInPolygon(m_nodes[i], polygon.m_nodes, polygon.m_numNodes - 1);
        if (inside == false)
        {
            isInPolygon = !isInPolygon;
        }
        if (isInPolygon)
        {
            m_nodeMask[i] = 1;
        }
    }

    return true;
}

bool GridGeom::Mesh::ComputeEdgeLengths()
{
    auto numEdges = GetNumEdges();
    m_edgeLengths.resize(numEdges, doubleMissingValue);
    for (int e = 0; e < numEdges; e++)
    {
        int first = m_edges[e].first;
        int second = m_edges[e].second;
        m_edgeLengths[e] = Distance(m_nodes[first], m_nodes[second], m_projection);
    }
    return true;
}

bool GridGeom::Mesh::IsFullFaceNotInPolygon(int faceIndex) const
{
    for (int n = 0; n < GetNumFaceEdges(faceIndex); n++)
    {
        if (m_nodeMask[m_facesNodes[faceIndex][n]] != 1)
        {
            return true;
        }
    }
    return false;
}

bool GridGeom::Mesh::FindCommonNode(int firstEdgeIndex, int secondEdgeIndex, int& node) const
{
    auto firstEdgeFirstNode = m_edges[firstEdgeIndex].first;
    auto firstEdgeEdgeSecondNode = m_edges[firstEdgeIndex].second;

    auto secondEdgeFirstNode = m_edges[secondEdgeIndex].first;
    auto secondEdgeSecondNode = m_edges[secondEdgeIndex].second;

    if (firstEdgeFirstNode < 0 || firstEdgeEdgeSecondNode < 0 || secondEdgeFirstNode < 0 || secondEdgeSecondNode < 0)
    {
        return false;
    }

    if (firstEdgeFirstNode == secondEdgeFirstNode || firstEdgeFirstNode == secondEdgeSecondNode)
    {
        node = firstEdgeFirstNode;
        return true;
    }

    if (firstEdgeEdgeSecondNode == secondEdgeFirstNode || firstEdgeEdgeSecondNode == secondEdgeSecondNode)
    {
        node = firstEdgeEdgeSecondNode;
        return true;
    }

    return true;
}

bool GridGeom::Mesh::FindEdge(int firstNodeIndex, int secondNodeIndex, int& edgeIndex) const
{
    if (firstNodeIndex < 0 || secondNodeIndex < 0)
    {
        return false;
    }

    edgeIndex = -1;
    for (auto n = 0; n < m_nodesNumEdges[firstNodeIndex]; n++)
    {
        int localEdgeIndex = m_nodesEdges[firstNodeIndex][n];
        auto firstEdgeOtherNode = m_edges[localEdgeIndex].first + m_edges[localEdgeIndex].second - firstNodeIndex;
        if (firstEdgeOtherNode == secondNodeIndex)
        {
            edgeIndex = localEdgeIndex;
            break;
        }
    }
    return true;
}

bool GridGeom::Mesh::GetBoundingBox(Point& lowerLeft, Point& upperRight) const
{

    double minx = std::numeric_limits<double>::max();
    double maxx = std::numeric_limits<double>::min();
    double miny = std::numeric_limits<double>::max();
    double maxy = std::numeric_limits<double>::min();
    for (int n = 0; n < GetNumNodes(); n++)
    {
        if (m_nodes[n].IsValid())
        {
            minx = std::min(minx, m_nodes[n].x);
            maxx = std::max(maxx, m_nodes[n].x);
            miny = std::min(miny, m_nodes[n].y);
            maxy = std::max(maxy, m_nodes[n].y);
        }
    }
    lowerLeft = { minx , miny };
    upperRight = { maxx , maxy };

    return true;
}

bool GridGeom::Mesh::OffsetSphericalCoordinates(double minx, double miny)
{
    if(m_projection==Projections::spherical)
    {
        for (int n = 0; n < GetNumNodes(); ++n)
        {
            if(m_nodes[n].x-360.0 >= minx)
            {
                m_nodes[n].x -= 360.0;
            }

            if (m_nodes[n].x < minx)
            {
                m_nodes[n].x += 360.0;
            }
        }
    }

    return true;
}

bool GridGeom::Mesh::BuildRTree()
{
    m_rtree.Clear();
    m_rtree.BuildTree(m_nodes, m_projection);
    return true;
}

bool GridGeom::Mesh::InsertMissingNodesInRTree()
{

    if (m_rtree.Size() < m_nodes.size())
    {
        //insert the missing nodes
        for (int i = m_rtree.Size(); i < m_nodes.size(); ++i)
        {
            m_rtree.InsertNode(m_nodes[i]);
        }
    }
    return true;
}


bool GridGeom::Mesh::GetNodeIndex(Point point, double searchRadius, int& vertexIndex)
{
    if (GetNumNodes() == 0)
    {
        return true;
    }

    if (m_rtree.Empty())
    {
        BuildRTree();
    }

    InsertMissingNodesInRTree();

    auto result = m_rtree.NearestNeighbours(point, searchRadius);
    if (!result.empty())
    {
        vertexIndex = result[0];
    }

    return true;
}