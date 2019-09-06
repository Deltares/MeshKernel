#ifndef ORTHOGONALIZATION_CPP
#define ORTHOGONALIZATION_CPP

#define _USE_MATH_DEFINES
#include <vector>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include "Operations.cpp"
#include "Mesh.hpp"

using namespace GridGeom;

template<typename Point>
class Orthogonalization
{
public:

    size_t m_maxNumNeighbours;
    std::vector< std::vector<size_t>> m_nodesNodes;
    std::vector<std::vector<double>>  m_weights;
    std::vector<std::vector<double>>  m_rightHandSide;
    std::vector<double> m_aspectRatios;
    std::vector<int> m_nodesTypes;                             //types of nodes,  1=internal, 2=on ring, 3=corner point, 0/-1=other (e.g. 1d)
    std::vector<int> m_faceNumNodes;                           //number of face nodes

    bool initialize(const Mesh<Point>& mesh)
    {
        m_maxNumNeighbours = *(std::max_element(mesh.m_nodesNumEdges.begin(), mesh.m_nodesNumEdges.end()));
        m_maxNumNeighbours += 1;
        m_nodesNodes.resize(mesh.m_nodes.size(), std::vector<size_t>(m_maxNumNeighbours, 0));
        m_weights.resize(mesh.m_nodes.size(), std::vector<double>(m_maxNumNeighbours, 0.0));
        m_rightHandSide.resize(mesh.m_nodes.size(), std::vector<double>(2, 0.0));
        m_aspectRatios.resize(mesh.m_edges.size(), 0.0);
        
        //for each node, determine the neighbouring nodes
        for (size_t n = 0; n < mesh.m_nodes.size(); n++)
        {
            for (size_t nn = 0; nn < mesh.m_nodesNumEdges[n]; nn++)
            {
                Edge edge = mesh.m_edges[mesh.m_nodesEdges[n][nn]];
                size_t neighbour = edge.first == n ? edge.second : edge.first;
                m_nodesNodes[n][nn] = neighbour;
            }
        }

        classifyNodes(mesh);

        return true;
    }

    //MAKENETNODESCODING
    bool classifyNodes(const Mesh<Point>& mesh)
    {
        m_nodesTypes.resize(mesh.m_nodes.size(), 0);

        for (int e = 0; e < mesh.m_edges.size(); e++)
        {
            size_t first = mesh.m_edges[e].first;
            size_t second = mesh.m_edges[e].second;

            if (mesh.m_edgesNumFaces[e] == 0)
            {
                m_nodesTypes[first] = -1;
                m_nodesTypes[second] = -1;
            }
            else if (mesh.m_edgesNumFaces[e] == 1)
            {
                m_nodesTypes[first] += 1;
                m_nodesTypes[second] += 1;
            }
        }

        for (int n = 0; n < mesh.m_nodes.size(); n++)
        {
            if (m_nodesTypes[n] == 1 || m_nodesTypes[n] == 2)
            {
                if (mesh.m_nodesNumEdges[n] == 2)
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
            
            if (mesh.m_nodesNumEdges[n] < 2)
            {
                //hanging node
                m_nodesTypes[n] = -1;
            }
        }
        return true;
    }


    //TODO: Unit test me please
    bool computeWeights(const Mesh<Point>& mesh)
    {
        double localOrthogonalizationToSmoothingFactor = 1.0;
        double localOrthogonalizationToSmoothingFactorSymmetric = 1.0 - localOrthogonalizationToSmoothingFactor;
        double mu = 1.0;
        for (size_t n = 0; n < mesh.m_nodes.size(); n++)
        {
            //the check for inside a polygon is now skipped
            for (size_t nn = 0; nn < mesh.m_nodesNumEdges[n]; nn++)
            {
                size_t edgeIndex = mesh.m_nodesEdges[n][nn];
                double aspectRatio = m_aspectRatios[edgeIndex];

                if (aspectRatio != doubleMissingValue)
                {
                    m_weights[n][nn] = localOrthogonalizationToSmoothingFactor * aspectRatio + localOrthogonalizationToSmoothingFactorSymmetric * mu;

                    if (mesh.m_edgesNumFaces[edgeIndex] == 1)
                    {
                        //boundary nodes
                        Point neighbouringNode = mesh.m_nodes[m_nodesNodes[n][nn]];
                        double neighbouringNodeDistance = Operations<Point>::distance(neighbouringNode, mesh.m_nodes[n]);
                        double aspectRatioByNodeDistance = aspectRatio * neighbouringNodeDistance;

                        size_t leftFace = mesh.m_edgesFaces[edgeIndex][0];
                        Point massCenterLeftFace = mesh.m_facesMasscenters[leftFace];
                        Point normal;
                        bool flippedNormal;
                        Operations<Point>::normalVectorInside(mesh.m_nodes[n], neighbouringNode, massCenterLeftFace, normal, flippedNormal);

                        m_rightHandSide[n][0] += localOrthogonalizationToSmoothingFactor * neighbouringNodeDistance * normal.x / 2.0 +
                            localOrthogonalizationToSmoothingFactorSymmetric * aspectRatioByNodeDistance * normal.x * 0.5 / mu;

                        m_rightHandSide[n][1] += localOrthogonalizationToSmoothingFactor * neighbouringNodeDistance * normal.y / 2.0 +
                            localOrthogonalizationToSmoothingFactorSymmetric * aspectRatioByNodeDistance * normal.y * 0.5 / mu;

                        m_weights[n][nn] = localOrthogonalizationToSmoothingFactor * 0.5 * aspectRatio +
                            localOrthogonalizationToSmoothingFactorSymmetric * 0.5 * mu;
                    }

                }
                else
                {
                    m_weights[n][nn] = 0.0;
                }

                // normalize
                double factor = std::accumulate(m_weights[n].begin(), m_weights[n].end(), 0.0);
                if (factor > 1e-14)
                {
                    factor = 1.0 / factor;
                    for (auto& w : m_weights[n]) w = w * factor;
                    m_rightHandSide[n][0] = factor * m_rightHandSide[n][0];
                    m_rightHandSide[n][1] = factor * m_rightHandSide[n][1];
                }
            }
        }
        return true;
    }

    bool solveWeights(const Mesh<Point>& mesh)
    {
        computeOperators(mesh);

        return true;
    }

    bool computeOperators(const Mesh<Point>& mesh)
    {
        //allocate small administration arrays only once
        std::vector<int> sharedFaces(maximumNumberOfEdgesPerNode,-1); //icell
        std::vector<size_t> connectedNodes(maximumNumberOfConnectedNodes, 0); //kk2
        std::vector<std::vector<size_t>> faceNodeMapping(maximumNumberOfConnectedNodes, std::vector<size_t>(maximumNumberOfNodesPerFace, 0));//kkc
        std::vector<double> xi(maximumNumberOfConnectedNodes,0.0);
        std::vector<double> eta(maximumNumberOfConnectedNodes,0.0);

        std::vector<std::vector<double>> nodeConnectedFaceNodes(mesh.m_nodes.size(), std::vector<double>(m_maxNumNeighbours, 0.0));
        std::vector<size_t> numNodeConnectedFaceNodes(mesh.m_nodes.size(), 0.0);
        
        computeFacesNumEdges(mesh);

        for (size_t n = 0; n < mesh.m_nodes.size(); n++)
        {
            int numSharedFaces = 0;
            int numConnectedNodes = 0;
            orthogonalizationAdministration(mesh, n, sharedFaces, numSharedFaces, connectedNodes, numConnectedNodes, faceNodeMapping);
            computeXiEta(mesh, n, sharedFaces, numSharedFaces, connectedNodes, numConnectedNodes, faceNodeMapping, xi,eta);
        }

        return true;
    }

    bool computeFacesNumEdges(const Mesh<Point>& mesh)
    {
        // Cache the result for later calls.
        m_faceNumNodes.resize(mesh.m_numFaces, false);
        for(int f=0; f< mesh.m_numFaces; f++)
        {
            m_faceNumNodes[f] = mesh.m_facesNodes.size();
        }
        return true;
    }

    // computes the shared faces and the connected nodes of a stencil node and the faceNodeMapping in the connectedNodes array for each shared face.
    bool orthogonalizationAdministration(const Mesh<Point>& mesh, const int currentNode, std::vector<int>& sharedFaces, int& numSharedFaces, std::vector<size_t>& connectedNodes, int& numConnectedNodes, std::vector<std::vector<size_t>>& faceNodeMapping)
    {
        for (auto& f : sharedFaces)  f = -1;

        if (mesh.m_nodesNumEdges[currentNode] < 2) return true;

        // 1. For the currentNode, find the shared faces
        int newFaceIndex = -999;
        numSharedFaces = 0;
        for (int e = 0; e < mesh.m_nodesNumEdges[currentNode]; e++)
        {
            size_t firstEdge = mesh.m_nodesEdges[currentNode][e]; 
            
            int secondEdgeIndex = e + 1;
            
            if (secondEdgeIndex >= mesh.m_nodesNumEdges[currentNode]) 
                secondEdgeIndex = 0;

            size_t secondEdge = mesh.m_nodesEdges[currentNode][secondEdgeIndex];

            if (mesh.m_edgesNumFaces[firstEdge] < 1 || mesh.m_edgesNumFaces[secondEdge] < 1) continue;

            // find the face that the first and the second edge share
            int firstFaceIndex = std::max(std::min(mesh.m_edgesNumFaces[firstEdge], size_t(2)), size_t(1)) - 1;
            int secondFaceIndex = std::max(std::min(mesh.m_edgesNumFaces[secondEdge], size_t(2)), size_t(1)) - 1;

            if (mesh.m_edgesFaces[firstEdge][0] != newFaceIndex &&
               (mesh.m_edgesFaces[firstEdge][0] == mesh.m_edgesFaces[secondEdge][0] || mesh.m_edgesFaces[firstEdge][0] == mesh.m_edgesFaces[secondEdge][secondFaceIndex]))
            {
                newFaceIndex = mesh.m_edgesFaces[firstEdge][0];
            }
            else if (mesh.m_edgesFaces[firstEdge][firstFaceIndex] != newFaceIndex &&
                (mesh.m_edgesFaces[firstEdge][firstFaceIndex] == mesh.m_edgesFaces[secondEdge][0] || mesh.m_edgesFaces[firstEdge][firstFaceIndex] == mesh.m_edgesFaces[secondEdge][secondFaceIndex]))
            {
                newFaceIndex = mesh.m_edgesFaces[firstEdge][firstFaceIndex];
            }
            else
            {
                newFaceIndex = -999;
            }

            //corner face (already found in the first iteration)
            if (mesh.m_nodesNumEdges[currentNode] == 2 && e == 1 && m_nodesTypes[currentNode]==3)
            {
                if (sharedFaces[0] == newFaceIndex) newFaceIndex = -999;
            }
            sharedFaces[numSharedFaces] = newFaceIndex;
            numSharedFaces += 1;
        }

        if (numSharedFaces < 1) return true;

        for (auto& v : connectedNodes) v = 0;
        int connectedNodesIndex = 0;
        connectedNodes[connectedNodesIndex] = currentNode;

        // edge connected nodes
        for (int e = 0; e < mesh.m_nodesNumEdges[currentNode]; e++)
        {
            size_t edgeIndex = mesh.m_nodesEdges[currentNode][e];
            size_t nodeIndex = mesh.m_edges[edgeIndex].first + mesh.m_edges[edgeIndex].second - currentNode;
            connectedNodesIndex++;
            connectedNodes[connectedNodesIndex] = nodeIndex;
        }

        // for each face store the positions of the its nodes in the connectedNodes (compressed array)
        if (faceNodeMapping.size() < numSharedFaces) 
            faceNodeMapping.resize(numSharedFaces);

        for (int f = 0; f < numSharedFaces; f++)
        {
            int faceIndex = sharedFaces[f];
            if (faceIndex < 0) continue;

            // find the stencil node position  in the current face
            size_t faceNodeIndex = 0;
            size_t numFaceNodes = m_faceNumNodes[faceIndex];
            for (int n = 0; n < numFaceNodes; n++)
            {
                if (mesh.m_facesNodes[faceIndex][n] == currentNode)
                {
                    faceNodeIndex = n;
                    break;
                }
            }

            for (int n = 0; n < numFaceNodes; n++)
            {

                if (faceNodeIndex >= numFaceNodes) 
                    faceNodeIndex -= numFaceNodes;
                
                int node = mesh.m_facesNodes[faceIndex][faceNodeIndex];


                bool isNewNode = true;
                for (int n = 0; n < connectedNodesIndex + 1; n++)
                {
                    if (node == connectedNodes[n])
                    {
                        isNewNode = false;
                        faceNodeMapping[f][faceNodeIndex] = n;
                        break;
                    }
                }

                if (isNewNode)
                {
                    connectedNodesIndex++;
                    connectedNodes[connectedNodesIndex] = node;
                    faceNodeMapping[f][faceNodeIndex] = connectedNodesIndex;
                }

                //update node index
                faceNodeIndex += 1;
            }
        }

        // compute the number of connected nodes
        numConnectedNodes = connectedNodesIndex + 1;

        return true;
    }


    bool computeXiEta(const Mesh<Point>& mesh,
        const int currentNode,
        const std::vector<int>& sharedFaces,
        const int& numSharedFaces,
        const std::vector<size_t>& connectedNodes,
        const size_t& numConnectedNodes,
        const std::vector<std::vector<size_t>>& faceNodeMapping,
        std::vector<double>& xi,
        std::vector<double>& eta)
    {
        std::fill(xi.begin(), xi.end(), 0.0);
        std::fill(eta.begin(), eta.end(), 0.0);
        // the angles for the squared nodes connected to the stencil nodes, first the ones directly connected, then the others
        std::vector<double> thetaSquare(numConnectedNodes, doubleMissingValue);
        // for each shared face, a bollean indicating if it is squared or not
        std::vector<bool> isSquareFace(numSharedFaces, false);

        int numNonStencilQuad = 0;
        //loop over the connected edges
        for (int f = 0; f < numSharedFaces; f++)
        {
            size_t edgeIndex = mesh.m_nodesEdges[currentNode][f];
            size_t nextNode = connectedNodes[f + 1]; // the first entry is always the stencil node 
            int faceLeft = mesh.m_edgesFaces[edgeIndex][0];
            int faceRigth = faceLeft;

            if (mesh.m_edgesNumFaces[edgeIndex] == 2)
                faceRigth = mesh.m_edgesFaces[edgeIndex][1];

            //check if it is a rectangular node (not currentNode itself) 
            bool isSquare = true;
            for (int e = 0; e < mesh.m_nodesNumEdges[nextNode]; e++)
            {
                size_t edge = mesh.m_nodesEdges[nextNode][e];
                for (int ff = 0; ff < mesh.m_edgesNumFaces[edge]; ff++)
                {
                    size_t face = mesh.m_edgesFaces[edge][ff]; 
                    if (face != faceLeft && face != faceRigth)
                    {
                        isSquare = isSquare && mesh.m_facesNodes[face].size() == 4;
                    }
                }
                if (!isSquare)
                {
                    break;
                }
            }

            //Compute optimal angle thetaSquare based on the node type
            if (isSquare)
            {
                if (m_nodesTypes[nextNode] == 1 || m_nodesTypes[nextNode] == 4)
                {
                    // Inner node
                    numNonStencilQuad = mesh.m_nodesNumEdges[nextNode] - 2;
                    thetaSquare[f + 1] = (2.0 - double(numNonStencilQuad) * 0.5) * M_PI;
                }
                if (m_nodesTypes[nextNode] == 2)
                {
                    // Inner node
                    numNonStencilQuad = mesh.m_nodesNumEdges[nextNode] - 1 - mesh.m_edgesNumFaces[edgeIndex];
                    thetaSquare[f + 1] = (1.0 - double(numNonStencilQuad) * 0.5) * M_PI;
                }
                if (m_nodesTypes[nextNode] == 3)
                {
                    // Inner node
                    numNonStencilQuad = mesh.m_nodesNumEdges[nextNode] - 2.0; 
                    thetaSquare[f + 1] = 0.5 * M_PI;
                }
            }

            int leftFaceIndex = f - 1; if (leftFaceIndex < 0) leftFaceIndex = leftFaceIndex + numSharedFaces;
            if (sharedFaces[f] > 1)
            {
                if (m_faceNumNodes[sharedFaces[f]]==4 ) numNonStencilQuad += 1;
            }
            if (sharedFaces[leftFaceIndex] > 1)
            {
                if (m_faceNumNodes[sharedFaces[leftFaceIndex]] == 4 ) numNonStencilQuad += 1;
            }
            if (numNonStencilQuad > 3)
            {
                isSquare = false;
            }

            isSquareFace[f] = isSquareFace[f] || isSquare;
            isSquareFace[leftFaceIndex] = isSquareFace[leftFaceIndex] || isSquare;
        }

        for (int f = 0; f < numSharedFaces; f++)
        {
            // boundary face
            if (sharedFaces[f] < 0) continue;

            // non boundary face 
            if (m_faceNumNodes[sharedFaces[f]] == 4)
            {
                for (int n = 0; n < m_faceNumNodes[sharedFaces[f]]; n++)
                {
                    if (faceNodeMapping[f][n] <= numSharedFaces) continue;
                    thetaSquare[faceNodeMapping[f][n]] = 0.5 * M_PI;
                }
            }
        }

        // Compute internal angle
        int numSquaredTriangles = 0.0;
        int numTriangles = 0;
        double phiSquaredTriangles = 0.0;
        double phiQuads = 0.0;
        double phiTriangles = 0.0;
        double phiTot = 0.0;
        numNonStencilQuad = 0;
        for (int f = 0; f < numSharedFaces; f++)
        {
            // boundary face
            if (sharedFaces[f] < 0) continue;

            int numFaceNodes = m_faceNumNodes[sharedFaces[f]];
            double phi = optimalEdgeAngle(numFaceNodes);

            if (isSquareFace[f] || numFaceNodes == 4)
            {
                int nextNode = f + 2; if (nextNode > numSharedFaces) nextNode = nextNode - numSharedFaces;
                bool isBoundaryEdge = mesh.m_edgesNumFaces[mesh.m_nodesEdges[currentNode][f]] == 1;
                phi = optimalEdgeAngle(numFaceNodes, thetaSquare[f + 1], thetaSquare[nextNode], isBoundaryEdge);
                if (numFaceNodes == 3)
                {
                    numSquaredTriangles += 1;
                    phiSquaredTriangles += phi;
                }
                else if (numFaceNodes == 4)
                {
                    numNonStencilQuad += 1;
                    phiQuads += phi;
                }
            }
            else
            {
                numTriangles += 1;
                phiTriangles += phi;
            }
            phiTot += phi;
        }


        double factor = 1.0;
        if (m_nodesTypes[currentNode] == 2) factor = 0.5;
        if (m_nodesTypes[currentNode] == 3) factor = 0.25;
        double mu = 0.0;
        double muSquaredTriangles = 0.0;
        double muTriangles = 0.0;
        double minPhi = 15.0 / 180.0 * M_PI;
        if (numTriangles > 0)
        {
            muTriangles = (factor * 2.0 * M_PI - (phiTot - phiTriangles)) / phiTriangles;
            muTriangles = std::max(muTriangles, double(numTriangles) * minPhi / phiTriangles);
        }
        else if (numSquaredTriangles > 0)
        {
            muSquaredTriangles = std::max(factor * 2.0 * M_PI - (phiTot - phiSquaredTriangles), double(numSquaredTriangles) * minPhi) / phiSquaredTriangles;
        }

        if (phiTot > 1e-18)
        {
            mu = factor * 2.0 * M_PI / (phiTot - (1.0 - muTriangles) * phiTriangles - (1.0 - muSquaredTriangles) * phiSquaredTriangles);
        }
        else if (numSharedFaces > 0)
        {
            //TODO: error
            return false;
        }

        double phi0 = 0.0;
        double dPhi0 = 0.0;
        double dPhi = 0.0;
        double dTheta = 0.0;
        for (int f = 0; f < numSharedFaces; f++)
        {
            phi0 = phi0 + 0.5 * dPhi;
            if (sharedFaces[f] < 0)
            {
                if (m_nodesTypes[currentNode] == 2)
                {
                    dPhi = M_PI;
                }
                else if (m_nodesTypes[currentNode] == 3)
                {
                    dPhi = 1.5 * M_PI;
                }
                else
                {
                    //TODO: error
                    return false;
                }
                phi0 = phi0 + 0.5 * dPhi;
                continue;
            }

            int numFaceNodes = m_faceNumNodes[sharedFaces[f]];
            if (numFaceNodes > maximumNumberOfEdgesPerNode)
            {
                //TODO: error
                return false;
            }

            dPhi0 = optimalEdgeAngle(numFaceNodes);
            if (isSquareFace[f])
            {
                int nextNode = f + 2; if (nextNode > numSharedFaces) nextNode = nextNode - numSharedFaces;
                bool isBoundaryEdge = mesh.m_edgesNumFaces[mesh.m_nodesEdges[currentNode][f]] == 1;
                dPhi0 = optimalEdgeAngle(numFaceNodes, thetaSquare[f + 1], thetaSquare[nextNode], isBoundaryEdge);
                if (numFaceNodes == 3)
                {
                    dPhi0 = muSquaredTriangles * dPhi0;
                }
            }
            else if (numFaceNodes == 3)
            {
                dPhi0 = muTriangles * dPhi0; 
            }

            dPhi = mu * dPhi0;
            phi0 = phi0 + 0.5 * dPhi;

            // determine the index of the current stencil node
            int nodeIndex = 0;
            for (int n = 0; n < numFaceNodes; n++)
            {
                if (mesh.m_facesNodes[sharedFaces[f]][n] == currentNode)
                {
                    nodeIndex = n;
                    break;
                }
            }

            // optimal angle
            dTheta = 2.0 * M_PI / double(numFaceNodes);

            // orientation of the face (necessary for folded cells)
            int previousNode = nodeIndex + 1; if (previousNode > numFaceNodes) previousNode -= numFaceNodes;
            int nextNode = nodeIndex - 1; if (nextNode < 0) nextNode += numFaceNodes;
            if (faceNodeMapping[f][nextNode] - faceNodeMapping[f][previousNode] == -1 ||
                faceNodeMapping[f][nextNode] - faceNodeMapping[f][previousNode] == mesh.m_nodesNumEdges[currentNode]-1)
            {
                dTheta = -dTheta;
            }

            double aspectRatio = (1.0 - std::cos(dTheta)) / std::sin(std::abs(dTheta)) * std::tan(0.5 * dPhi);
            double radius = std::cos(0.5 * dPhi) / (1.0-cos(dTheta));

            for (int n = 0; n < numFaceNodes; n++)
            {
                double theta = dTheta * (n - nodeIndex);
                double xip = radius - radius * std::cos(theta);
                double ethap = -radius * std::sin(theta);

                xi[faceNodeMapping[f][n]] = xip * std::cos(phi0) - aspectRatio * ethap * std::sin(phi0);
                eta[faceNodeMapping[f][n]] = xip * std::sin(phi0) + aspectRatio * ethap * std::cos(phi0);
            }
        }

        return true;
    }

    double optimalEdgeAngle(const int numFaceNodes, const double theta1 = doubleMissingValue, const double theta2 = doubleMissingValue, bool isBoundaryEdge = false)
    {
        double angle = M_PI * (1 - 2.0 / double(numFaceNodes));

        if (theta1 != doubleMissingValue && theta2 != doubleMissingValue && numFaceNodes == 3)
        {
            angle = 0.25 * M_PI;
            if (theta1 + theta2 == M_PI && !isBoundaryEdge)
            {
                angle = 0.5 * M_PI;
            }
        }
        return angle;
    }

    bool aspectRatio(const Mesh<Point>& mesh)
    {
        std::vector<std::vector<double>> averageEdgesLength(mesh.m_edges.size(), std::vector<double>(2, doubleMissingValue));
        std::vector<double> averageFlowEdgesLength(mesh.m_edges.size(), doubleMissingValue);
        std::vector<bool> curvilinearGridIndicator(mesh.m_nodes.size(), true);
        std::vector<double> edgesLength(mesh.m_edges.size(), 0.0);

        for (size_t e = 0; e < mesh.m_edges.size(); e++)
        {
            size_t first = mesh.m_edges[e].first;
            size_t second = mesh.m_edges[e].second;

            if (first == second) continue; 
            double edgeLength = Operations<Point>::distance(mesh.m_nodes[first], mesh.m_nodes[second]);
            edgesLength[e] = edgeLength;

            Point leftCenter;
            Point rightCenter;
            if (mesh.m_edgesNumFaces[e] > 0)
            {
                leftCenter = mesh.m_facesCircumcenters[mesh.m_edgesFaces[e][0]];
            }
            else
            {
                leftCenter = mesh.m_nodes[first];
            }

            //find right cell center, if it exists
            if (mesh.m_edgesNumFaces[e] == 2)
            {
                rightCenter = mesh.m_facesCircumcenters[mesh.m_edgesFaces[e][1]];
            }
            else
            {
                //otherwise, make ghost node by imposing boundary condition
                double dinry = Operations<Point>::outerProductTwoSegments(mesh.m_nodes[first], mesh.m_nodes[second], mesh.m_nodes[first], leftCenter);
                dinry = dinry / std::max(edgeLength * edgeLength, minimumEdgeLength);

                double x0_bc = (1.0 - dinry) * mesh.m_nodes[first].x + dinry * mesh.m_nodes[second].x;
                double y0_bc = (1.0 - dinry) * mesh.m_nodes[first].y + dinry * mesh.m_nodes[second].y;
                rightCenter.x = 2.0 * x0_bc - leftCenter.x;
                rightCenter.y = 2.0 * y0_bc - leftCenter.y;
            }

            averageFlowEdgesLength[e] = Operations<Point>::distance(leftCenter, rightCenter);
        }

        // Compute normal length
        for (int f = 0; f < mesh.m_facesNodes.size(); f++)
        {
            size_t numberOfFaceNodes = mesh.m_facesNodes[f].size();
            if (numberOfFaceNodes < 3) continue;

            for (int n = 0; n < numberOfFaceNodes; n++)
            {
                if (numberOfFaceNodes != 4) curvilinearGridIndicator[mesh.m_facesNodes[f][n]] = false;
                size_t edgeIndex = mesh.m_facesEdges[f][n];

                if (mesh.m_edgesNumFaces[edgeIndex] < 1) continue;

                //get the other links in the right numbering
                //TODO: ask why only 3 are requested, why an average lenght stored in averageEdgesLength is needed?
                //int kkm1 = n - 1; if (kkm1 < 0) kkm1 = kkm1 + numberOfFaceNodes;
                //int kkp1 = n + 1; if (kkp1 >= numberOfFaceNodes) kkp1 = kkp1 - numberOfFaceNodes;
                //
                //size_t klinkm1 = mesh.m_facesEdges[f][kkm1];
                //size_t klinkp1 = mesh.m_facesEdges[f][kkp1];
                //

                double edgeLength = edgesLength[edgeIndex];
                if (edgeLength != 0.0)
                {
                    m_aspectRatios[edgeIndex] = averageFlowEdgesLength[edgeIndex] / edgeLength;
                }

                //quads
                if (numberOfFaceNodes == 4)
                {
                    int kkp2 = n + 2; if (kkp2 >= numberOfFaceNodes) kkp2 = kkp2 - numberOfFaceNodes;
                    size_t klinkp2 = mesh.m_facesEdges[f][kkp2];
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

        for (size_t e = 0; e < mesh.m_edges.size(); e++)
        {
            size_t first = mesh.m_edges[e].first;
            size_t second = mesh.m_edges[e].second;

            if (first < 0 || second < 0) continue;
            if (mesh.m_edgesNumFaces[e] < 1) continue;
            // Consider only quads
            if (!curvilinearGridIndicator[first] || !curvilinearGridIndicator[second]) continue;

            if (mesh.m_edgesNumFaces[e] == 1)
            {
                if (averageEdgesLength[e][0] != 0.0 && averageEdgesLength[e][0] != doubleMissingValue)
                {
                    m_aspectRatios[e] = averageFlowEdgesLength[e] / averageEdgesLength[e][0];
                }
            }
            else
            {
                if (averageEdgesLength[e][0] != 0.0 && averageEdgesLength[e][1] != 0.0 &&
                    averageEdgesLength[e][0] != doubleMissingValue && averageEdgesLength[e][1] != doubleMissingValue)
                {
                    m_aspectRatios[e] = curvilinearToOrthogonalRatio * m_aspectRatios[e] +
                        (1.0 - curvilinearToOrthogonalRatio) * averageFlowEdgesLength[e] / (0.5 * (averageEdgesLength[e][0] + averageEdgesLength[e][1]));
                }
            }
        }
        return true;
    }

};

#endif
