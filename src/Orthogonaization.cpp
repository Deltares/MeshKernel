#ifndef ORTHOGONALIZATION_CPP
#define ORTHOGONALIZATION_CPP

#define _USE_MATH_DEFINES
#include <vector>
#include <algorithm>
#include <numeric>
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

    bool initialize(const Mesh<Point>& mesh)
    {
        size_t m_maxNumNeighbours = *(std::max_element(mesh.m_nodesNumEdges.begin(), mesh.m_nodesNumEdges.end()));
        m_maxNumNeighbours += 1;
        m_nodesNodes.resize(mesh.m_nodes.size(), std::vector<size_t>(m_maxNumNeighbours, 0));
        m_weights.resize(mesh.m_nodes.size(), std::vector<double>(m_maxNumNeighbours, 0.0));
        m_rightHandSide.resize(mesh.m_nodes.size(), std::vector<double>(2, 0.0));
        m_aspectRatios.resize(mesh.m_edges.size(), 0.0);
        //for each node, determine the neighbours
        for (size_t n = 0; n < mesh.m_nodes.size(); n++)
        {
            for (size_t nn = 0; nn < mesh.m_nodesNumEdges[n]; nn++)
            {
                Edge edge = mesh.m_edges[mesh.m_nodesEdges[n][nn]];
                size_t neighbour = edge.first == n ? edge.second : edge.first;
                m_nodesNodes[n][nn] = neighbour;
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
        //allocate small arrays in the stack for speed
        std::vector<size_t> connectedFaces(maximumNumberOfEdgesPerNode,0); //icell
        std::vector<size_t> localNodeAdministration(maximumNumberOfNodesInStencil, 0); //kk2
        std::vector<std::vector<size_t>> positionOfFacesInLocalNodeAdministration(maximumNumberOfNodesPerFace,std::vector<size_t>(connectedFacesNumber, 0));//kkc
        double admXi[maximumNumberOfNodesInStencil]{0};  //xi
        double admEta[maximumNumberOfNodesInStencil]{0}; //eta

        size_t stencilEdge;
        size_t stencilFace;
        std::vector<double> topXi(m_maxNumNeighbours,0.0);
        std::vector<double> topEta(m_maxNumNeighbours,0.0);
        std::vector<std::vector<double>> admNode(mesh.m_nodes.size(), std::vector<double>(m_maxNumNeighbours, 0.0));

        size_t localNodeIndex;
        int numConnectedFaces;
        for (size_t n = 0; n < mesh.m_nodes.size(); n++)
        {
            if (mesh.m_nodesNumEdges[n] < 2)continue;

            orthogonalizationAdministration(mesh, n, connectedFaces, numConnectedFaces, localNodeAdministration, localNodeIndex, positionOfFacesInLocalNodeAdministration);

        }

        return true;
    }

    bool orthogonalizationAdministration(const Mesh<Point>& mesh, int n, std::vector<size_t>& connectedFaces, int & numConnectedFaces, std::vector<size_t>& localNodeAdministration, size_t & localNodeIndex, std::vector<std::vector<size_t>> & positionOfFacesInLocalNodeAdministration)
    {
        int Nsize = connectedFacesNumber;
        int KSize = maximumNumberOfNodesInStencil;

        for (auto& e : connectedFaces)  e = 0;
        

        if (mesh.m_nodesNumEdges[n] < 2) return true;

        int newFaceIndex = -999;
        numConnectedFaces = 0;
        for (int nn = 0; nn < mesh.m_nodesNumEdges[n]; nn++)
        {
            size_t firstEdge = mesh.m_nodesEdges[n][nn];
            int secondIndex = nn + 1;
            if (secondIndex >= mesh.m_nodesNumEdges[n]) secondIndex = 0;

            size_t secondEdge = mesh.m_nodesEdges[n][secondIndex];

            if (mesh.m_edgesNumFaces[firstEdge] < 1 || mesh.m_edgesNumFaces[secondEdge] < 1) continue;

            int firstFaceIndex = std::max(std::min(mesh.m_edgesNumFaces[firstEdge], (size_t)2), (size_t)1);
            int secondFaceIndex = std::max(std::min(mesh.m_edgesNumFaces[secondEdge], (size_t)2), (size_t)1);

            if ((mesh.m_edgesFaces[firstEdge][0] == mesh.m_edgesFaces[secondEdge][0] || mesh.m_edgesFaces[firstEdge][0] == mesh.m_edgesFaces[secondEdge][secondFaceIndex]))
            {
                newFaceIndex = mesh.m_edgesFaces[firstEdge][0];
            }
            else if ((mesh.m_edgesFaces[firstEdge][firstFaceIndex] == mesh.m_edgesFaces[secondEdge][0] || mesh.m_edgesFaces[firstEdge][firstFaceIndex] == mesh.m_edgesFaces[secondEdge][secondFaceIndex]))
            {
                newFaceIndex = mesh.m_edgesFaces[firstEdge][firstFaceIndex];
            }

            if (mesh.m_nodesEdges[n][nn]==2 && nn==2)
            {
                if (connectedFaces[0]) newFaceIndex = -999;
            }
            connectedFaces[numConnectedFaces] = newFaceIndex;
            numConnectedFaces += numConnectedFaces + 1;
        }

        if (numConnectedFaces <= 0) return true;

        for (auto& v : localNodeAdministration) v = 0;
        localNodeIndex = 0;
        localNodeAdministration[localNodeIndex] = n;

        // edge - connected nodes
        for (int nn = 0; nn < mesh.m_nodesNumEdges[n]; nn++)
        {
            size_t edge = mesh.m_nodesEdges[n][nn];
            size_t nodeIndex = mesh.m_edges[edge].first + mesh.m_edges[edge].second - n;
            localNodeIndex++;
            localNodeAdministration[localNodeIndex] = nodeIndex;
        }

        // other - connected nodes
        for (int nn = 0; nn < numConnectedFaces; nn++)
        {
            int faceIndex = connectedFaces[numConnectedFaces];
            if (faceIndex < 0) continue;

            size_t nodeIndex = 0;
            for( int nnn = 0; nnn < mesh.m_facesNodes[faceIndex].size(); nnn++)
            {
                if(mesh.m_facesNodes[faceIndex][nnn] == n)
                {
                    nodeIndex = nnn;
                    break;
                }
            }

            size_t numOtherNodes = mesh.m_facesEdges[faceIndex].size();
            for (int nnn = 0; nnn < numOtherNodes; nnn++)
            {
                nodeIndex += 1;
                if (nodeIndex >= numOtherNodes) nodeIndex -= numOtherNodes;

                int otherFaceIndex = connectedFaces[numConnectedFaces];
                bool isNewNode = true;

                for (int nnnn = 0; nnnn < localNodeIndex; nnnn++)
                {
                    if(otherFaceIndex== localNodeAdministration[nnnn])
                    {
                        isNewNode = false;
                        positionOfFacesInLocalNodeAdministration[nnn][nn] = nnnn;
                    }
                }

                if(isNewNode)
                {
                    localNodeIndex++;
                    localNodeAdministration[localNodeIndex] = otherFaceIndex;
                    positionOfFacesInLocalNodeAdministration[nnn][nn] = localNodeIndex;
                }
            }
        }

        return true;
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
