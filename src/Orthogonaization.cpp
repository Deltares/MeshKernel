#ifndef ORTHOGONALIZATION_CPP
#define ORTHOGONALIZATION_CPP

#define _USE_MATH_DEFINES
#include <vector>
#include <algorithm>
#include "Operations.cpp"
#include "Mesh.hpp"

using namespace GridGeom;

template<typename Point>
class Orthogonalization
{

public:

    std::vector<double> m_aspectRatio;

    bool orthogonalize(Mesh<Point>& mesh)
    {
        //size_t maxNumNeighbours = std::max_element(mesh.m_nodesNumEdges.begin(), mesh.m_nodesNumEdges.end());
        maxNumNeighbours += 1;

        std::vector< std::vector<size_t>> nodesNodes(mesh.m_nodes.size(), std::vector<size_t>(maxNumNeighbours, 0));
        //for each node, determine the neighbours
        for (size_t n = 0; n < mesh.m_nodes.size(); n++)
        {
            for (size_t nn = 0; nn < mesh.m_nodesNumEdges[n]; nn++)
            {
                Edge edge = mesh.m_edges[mesh.m_nodesEdges[n][nn]];
                size_t neighbour = edge.first == n ? edge.second : edge.first;
                //nodesNodes[n] = neighbour;
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

    bool aspectRatio(const Mesh<Point>& mesh)
    {
        std::vector<std::vector<double>> averageEdgesLength(mesh.m_edges.size(), std::vector<double>(2, doubleMissingValue));
        std::vector<double> averageFlowEdgesLength(mesh.m_edges.size(), doubleMissingValue);
        std::vector<bool> curvilinearGridIndicator(mesh.m_nodes.size(), true);
        std::vector<double> edgesLength(mesh.m_edges.size(), 0.0);
        m_aspectRatio.resize(mesh.m_edges.size(), 0.0);

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
                double dinry = Operations<Point>::innerProductTwoSegments(mesh.m_nodes[first], mesh.m_nodes[second], mesh.m_nodes[first], leftCenter);
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
                    m_aspectRatio[edgeIndex] = averageFlowEdgesLength[edgeIndex] / edgeLength;
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
};

#endif
