#pragma once

#include <vector>
#include <algorithm>
#include <numeric>
#include <string>

#include "Operations.cpp"
#include "Constants.cpp"
#include "Mesh.hpp"
#include "Entities.hpp"
#include "Orthogonalizer.hpp"

GridGeom::Orthogonalizer::Orthogonalizer()
{
}

bool GridGeom::Orthogonalizer::Compute(Mesh& mesh)
{
    // Compute mesh aspect ratios
    AspectRatio(mesh);

    mesh.ComputeNodeNeighbours();
    m_weights.resize(mesh.GetNumNodes(), std::vector<double>(mesh.m_maxNumNeighbours, 0.0));
    m_rhs.resize(mesh.GetNumNodes(), std::vector<double>(2, 0.0));
    std::fill(m_rhs.begin(), m_rhs.end(), std::vector<double>(2, 0.0));
    m_aspectRatios.resize(mesh.GetNumEdges(), 0.0);

    for (auto n = 0; n < mesh.GetNumNodes(); n++)
    {
        if (mesh.m_nodesTypes[n] != 1 && mesh.m_nodesTypes[n] != 2)
        {
            continue;
        }

        for (auto nn = 0; nn < mesh.m_nodesNumEdges[n]; nn++)
        {

            const auto edgeIndex = mesh.m_nodesEdges[n][nn];
            double aspectRatio = m_aspectRatios[edgeIndex];
            m_weights[n][nn] = 0.0;

            if (aspectRatio != doubleMissingValue)
            {
                // internal nodes
                m_weights[n][nn] = aspectRatio;

                if (mesh.m_edgesNumFaces[edgeIndex] == 1)
                {
                    // boundary nodes
                    m_weights[n][nn] = 0.5 * aspectRatio;

                    // compute the edge length
                    Point neighbouringNode = mesh.m_nodes[mesh.m_nodesNodes[n][nn]];
                    double neighbouringNodeDistance = Distance(neighbouringNode, mesh.m_nodes[n], mesh.m_projection);
                    double aspectRatioByNodeDistance = aspectRatio * neighbouringNodeDistance;

                    auto leftFace = mesh.m_edgesFaces[edgeIndex][0];
                    bool flippedNormal;
                    Point normal;
                    NormalVectorInside(mesh.m_nodes[n], neighbouringNode, mesh.m_facesMassCenters[leftFace], normal, flippedNormal, mesh.m_projection);

                    if (mesh.m_projection == Projections::spherical && mesh.m_projection != Projections::sphericalAccurate)
                    {
                        normal.x = normal.x * std::cos(degrad_hp * 0.5 * (mesh.m_nodes[n].y + neighbouringNode.y));
                    }

                    m_rhs[n][0] += neighbouringNodeDistance * normal.x * 0.5;
                    m_rhs[n][1] += neighbouringNodeDistance * normal.y * 0.5;
                }

            }
        }

        // normalize
        double factor = std::accumulate(m_weights[n].begin(), m_weights[n].end(), 0.0);
        if (std::abs(factor) > 1e-14)
        {
            factor = 1.0 / factor;
            for (auto& w : m_weights[n]) w = w * factor;
            m_rhs[n][0] = factor * m_rhs[n][0];
            m_rhs[n][1] = factor * m_rhs[n][1];
        }

    }
    return true;
}

bool GridGeom::Orthogonalizer::AspectRatio(const Mesh& mesh)
{
    std::vector<std::vector<double>> averageEdgesLength(mesh.GetNumEdges(), std::vector<double>(2, doubleMissingValue));
    std::vector<double> averageFlowEdgesLength(mesh.GetNumEdges(), doubleMissingValue);
    std::vector<bool> curvilinearGridIndicator(mesh.GetNumNodes(), true);
    std::vector<double> edgesLength(mesh.GetNumEdges(), 0.0);
    m_aspectRatios.resize(mesh.GetNumEdges(), 0.0);

    for (auto e = 0; e < mesh.GetNumEdges(); e++)
    {
        auto first = mesh.m_edges[e].first;
        auto second = mesh.m_edges[e].second;

        if (first == second) continue;
        double edgeLength = Distance(mesh.m_nodes[first], mesh.m_nodes[second], mesh.m_projection);
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
            double dinry = InnerProductTwoSegments(mesh.m_nodes[first], mesh.m_nodes[second], mesh.m_nodes[first], leftCenter, mesh.m_projection);
            dinry = dinry / std::max(edgeLength * edgeLength, minimumEdgeLength);

            double x0_bc = (1.0 - dinry) * mesh.m_nodes[first].x + dinry * mesh.m_nodes[second].x;
            double y0_bc = (1.0 - dinry) * mesh.m_nodes[first].y + dinry * mesh.m_nodes[second].y;
            rightCenter.x = 2.0 * x0_bc - leftCenter.x;
            rightCenter.y = 2.0 * y0_bc - leftCenter.y;
        }

        averageFlowEdgesLength[e] = Distance(leftCenter, rightCenter, mesh.m_projection);
    }

    // Compute normal length
    for (int f = 0; f < mesh.GetNumFaces(); f++)
    {
        auto numberOfFaceNodes = mesh.GetNumFaceEdges(f);
        if (numberOfFaceNodes < 3) continue;

        for (int n = 0; n < numberOfFaceNodes; n++)
        {
            if (numberOfFaceNodes != 4) curvilinearGridIndicator[mesh.m_facesNodes[f][n]] = false;
            auto edgeIndex = mesh.m_facesEdges[f][n];

            if (mesh.m_edgesNumFaces[edgeIndex] < 1) continue;

            //get the other links in the right numbering
            //TODO: ask why only 3 are requested, why an average length stored in averageEdgesLength is needed?
            //int kkm1 = n - 1; if (kkm1 < 0) kkm1 = kkm1 + numberOfFaceNodes;
            //int kkp1 = n + 1; if (kkp1 >= numberOfFaceNodes) kkp1 = kkp1 - numberOfFaceNodes;
            //
            //std::size_t klinkm1 = mesh.m_facesEdges[f][kkm1];
            //std::size_t klinkp1 = mesh.m_facesEdges[f][kkp1];
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
                auto klinkp2 = mesh.m_facesEdges[f][kkp2];
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

    for (auto e = 0; e < mesh.GetNumEdges(); e++)
    {
        auto first = mesh.m_edges[e].first;
        auto second = mesh.m_edges[e].second;

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