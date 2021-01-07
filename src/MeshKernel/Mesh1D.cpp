#pragma once
#include "MeshKernel/Mesh1D.hpp"

#include <MeshKernel/Entities.hpp>
#include <vector>

meshkernel::Mesh1D::Mesh1D(const std::vector<Edge>& edges,
                           const std::vector<Point>& nodes,
                           Projection projection) : Mesh(edges, nodes, m_projection){};

void meshkernel::Mesh1D::Administrate()
{
    DeleteInvalidNodesAndEdges();
    NodeAdministration();
}