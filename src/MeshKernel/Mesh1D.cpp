#pragma once
#include "MeshKernel/Mesh1D.hpp"

#include <MeshKernel/Entities.hpp>
#include <vector>

meshkernel::Mesh1D::Mesh1D(const std::vector<Edge>& edges,
                           const std::vector<Point>& nodes,
                           const std::vector<bool>& nodeMask,
                           Projection projection) : m_edges(edges), m_nodes(nodes), m_nodeMask(nodeMask), m_projection(projection){};

void meshkernel::Mesh1D::FindFaces()
{
}