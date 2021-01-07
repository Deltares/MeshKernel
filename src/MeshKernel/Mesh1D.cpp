#pragma once
#include "MeshKernel/Mesh1D.hpp"

#include <MeshKernel/Entities.hpp>
#include <vector>

meshkernel::Mesh1D::Mesh1D(const std::vector<Edge>& edges,
                           const std::vector<Point>& nodes,
                           Projection projection) : m_edges(edges), m_nodes(nodes), m_projection(projection){};

void meshkernel::Mesh1D::FindFaces()
{
}