//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2025.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// contact: delft3d.support@deltares.nl
// Stichting Deltares
// P.O. Box 177
// 2600 MH Delft, The Netherlands
//
// All indications and logos of, and references to, "Delft3D" and "Deltares"
// are registered trademarks of Stichting Deltares, and remain the property of
// Stichting Deltares. All rights reserved.
//
//------------------------------------------------------------------------------

#include "TestMeshGeneration.hpp"

#include "MeshKernelApi/MeshKernel.hpp"

int GenerateCurvilinearMesh(const int meshKernelId,
                            const int nodesX, const int nodesY,
                            const double deltaX, const double deltaY,
                            const double originX, const double originY)
{
    meshkernel::MakeGridParameters makeGridParameters;

    // num_columns and num_rows indicate number of elements in each direction, so value = nodes - 1
    makeGridParameters.num_columns = nodesX - 1;
    makeGridParameters.num_rows = nodesY - 1;
    makeGridParameters.angle = 0.0;
    makeGridParameters.origin_x = originX;
    makeGridParameters.origin_y = originY;
    makeGridParameters.block_size_x = deltaX;
    makeGridParameters.block_size_y = deltaY;

    // Generate curvilinear grid.
    int errorCode = meshkernelapi::mkernel_curvilinear_compute_rectangular_grid(meshKernelId, makeGridParameters);
    return errorCode;
}

int GenerateUnstructuredMesh(int meshKernelId, meshkernel::UInt numRows, meshkernel::UInt numColumns, double delta)
{
    // Set-up new mesh
    auto [num_nodes, num_edges, node_x, node_y, edge_nodes] = MakeRectangularMeshForApiTesting(numRows, numColumns, delta);
    meshkernelapi::Mesh2D mesh2d{};
    mesh2d.num_edges = static_cast<int>(num_edges);
    mesh2d.num_nodes = static_cast<int>(num_nodes);
    mesh2d.node_x = node_x.data();
    mesh2d.node_y = node_y.data();
    mesh2d.edge_nodes = edge_nodes.data();
    int errorCode = meshkernelapi::mkernel_mesh2d_set(meshKernelId, mesh2d);

    return errorCode;
}
