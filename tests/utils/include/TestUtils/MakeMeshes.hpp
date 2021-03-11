//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2021.
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

#pragma once
#include <memory>
#include <string>

#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernelApi/Mesh2D.hpp>

/// @brief This convenience struct collects the unique_ptr for \ref AllocateMesh2dData
struct Mesh2dPointers
{
    std::unique_ptr<int> edge_nodes;
    std::unique_ptr<int> face_nodes;
    std::unique_ptr<int> face_edges;
    std::unique_ptr<int> nodes_per_face;
    std::unique_ptr<double> node_x;
    std::unique_ptr<double> node_y;
    std::unique_ptr<double> edge_x;
    std::unique_ptr<double> edge_y;
    std::unique_ptr<double> face_x;
    std::unique_ptr<double> face_y;
};

/// @brief This function allocates the space reported by \ref meshkernelapi::mkernel_get_mesh2d_dimensions
/// The function should be called between \ref meshkernelapi::mkernel_get_mesh2d_dimensions and \ref meshkernelapi::mkernel_get_mesh2d_data
///
/// @param[in,out] mesh2d The dimensions in the mesh2d are used to allocate memory for the pointers in the mesh2d
/// @returns The Mesh2dPointers instance is only returned so that the allocated memory will be owned by the caller
Mesh2dPointers AllocateMesh2dData(meshkernelapi::Mesh2D mesh2d);

meshkernelapi::Mesh2D ReadLegacyMeshFromFileForApiTesting(std::string filePath);

std::shared_ptr<meshkernel::Mesh2D> ReadLegacyMeshFromFile(std::string filePath, meshkernel::Projection projection = meshkernel::Projection::cartesian);

std::shared_ptr<meshkernel::Mesh2D> MakeRectangularMeshForTesting(int n, int m, double delta, meshkernel::Projection projection, meshkernel::Point origin = {0.0, 0.0});

meshkernelapi::Mesh2D MakeRectangularMeshForApiTesting(int n, int m, double delta);

void DeleteRectangularMeshForApiTesting(const meshkernelapi::Mesh2D& mesh2d);

std::shared_ptr<meshkernel::Mesh2D> MakeSmallSizeTriangularMeshForTestingAsNcFile();

std::shared_ptr<meshkernel::Mesh2D> MakeCurvilinearGridForTesting();
