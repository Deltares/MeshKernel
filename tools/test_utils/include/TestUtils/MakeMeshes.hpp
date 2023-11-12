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
#include <filesystem>
#include <memory>
#include <string>

#include <MeshKernel/Definitions.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh1D.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Point.hpp>

std::tuple<size_t,
           size_t,
           std::vector<double>,
           std::vector<double>,
           std::vector<int>,
           std::vector<int>,
           std::vector<int>>
ReadLegacyMeshFile(std::filesystem::path const& file_path);

std::tuple<std::vector<meshkernel::Point>,
           std::vector<meshkernel::Edge>>
ComputeEdgesAndNodes(
    std::filesystem::path const& file_path,
    meshkernel::Mesh::Type meshType);

std::shared_ptr<meshkernel::Mesh2D> MakeRectangularMeshForTesting(
    meshkernel::UInt n,
    meshkernel::UInt m,
    double dim_x,
    double dim_y,
    meshkernel::Projection projection,
    meshkernel::Point const& origin = {0.0, 0.0});

std::shared_ptr<meshkernel::Mesh2D> MakeRectangularMeshForTesting(
    meshkernel::UInt n,
    meshkernel::UInt m,
    double delta,
    meshkernel::Projection projection,
    meshkernel::Point const& origin = {0.0, 0.0});

std::shared_ptr<meshkernel::Mesh2D> ReadLegacyMesh2DFromFile(
    std::filesystem::path const& file_path,
    meshkernel::Projection projection = meshkernel::Projection::cartesian);

std::shared_ptr<meshkernel::Mesh1D> ReadLegacyMesh1DFromFile(
    std::filesystem::path const& file_path,
    meshkernel::Projection projection = meshkernel::Projection::cartesian);

std::tuple<meshkernel::UInt,
           meshkernel::UInt,
           std::vector<double>,
           std::vector<double>,
           std::vector<int>>
MakeRectangularMeshForApiTesting(
    meshkernel::UInt numRows,
    meshkernel::UInt numColumns,
    double delta);

std::shared_ptr<meshkernel::Mesh2D> MakeSmallSizeTriangularMeshForTestingAsNcFile();

std::shared_ptr<meshkernel::Mesh2D> MakeCurvilinearGridForTesting();

std::tuple<std::vector<double>,
           std::vector<double>,
           std::vector<int>,
           std::vector<int>,
           std::vector<int>>
MakeMeshWithFaceNodesForApiTesting();

std::tuple<std::vector<meshkernel::Point>,
           std::vector<meshkernel::Edge>,
           std::vector<std::vector<meshkernel::UInt>>,
           std::vector<meshkernel::UInt>>
MakeMeshWithFaceNodes();
