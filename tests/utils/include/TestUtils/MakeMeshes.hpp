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

#include <MeshKernel/Mesh1D.hpp>
#include <MeshKernel/Mesh2D.hpp>

std::tuple<size_t, size_t, std::shared_ptr<double>, std::shared_ptr<double>, std::vector<int>, std::shared_ptr<int>, std::shared_ptr<int>> ReadLegacyMeshFile(std::string const& filePath);

std::tuple<std::vector<meshkernel::Point>, std::vector<meshkernel::Edge>> ComputeEdgesAndNodes(std::string const& filePath, meshkernel::Mesh::Type meshType);

std::shared_ptr<meshkernel::Mesh2D> MakeRectangularMeshForTesting(int n, int m, double delta, meshkernel::Projection projection, meshkernel::Point origin = {0.0, 0.0});

std::shared_ptr<meshkernel::Mesh2D> ReadLegacyMesh2DFromFile(std::string const& filePath, meshkernel::Projection projection = meshkernel::Projection::cartesian);

std::shared_ptr<meshkernel::Mesh1D> ReadLegacyMesh1DFromFile(std::string const& filePath, meshkernel::Projection projection = meshkernel::Projection::cartesian);

std::tuple<size_t, size_t, std::shared_ptr<double>, std::shared_ptr<double>, std::shared_ptr<int>> MakeRectangularMeshForApiTesting(int n, int m, double delta);

std::shared_ptr<meshkernel::Mesh2D> MakeSmallSizeTriangularMeshForTestingAsNcFile();

std::shared_ptr<meshkernel::Mesh2D> MakeCurvilinearGridForTesting();

std::tuple<std::vector<double>,
           std::vector<double>,
           std::vector<int>,
           std::vector<int>,
           std::vector<int>>
MakeMeshWithFaceNodes();