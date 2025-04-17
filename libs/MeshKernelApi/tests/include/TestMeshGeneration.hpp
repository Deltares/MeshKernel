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

#pragma once

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <MeshKernel/Definitions.hpp>
#include <MeshKernel/Parameters.hpp>

#include <MeshKernelApi/MeshKernel.hpp>

#include <TestUtils/MakeMeshes.hpp>

/// @brief Generate an unstructured mesh
int GenerateUnstructuredMesh(int meshKernelId, meshkernel::UInt numRows = 2, meshkernel::UInt numColumns = 3, double delta = 1.0);

/// @brief Generate an structured curvilinear mesh
int GenerateCurvilinearMesh(const int meshKernelId,
                            const int nodesX, const int nodesY,
                            const double deltaX, const double deltaY,
                            const double originX, const double originY);
