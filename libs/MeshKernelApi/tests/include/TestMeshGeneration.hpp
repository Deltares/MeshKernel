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
