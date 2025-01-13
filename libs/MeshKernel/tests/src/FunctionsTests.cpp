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

#include <chrono>
#include <gtest/gtest.h>
#include <random>

#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Operations.hpp>

TEST(FunctionsTest, NormalVectorInsideTestCartesian)
{
    // 1 Setup
    meshkernel::Point firstPoint{0.0, 0.0};
    meshkernel::Point secondPoint{1.0, 0.0};
    meshkernel::Point insidePoint{0.0, -1.0};
    meshkernel::Point normal;
    bool flippedNormal;

    // 2 Execute
    NormalVectorInside(firstPoint, secondPoint, insidePoint, normal, flippedNormal, meshkernel::Projection::cartesian);

    // 3 Validation
    ASSERT_EQ(normal.x, 0.0);
    ASSERT_EQ(normal.y, 1.0);
    ASSERT_EQ(flippedNormal, true);
}
