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
#include <memory>
#include <random>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/MeshTransformation.hpp"
#include "MeshKernel/Polygons.hpp"
#include "MeshKernel/UndoActions/UndoAction.hpp"

#include "TestUtils/Definitions.hpp"
#include "TestUtils/MakeMeshes.hpp"

namespace mk = meshkernel;

TEST(MeshTransformationTest, BasicTranslationTest)
{
    // Test basic functionality of the translation class
    double translationX = 1.0;
    double translationY = 2.0;

    mk::Vector vec(translationX, translationY);
    mk::Translation translation(vec);

    EXPECT_EQ(translationX, translation.vector().x());
    EXPECT_EQ(translationY, translation.vector().y());

    translation.identity();
    EXPECT_EQ(0.0, translation.vector().x());
    EXPECT_EQ(0.0, translation.vector().y());

    std::swap(vec.x(), vec.y());
    vec.x() = -vec.x();
    vec.y() = -vec.y();
    translation.reset(vec);

    EXPECT_EQ(-translationY, translation.vector().x());
    EXPECT_EQ(-translationX, translation.vector().y());
}

TEST(MeshTransformationTest, BasicRotationTest)
{
    // Test basic functionality of the rotation class
    double theta = 45.0;

    mk::Rotation rotation(theta);

    EXPECT_EQ(theta, rotation.angle());

    theta *= 1.5;
    rotation.reset(theta);
    EXPECT_EQ(theta, rotation.angle());

    rotation.identity();
    EXPECT_EQ(0.0, rotation.angle());
}

TEST(MeshTransformationTest, PointTranslationTest)
{
    // Test a simple translation of a point in Cartesian coordinates.
    double originX = 12.0;
    double originY = -19.0;
    double translationX = 1.0;
    double translationY = 2.0;

    mk::Vector vec(translationX, translationY);
    mk::Translation translation(vec);
    mk::Point pnt(originX, originY);
    pnt = translation(pnt);

    EXPECT_EQ(originX + translationX, pnt.x);
    EXPECT_EQ(originY + translationY, pnt.y);
}

TEST(MeshTransformationTest, BasicRigidBodyTransformationTest)
{
    // Test the rigid body transformation
    // The test will rotate a series of points about a point that is not the origin
    //
    // Step 1: create composition of all simple transformations
    //    1. translate the rotation point (-1,1) to origin (0,0)
    //    2. rotate by 90 degrees
    //    3. translate back to original rotation point (-1,1)
    //
    // Step 2 apply composite transformation to a series of points
    constexpr double tolerance = 1.0e-10;

    double theta = -90.0;
    mk::Point rotationPoint(-1.0, 1.0);

    mk::RigidBodyTransformation transformation;

    // Expect initial state to be the identity
    EXPECT_EQ(0.0, transformation.rotation().angle());
    EXPECT_EQ(0.0, transformation.translation().vector().x());
    EXPECT_EQ(0.0, transformation.translation().vector().y());

    // First translate
    mk::Vector vec(-rotationPoint.x, -rotationPoint.y);
    mk::Translation translation(vec);
    transformation.compose(translation);

    EXPECT_EQ(0.0, transformation.rotation().angle());
    EXPECT_EQ(-rotationPoint.x, transformation.translation().vector().x());
    EXPECT_EQ(-rotationPoint.y, transformation.translation().vector().y());

    // Second rotate
    mk::Rotation rotation(theta);
    transformation.compose(rotation);

    mk::Vector rotatedVector = rotation(vec);

    EXPECT_EQ(theta, transformation.rotation().angle());
    EXPECT_EQ(rotatedVector.x(), transformation.translation().vector().x());
    EXPECT_EQ(rotatedVector.y(), transformation.translation().vector().y());

    // Third translate back
    vec = {rotationPoint.x, rotationPoint.y};
    translation.reset(vec);
    transformation.compose(translation);

    EXPECT_EQ(theta, transformation.rotation().angle());

    double expectedTranslationX = -2.0;
    double expectedTranslationY = 0.0;

    EXPECT_NEAR(expectedTranslationX, transformation.translation().vector().x(), tolerance);
    EXPECT_NEAR(expectedTranslationY, transformation.translation().vector().y(), tolerance);

    //------------------------------------------------------------
    // Apply the composite transformation to a series of points

    mk::Point pnt1{-3.0, 4.0};
    mk::Point pnt2{-2.0, 2.0};
    mk::Point pnt3{-3.0, 1.0};

    mk::Point expected1{2.0, 3.0};
    mk::Point expected2{0.0, 2.0};
    mk::Point expected3{-1.0, 3.0};

    pnt1 = transformation(pnt1);
    EXPECT_NEAR(expected1.x, pnt1.x, tolerance);
    EXPECT_NEAR(expected1.y, pnt1.y, tolerance);

    pnt2 = transformation(pnt2);
    EXPECT_NEAR(expected2.x, pnt2.x, tolerance);
    EXPECT_NEAR(expected2.y, pnt2.y, tolerance);

    pnt3 = transformation(pnt3);
    EXPECT_NEAR(expected3.x, pnt3.x, tolerance);
    EXPECT_NEAR(expected3.y, pnt3.y, tolerance);
}

TEST(MeshTransformationTest, PointRotationTest)
{
    // Test a simple rotation of a point in Cartesian coordinates.
    double originX = 12.0;
    double originY = -19.0;

    double theta = 45.0;
    double cosTheta = std::cos(theta * M_PI / 180.0);
    double sinTheta = std::sin(theta * M_PI / 180.0);

    mk::Rotation rotation(theta);
    mk::Point pnt(originX, originY);
    pnt = rotation(pnt);

    EXPECT_EQ(originX * cosTheta - originY * sinTheta, pnt.x);
    EXPECT_EQ(originX * sinTheta + originY * cosTheta, pnt.y);
}

TEST(MeshTransformationTest, MeshTranslationTest)
{
    // Test the translation on the entire mesh.

    mk::UInt nx = 11;
    mk::UInt ny = 11;

    double delta = 10.0;

    // Make two copies of the same mesh so that they can be compared in the test.
    std::shared_ptr<mk::Mesh2D> originalMesh = MakeRectangularMeshForTesting(nx, ny, delta, mk::Projection::cartesian);
    std::shared_ptr<mk::Mesh2D> mesh = MakeRectangularMeshForTesting(nx, ny, delta, mk::Projection::cartesian);

    mk::Vector vec(1.0, 2.0);

    mk::Translation translation(vec);

    const std::vector<mk::Point> meshPoints(mesh->Nodes());

    std::unique_ptr<mk::UndoAction> undoAction = mk::MeshTransformation::Compute(*mesh, translation);

    for (mk::UInt i = 0; i < mesh->GetNumNodes(); ++i)
    {
        EXPECT_EQ(originalMesh->Node(i).x + vec.x(), mesh->Node(i).x);
        EXPECT_EQ(originalMesh->Node(i).y + vec.y(), mesh->Node(i).y);
    }

    undoAction->Restore();

    ASSERT_EQ(mesh->GetNumNodes(), static_cast<mk::UInt>(meshPoints.size()));

    constexpr double tolerance = 1.0e-10;

    for (mk::UInt i = 0; i < meshPoints.size(); ++i)
    {
        EXPECT_NEAR(meshPoints[i].x, mesh->Node(i).x, tolerance);
        EXPECT_NEAR(meshPoints[i].y, mesh->Node(i).y, tolerance);
    }
}

TEST(MeshTransformationTest, MeshRotationTest)
{
    // Test the rotation on the entire mesh.

    mk::UInt nx = 11;
    mk::UInt ny = 11;

    double delta = 10.0;

    // Make two copies of the same mesh so that they can be compared in the test.
    std::shared_ptr<mk::Mesh2D> originalMesh = MakeRectangularMeshForTesting(nx, ny, delta, mk::Projection::cartesian);
    std::shared_ptr<mk::Mesh2D> mesh = MakeRectangularMeshForTesting(nx, ny, delta, mk::Projection::cartesian);

    double theta = 45.0;

    mk::Rotation rotation(theta);
    double cosTheta = std::cos(theta * M_PI / 180.0);
    double sinTheta = std::sin(theta * M_PI / 180.0);

    const std::vector<mk::Point> meshPoints(mesh->Nodes());

    std::unique_ptr<mk::UndoAction> undoAction = mk::MeshTransformation::Compute(*mesh, rotation);

    for (mk::UInt i = 0; i < mesh->GetNumNodes(); ++i)
    {
        mk::Point expected{originalMesh->Node(i).x * cosTheta - originalMesh->Node(i).y * sinTheta,
                           originalMesh->Node(i).x * sinTheta + originalMesh->Node(i).y * cosTheta};
        EXPECT_EQ(expected.x, mesh->Node(i).x);
        EXPECT_EQ(expected.y, mesh->Node(i).y);
    }

    undoAction->Restore();

    ASSERT_EQ(mesh->GetNumNodes(), static_cast<mk::UInt>(meshPoints.size()));

    constexpr double tolerance = 1.0e-10;

    for (mk::UInt i = 0; i < meshPoints.size(); ++i)
    {
        EXPECT_NEAR(meshPoints[i].x, mesh->Node(i).x, tolerance);
        EXPECT_NEAR(meshPoints[i].y, mesh->Node(i).y, tolerance);
    }
}

TEST(MeshTransformationTest, IncorrectProjectionTest)
{
    // Test correct failure with non Cartesian projection.

    // Generate mesh in spherical coordinate system
    std::shared_ptr<mk::Mesh2D> mesh = MakeRectangularMeshForTesting(11, 11, 10.0, mk::Projection::spherical);

    mk::Rotation rotation(45.0);

    // Should throw an exception with spherical coordinate system
    EXPECT_THROW([[maybe_unused]] auto action = mk::MeshTransformation::Compute(*mesh, rotation), mk::MeshKernelError);

    // Change projection to Projection::sphericalAccurate
    mesh->m_projection = mk::Projection::sphericalAccurate;

    // Should throw an exception with spherical-accurate coordinate system
    EXPECT_THROW([[maybe_unused]] auto action = mk::MeshTransformation::Compute(*mesh, rotation), mk::MeshKernelError);
}
