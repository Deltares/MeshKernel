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

#include <gtest/gtest.h>

#include <memory>

#include "TestUtils/MakeMeshes.hpp"

#include <MeshKernel/Contacts.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh1D.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/Utilities/Utilities.hpp>
#include <gmock/gmock-matchers.h>

using namespace meshkernel;

class ContactsTests : public testing::Test
{
public:
    std::unique_ptr<Mesh1D> MakeMesh1D() const
    {
        // Create 1d mesh
        std::vector<Point> nodes{
            {1.73493900000000, -7.6626510000000},
            {2.35659313023165, 1.67281447902331},
            {5.38347452702839, 10.3513746546384},
            {14.2980910429074, 12.4797224193970},
            {22.9324017677239, 15.3007317677239},
            {25.3723169493137, 24.1623588554512},
            {25.8072280000000, 33.5111870000000}};
        std::vector<Edge> edges{{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}};
        return std::make_unique<Mesh1D>(edges, nodes, Projection::cartesian);
    }

    std::unique_ptr<Mesh1D> MakeMesh1DOutsideMesh2D() const
    {
        // Create 1d mesh
        std::vector<Point> nodes{
            {-10.0, 5.0},
            {-10.0, -10.0},
            {15.0, 15.0},
            {35.0, 35.0}};
        std::vector<Edge> edges{{0, 1}, {1, 2}, {2, 3}};
        return std::make_unique<Mesh1D>(edges, nodes, Projection::cartesian);
    }
};

TEST_F(ContactsTests, ComputeSingleContacts1dMeshInside2dMesh)
{

    // Create 1d mesh
    const auto mesh1d = MakeMesh1D();

    // Create 2d mesh
    const auto mesh2d = MakeRectangularMeshForTesting(4, 4, 10, Projection::cartesian, {0.0, 0.0});

    // Create contacts
    std::vector<bool> onedNodeMask(mesh1d->GetNumNodes(), true);
    Contacts contacts(*mesh1d, *mesh2d);

    // Set the polygon where to generate the contacts
    std::vector<Point> polygonPoints{{-30, -20}, {40, -20}, {40, 50}, {-40, 50}, {-30, -20}};
    Polygons polygon(polygonPoints, Projection::cartesian);

    // Execute
    contacts.ComputeSingleContacts(onedNodeMask, polygon, 0.0);

    // Assert
    ASSERT_THAT(contacts.Mesh1dIndices(), ::testing::ElementsAre(1, 2, 3, 4, 5));
    ASSERT_THAT(contacts.Mesh2dIndices(), ::testing::ElementsAre(0, 1, 4, 7, 8));
}

TEST_F(ContactsTests, ComputeSingleContacts_onAMesh2DWithLaterals_ShouldComputeContats)
{

    // Create 1d mesh
    const auto mesh1d = MakeMesh1DOutsideMesh2D();

    // Create 2d mesh
    const auto mesh2d = MakeRectangularMeshForTesting(4, 4, 10, Projection::cartesian, {0.0, 0.0});

    // Create contacts
    std::vector<bool> onedNodeMask(mesh1d->GetNumNodes(), true);
    Contacts contacts(*mesh1d, *mesh2d);

    // Set the polygon where to generate the contacts
    std::vector<Point> polygonPoints{{-30, -20}, {40, -20}, {40, 50}, {-40, 50}, {-30, -20}};
    Polygons polygon(polygonPoints, Projection::cartesian);

    // Execute
    contacts.ComputeSingleContacts(onedNodeMask, polygon, 5.0);

    // Assert
    ASSERT_THAT(contacts.Mesh1dIndices(), ::testing::ElementsAre(0, 2));
    ASSERT_THAT(contacts.Mesh2dIndices(), ::testing::ElementsAre(0, 4));
}

TEST_F(ContactsTests, ComputeSingleContactsGivenFaceNodesofMesh2D)
{
    // Create 1d mesh
    const auto mesh1d = MakeMesh1D();

    // create 2d mesh
    auto const [nodes, edges, faceNodes, numFaceNodes] = MakeMeshWithFaceNodes();
    auto const mesh2d = std::make_shared<Mesh2D>(edges,
                                                 nodes,
                                                 faceNodes,
                                                 numFaceNodes,
                                                 Projection::cartesian);
    // Create 1D node mask
    std::vector<bool> onedNodeMask(mesh1d->GetNumNodes(), true);

    // Create contacts
    Contacts contacts(*mesh1d, *mesh2d);

    // Set the polygon where to generate the contacts
    std::vector<Point> const polygonPoints{{-30, -20}, {40, -20}, {40, 50}, {-40, 50}, {-30, -20}};
    Polygons const polygon(polygonPoints, Projection::cartesian);

    // Execute
    ASSERT_NO_THROW(contacts.ComputeSingleContacts(onedNodeMask, polygon, 0.0));

    auto const& mesh1dIndices = contacts.Mesh1dIndices();
    auto const& mesh2dIndices = contacts.Mesh2dIndices();

    // Assert
    ASSERT_EQ(5, mesh1dIndices.size());

    ASSERT_EQ(1, mesh1dIndices[0]);
    ASSERT_EQ(2, mesh1dIndices[1]);
    ASSERT_EQ(3, mesh1dIndices[2]);
    ASSERT_EQ(4, mesh1dIndices[3]);
    ASSERT_EQ(5, mesh1dIndices[4]);

    ASSERT_EQ(0, mesh2dIndices[0]);
    ASSERT_EQ(3, mesh2dIndices[1]);
    ASSERT_EQ(4, mesh2dIndices[2]);
    ASSERT_EQ(5, mesh2dIndices[3]);
    ASSERT_EQ(8, mesh2dIndices[4]);
}

TEST_F(ContactsTests, ComputeMultipleContacts1dMeshInside2dMesh)
{

    // Create 1d mesh
    const auto mesh1d = MakeMesh1D();

    // Create 2d mesh
    const auto mesh2d = MakeRectangularMeshForTesting(4, 4, 10, Projection::cartesian, {0.0, 0.0});

    // Create contacts
    std::vector<bool> onedNodeMask(mesh1d->GetNumNodes(), true);
    Contacts contacts(*mesh1d, *mesh2d);

    // Execute
    contacts.ComputeMultipleContacts(onedNodeMask);

    // Assert
    ASSERT_THAT(contacts.Mesh1dIndices(), ::testing::ElementsAre(1, 2, 3, 4, 5));
    ASSERT_THAT(contacts.Mesh2dIndices(), ::testing::ElementsAre(0, 1, 4, 7, 8));
}

TEST_F(ContactsTests, ComputeContactsWithPoints)
{
    // Create 1d mesh
    const auto mesh1d = MakeMesh1D();

    // Create 2d mesh
    const auto mesh2d = MakeRectangularMeshForTesting(4, 4, 10, Projection::cartesian, {0.0, 0.0});

    // Create contacts
    std::vector<bool> onedNodeMask(mesh1d->GetNumNodes(), true);
    Contacts contacts(*mesh1d, *mesh2d);

    // Create points to connect
    std::vector<Point> pointsToConnect{
        {2.9225159, 15.9190083},
        {16.3976765, 5.9373722},
        {27.6269741, 17.7656116},
        {22.3367290, 21.4588184}};

    // Execute
    contacts.ComputeContactsWithPoints(onedNodeMask, pointsToConnect);

    // Assert
    ASSERT_THAT(contacts.Mesh1dIndices(), ::testing::ElementsAre(2, 3, 4, 5));
    ASSERT_THAT(contacts.Mesh2dIndices(), ::testing::ElementsAre(1, 3, 7, 8));
}

TEST(Contacts, ComputeContactsWithPolygons)
{
    using namespace constants;

    // Create 1d mesh
    std::vector<Point> nodes{
        {144.578313, 230.636833},
        {148.479997, 250.145254},
        {152.326276, 269.66441},
        {155.992653, 289.21842},
        {163.108569, 307.723781},
        {171.485584, 325.759464},
        {180.71473, 343.361252},
        {191.75036, 359.914696},
        {203.615349, 375.861699},
        {216.043509, 391.3969},
        {232.539143, 402.399733},
        {249.32359, 413.080745},
        {266.401484, 423.278547},
        {283.618936, 433.246546},
        {302.83515, 437.07713},
        {322.63432, 438.89845},
        {342.529084, 438.89845},
        {362.423848, 438.89845},
        {382.318612, 438.89845},
        {402.213377, 438.89845},
        {422.108141, 438.89845},
        {442.002905, 438.89845},
        {461.654996, 441.830623},
        {481.316269, 444.846117},
        {501.057403, 447.313758},
        {520.392165, 451.956905},
        {539.189727, 458.275803},
        {557.4759, 466.112735},
        {574.862173, 475.610926},
        {591.415617, 486.646556},
        {607.969062, 497.682185},
        {624.522506, 508.717815},
        {641.145408, 519.647697},
        {657.848458, 530.455553},
        {670.513477, 545.357135},
        {681.549107, 561.910579},
        {692.85851, 578.277567},
        {704.219703, 594.609283},
        {713.16624, 612.094922},
        {718.736774, 631.193896},
        {724.244832, 650.310657},
        {729.569544, 669.479619},
        {734.250619, 688.809373},
        {738.557122, 708.232391},
        {742.655652, 727.700409},
        {746.521516, 747.213434},
        {749.978934, 766.805471},
        {752.151462, 786.57487}};

    std::vector<Edge> edges;
    for (UInt index = 0; index < nodes.size() - 1; ++index)
    {
        edges.emplace_back(Edge{index, index + 1});
    }

    const auto mesh1d = std::make_shared<Mesh1D>(edges, nodes, Projection::cartesian);

    // Create 2d mesh
    const auto mesh2d = MakeRectangularMeshForTesting(100, 100, 10, Projection::cartesian, {0.0, 0.0});

    // Create contacts
    std::vector<bool> onedNodeMask(nodes.size(), true);
    Contacts contacts(*mesh1d, *mesh2d);

    // Set the polygon where to generate the contacts
    std::vector<Point> polygonPoints{
        {260.004578, 180.002838},
        {170.307068, 180.002838},
        {169.094971, 272.124969},
        {283.035034, 284.246307},
        {300.004883, 205.457642},
        {260.004578, 180.002838},
        {missing::doubleValue, missing::doubleValue},
        {212.731567, 422.429504},
        {153.337280, 424.853821},
        {149.700867, 520.612366},
        {245.459045, 525.460876},
        {268.489502, 474.551270},
        {212.731567, 422.429504},
        {missing::doubleValue, missing::doubleValue},
        {483.036316, 307.276855},
        {378.793213, 329.095245},
        {404.247925, 403.035400},
        {510.915283, 394.550476},
        {483.036316, 307.276855},
        {missing::doubleValue, missing::doubleValue},
        {476.975647, 498.793945},
        {434.551147, 495.157532},
        {416.369202, 587.279663},
        {526.672974, 589.703979},
        {476.975647, 498.793945},
        {missing::doubleValue, missing::doubleValue},
        {719.401428, 370.307800},
        {655.158630, 395.762604},
        {652.734314, 475.763428},
        {753.341003, 487.884766},
        {719.401428, 370.307800},
        {missing::doubleValue, missing::doubleValue},
        {632.128113, 667.280518},
        {575.158142, 732.735718},
        {583.643005, 778.796753},
        {645.461609, 792.130249},
        {691.522522, 732.735718},
        {632.128113, 667.280518},
        {missing::doubleValue, missing::doubleValue},
    };
    Polygons polygon(polygonPoints, Projection::cartesian);

    // Execute
    contacts.ComputeContactsWithPolygons(onedNodeMask, polygon);

    // Assert
    ASSERT_THAT(contacts.Mesh1dIndices(), ::testing::ElementsAre(2, 10, 19, 23, 31, 44));
    ASSERT_THAT(contacts.Mesh2dIndices(), ::testing::ElementsAre(1709, 2121, 3999, 4703, 6482, 6805));
}

TEST(Contacts, ComputeBoundaryContacts)
{
    // Create 1d mesh
    std::vector<Point> nodes{
        {-16.1886410000000, 0.89018900000000},
        {-16.1464995876014, 9.78201442138723},
        {-16.1043581752028, 18.6738398427745},
        {-16.0622167628042, 27.5656652641617},
        {-15.7539488236928, 36.1966603330179},
        {-6.86476658679268, 36.4175095626911},
        {2.02441565010741, 36.6383587923643},
        {10.9135970000000, 36.8592080000000}};
    std::vector<Edge> edges{{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}};
    const auto mesh1d = std::make_shared<Mesh1D>(edges, nodes, Projection::cartesian);

    // Create 2d mesh
    const auto mesh2d = MakeRectangularMeshForTesting(4, 4, 10, Projection::cartesian, {0.0, 0.0});

    // Create contacts
    std::vector<bool> onedNodeMask(nodes.size(), true);
    Contacts contacts(*mesh1d, *mesh2d);

    // Set the polygon where to generate the contacts
    std::vector<Point> polygonPoints{{-30, -20}, {40, -20}, {40, 50}, {-40, 50}, {-30, -20}};
    Polygons polygon(polygonPoints, Projection::cartesian);

    // Execute
    contacts.ComputeBoundaryContacts(onedNodeMask, polygon, 200.0);

    // Assert
    ASSERT_THAT(contacts.Mesh1dIndices(), ::testing::ElementsAre(0, 2, 6, 0, 7, 7, 7, 7));
    ASSERT_THAT(contacts.Mesh2dIndices(), ::testing::ElementsAre(0, 1, 2, 3, 5, 6, 7, 8));
}

TEST_F(ContactsTests, SetIndices)
{
    // Create 1d mesh
    const auto mesh1d = MakeMesh1D();

    // Create 2d mesh
    const auto mesh2d = MakeRectangularMeshForTesting(4, 4, 10, Projection::cartesian, {0.0, 0.0});

    Contacts contacts(*mesh1d, *mesh2d);
    std::vector<UInt> mesh1dIndices;
    std::vector<UInt> mesh2dIndices;

    // Contacts should accept empty mesh1dIndices and mesh2dIndices arrays

    mesh1dIndices = {1, 2, 3};
    EXPECT_THROW(contacts.SetIndices(mesh1dIndices, mesh2dIndices), AlgorithmError);

    mesh2dIndices = {1, 2, 3, 4};
    EXPECT_THROW(contacts.SetIndices(mesh1dIndices, mesh2dIndices), AlgorithmError);

    mesh2dIndices = {4, 5, 6};
    EXPECT_NO_THROW(contacts.SetIndices(mesh1dIndices, mesh2dIndices));

    EXPECT_EQ(mesh1dIndices, contacts.Mesh1dIndices());
    EXPECT_EQ(mesh2dIndices, contacts.Mesh2dIndices());
}

TEST_F(ContactsTests, AllFalseNodeMask)
{

    // Create 1d mesh

    std::vector<Point> nodes{
        {0.5, 0.5},
        {1.5, 1.5},
        {2.5, 2.5},
        {3.5, 3.5},
        {4.5, 4.5}};
    std::vector<Edge> edges{{0, 1}, {1, 2}, {2, 3}, {3, 4}};

    auto mesh1d = std::make_unique<Mesh1D>(edges, nodes, Projection::cartesian);

    mesh1d->Administrate();

    // Create 2d mesh
    auto mesh2d = MakeRectangularMeshForTesting(6, 6, 1.0, Projection::cartesian, {0.0, 0.0}, true, true);
    mesh2d->Administrate();

    // Create contacts
    std::vector<bool> onedNodeMask(mesh1d->GetNumNodes(), false);
    Contacts contacts(*mesh1d, *mesh2d);

    // Set the polygon where to generate the contacts
    std::vector<Point> points{{0.5, 2.5}, {3.5, 1.5}, {4.5, 2.5}};

    // Execute
    contacts.ComputeContactsWithPoints(onedNodeMask, points);

    auto m1dIndices = contacts.Mesh1dIndices();
    auto m2dIndices = contacts.Mesh2dIndices();

    ASSERT_EQ(m1dIndices.size(), 0);
    ASSERT_EQ(m2dIndices.size(), 0);
}

TEST_F(ContactsTests, AllTrueNodeMask)
{

    // Create 1d mesh

    std::vector<Point> nodes{
        {0.5, 0.5},
        {1.5, 1.5},
        {2.5, 2.5},
        {3.5, 3.5},
        {4.5, 4.5}};
    std::vector<Edge> edges{{0, 1}, {1, 2}, {2, 3}, {3, 4}};

    auto mesh1d = std::make_unique<Mesh1D>(edges, nodes, Projection::cartesian);

    mesh1d->Administrate();

    // Create 2d mesh
    auto mesh2d = MakeRectangularMeshForTesting(6, 6, 1.0, Projection::cartesian, {0.0, 0.0}, true, true);
    mesh2d->Administrate();

    // Create contacts
    std::vector<bool> onedNodeMask(mesh1d->GetNumNodes(), true);
    Contacts contacts(*mesh1d, *mesh2d);

    // Set the polygon where to generate the contacts
    std::vector<Point> points{{0.5, 2.5}, {3.5, 1.5}, {4.5, 2.5}};

    // Execute
    contacts.ComputeContactsWithPoints(onedNodeMask, points);

    auto m1dIndices = contacts.Mesh1dIndices();
    auto m2dIndices = contacts.Mesh2dIndices();

    ASSERT_EQ(m1dIndices.size(), 3);
    ASSERT_EQ(m2dIndices.size(), 3);

    EXPECT_EQ(m1dIndices[0], 1);
    EXPECT_EQ(m1dIndices[1], 2);
    EXPECT_EQ(m1dIndices[2], 3);

    EXPECT_EQ(m2dIndices[0], 2);
    EXPECT_EQ(m2dIndices[1], 16);
    EXPECT_EQ(m2dIndices[2], 22);
}

TEST_F(ContactsTests, AllMixedBooleanNodeMask)
{

    // Create 1d mesh

    std::vector<Point> nodes{
        {0.5, 0.5},
        {1.5, 1.5},
        {2.5, 2.5},
        {3.5, 3.5},
        {4.5, 4.5},
        {5.5, 5.5},
        {6.5, 6.5},
        {8.5, 8.5},
        {9.5, 9.5}};
    std::vector<Edge> edges{{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {7, 8}};

    auto mesh1d = std::make_unique<Mesh1D>(edges, nodes, Projection::cartesian);

    mesh1d->Administrate();

    // Create 2d mesh
    auto mesh2d = MakeRectangularMeshForTesting(11, 11, 1.0, Projection::cartesian, {0.0, 0.0}, true, true);
    mesh2d->Administrate();

    // Create contacts
    std::vector<bool> onedNodeMask(mesh1d->GetNumNodes(), true);
    onedNodeMask[1] = false;
    onedNodeMask[2] = false;
    Contacts contacts(*mesh1d, *mesh2d);

    // Set the polygon where to generate the contacts
    std::vector<Point> points{{0.5, 2.5}, {3.5, 1.5}, {4.5, 2.5}, {7.5, 5.5}};

    // Execute
    contacts.ComputeContactsWithPoints(onedNodeMask, points);

    auto m1dIndices = contacts.Mesh1dIndices();
    auto m2dIndices = contacts.Mesh2dIndices();

    ASSERT_EQ(m1dIndices.size(), 2);
    ASSERT_EQ(m2dIndices.size(), 2);

    EXPECT_EQ(m1dIndices[0], 3);
    EXPECT_EQ(m1dIndices[1], 6);

    EXPECT_EQ(m2dIndices[0], 42);
    EXPECT_EQ(m2dIndices[1], 75);
}

TEST(Contacts, GenerateContactsMorePointsThanOneDMesh)
{

    // This test checks that the contacts generated are correct when the
    // number of points is greater than the number of nodes in the oned mesh

    // Create 1d mesh
    std::vector<Point> nodes{{-9.0, -3.0},
                             {-8.18181818, -2.36363636},
                             {-7.36363636, -1.72727273},
                             {-6.54545455, -1.09090909},
                             {-5.72727273, -0.45454545},
                             {-4.90909091, 0.18181818},
                             {-4.09090909, 0.81818182},
                             {-3.27272727, 1.45454545},
                             {-2.45454545, 2.09090909},
                             {-1.63636364, 2.72727273},
                             {-0.81818182, 3.36363636},
                             {0.0, 4.0},
                             {1.0, 4.11111111},
                             {2.0, 4.22222222},
                             {3.0, 4.33333333},
                             {4.0, 4.44444444},
                             {5.0, 4.55555556},
                             {6.0, 4.66666667},
                             {7.0, 4.77777778},
                             {8.0, 4.88888889},
                             {9.0, 5.0}};
    std::vector<Edge> edges{nodes.size() - 1};

    for (meshkernel::UInt i = 0; i < edges.size(); ++i)
    {
        edges[i].first = i;
        edges[i].second = i + 1;
    }

    const auto mesh1d = std::make_shared<Mesh1D>(edges, nodes, Projection::cartesian);

    // Create 2d mesh
    const auto mesh2d = MakeRectangularMeshForTesting(32, 8, 0.5, Projection::cartesian, {-8.0, 2.5});

    // Create contacts
    std::vector<bool> onedNodeMask(nodes.size(), true);
    Contacts contacts(*mesh1d, *mesh2d);

    std::vector<Point> points{{-1.75, 2.75},
                              {-1.25, 2.75},
                              {-1.25, 3.25},
                              {-0.75, 3.25},
                              {-0.75, 3.75},
                              {-0.25, 3.75},
                              {0.25, 3.75},
                              {-0.25, 4.25},
                              {0.25, 4.25},
                              {0.75, 4.25},
                              {1.25, 4.25},
                              {1.75, 4.25},
                              {2.25, 4.25},
                              {2.75, 4.25},
                              {3.25, 4.25},
                              {3.75, 4.25},
                              {4.25, 4.25},
                              {4.25, 4.75},
                              {4.75, 4.75},
                              {5.25, 4.75},
                              {5.75, 4.75},
                              {6.25, 4.75},
                              {6.75, 4.75},
                              {7.25, 4.75}};

    // Execute
    contacts.ComputeContactsWithPoints(onedNodeMask, points);

    std::vector<meshkernel::UInt> expected1dIndices{9, 9, 10, 10, 10, 11, 11, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 15, 16, 16, 17, 17, 18, 18};
    std::vector<meshkernel::UInt> expected2dIndices{84, 91, 92, 99, 100, 107, 114, 108, 115, 122, 129, 136, 143, 150, 157, 164, 171, 172, 179, 186, 193, 200, 207, 214};

    auto m1dIndices = contacts.Mesh1dIndices();
    auto m2dIndices = contacts.Mesh2dIndices();

    ASSERT_EQ(m1dIndices.size(), expected1dIndices.size());
    ASSERT_EQ(m2dIndices.size(), expected2dIndices.size());

    for (size_t i = 0; i < m1dIndices.size(); ++i)
    {
        EXPECT_EQ(expected1dIndices[i], m1dIndices[i]);
    }

    for (size_t i = 0; i < m2dIndices.size(); ++i)
    {
        EXPECT_EQ(expected2dIndices[i], m2dIndices[i]);
    }
}
