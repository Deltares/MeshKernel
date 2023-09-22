
#include <gtest/gtest.h>

#include <memory>

#include "TestUtils/MakeMeshes.hpp"

#include <MeshKernel/Contacts.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh1D.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Polygons.hpp>
#include <gmock/gmock-matchers.h>

class ContactsTests : public testing::Test
{
public:
    std::shared_ptr<meshkernel::Mesh1D> MakeMesh1D() const
    {
        // Create 1d mesh
        std::vector<meshkernel::Point> nodes{
            {1.73493900000000, -7.6626510000000},
            {2.35659313023165, 1.67281447902331},
            {5.38347452702839, 10.3513746546384},
            {14.2980910429074, 12.4797224193970},
            {22.9324017677239, 15.3007317677239},
            {25.3723169493137, 24.1623588554512},
            {25.8072280000000, 33.5111870000000}};
        std::vector<meshkernel::Edge> edges{{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}};
        return std::make_shared<meshkernel::Mesh1D>(edges, nodes, meshkernel::Projection::Type::Cartesian);
    }

    std::shared_ptr<meshkernel::Mesh1D> MakeMesh1DOutsideMesh2D() const
    {
        // Create 1d mesh
        std::vector<meshkernel::Point> nodes{
            {-10.0, 5.0},
            {-10.0, -10.0},
            {15.0, 15.0},
            {35.0, 35.0}};
        std::vector<meshkernel::Edge> edges{{0, 1}, {1, 2}, {2, 3}};
        return std::make_shared<meshkernel::Mesh1D>(edges, nodes, meshkernel::Projection::Type::Cartesian);
    }
};

TEST_F(ContactsTests, ComputeSingleContacts1dMeshInside2dMesh)
{

    // Create 1d mesh
    const auto mesh1d = MakeMesh1D();

    // Create 2d mesh
    const auto mesh2d = MakeRectangularMeshForTesting(4, 4, 10, meshkernel::Projection::Type::Cartesian, {0.0, 0.0});

    // Create contacts
    std::vector<bool> onedNodeMask(mesh1d->GetNumNodes(), true);
    meshkernel::Contacts contacts(mesh1d, mesh2d);

    // Set the polygon where to generate the contacts
    std::vector<meshkernel::Point> polygonPoints{{-30, -20}, {40, -20}, {40, 50}, {-40, 50}, {-30, -20}};
    meshkernel::Polygons polygon(polygonPoints, meshkernel::Projection::Type::Cartesian);

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
    const auto mesh2d = MakeRectangularMeshForTesting(4, 4, 10, meshkernel::Projection::Type::Cartesian, {0.0, 0.0});

    // Create contacts
    std::vector<bool> onedNodeMask(mesh1d->GetNumNodes(), true);
    meshkernel::Contacts contacts(mesh1d, mesh2d);

    // Set the polygon where to generate the contacts
    std::vector<meshkernel::Point> polygonPoints{{-30, -20}, {40, -20}, {40, 50}, {-40, 50}, {-30, -20}};
    meshkernel::Polygons polygon(polygonPoints, meshkernel::Projection::Type::Cartesian);

    // Execute
    contacts.ComputeSingleContacts(onedNodeMask, polygon, 5.0);

    // Assert
    ASSERT_THAT(contacts.Mesh1dIndices(), ::testing::ElementsAre(0, 2));
    ASSERT_THAT(contacts.Mesh2dIndices(), ::testing::ElementsAre(0, 4));
}

TEST_F(ContactsTests, ComputeMultipleContacts1dMeshInside2dMesh)
{

    // Create 1d mesh
    const auto mesh1d = MakeMesh1D();

    // Create 2d mesh
    const auto mesh2d = MakeRectangularMeshForTesting(4, 4, 10, meshkernel::Projection::Type::Cartesian, {0.0, 0.0});

    // Create contacts
    std::vector<bool> onedNodeMask(mesh1d->GetNumNodes(), true);
    meshkernel::Contacts contacts(mesh1d, mesh2d);

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
    const auto mesh2d = MakeRectangularMeshForTesting(4, 4, 10, meshkernel::Projection::Type::Cartesian, {0.0, 0.0});

    // Create contacts
    std::vector<bool> onedNodeMask(mesh1d->GetNumNodes(), true);
    meshkernel::Contacts contacts(mesh1d, mesh2d);

    // Create points to connect
    std::vector<meshkernel::Point> pointsToConnect{
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
    using namespace meshkernel::constants;

    // Create 1d mesh
    std::vector<meshkernel::Point> nodes{
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

    std::vector<meshkernel::Edge> edges;
    for (meshkernel::UInt index = 0; index < nodes.size() - 1; ++index)
    {
        edges.emplace_back(meshkernel::Edge{index, index + 1});
    }

    const auto mesh1d = std::make_shared<meshkernel::Mesh1D>(edges, nodes, meshkernel::Projection::Type::Cartesian);

    // Create 2d mesh
    const auto mesh2d = MakeRectangularMeshForTesting(100, 100, 10, meshkernel::Projection::Type::Cartesian, {0.0, 0.0});

    // Create contacts
    std::vector<bool> onedNodeMask(nodes.size(), true);
    meshkernel::Contacts contacts(mesh1d, mesh2d);

    // Set the polygon where to generate the contacts
    std::vector<meshkernel::Point> polygonPoints{
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
    meshkernel::Polygons polygon(polygonPoints, meshkernel::Projection::Type::Cartesian);

    // Execute
    contacts.ComputeContactsWithPolygons(onedNodeMask, polygon);

    // Assert
    ASSERT_THAT(contacts.Mesh1dIndices(), ::testing::ElementsAre(2, 10, 19, 23, 31, 44));
    ASSERT_THAT(contacts.Mesh2dIndices(), ::testing::ElementsAre(1709, 2121, 3999, 4703, 6482, 6805));
}

TEST(Contacts, ComputeBoundaryContacts)
{
    // Create 1d mesh
    std::vector<meshkernel::Point> nodes{
        {-16.1886410000000, 0.89018900000000},
        {-16.1464995876014, 9.78201442138723},
        {-16.1043581752028, 18.6738398427745},
        {-16.0622167628042, 27.5656652641617},
        {-15.7539488236928, 36.1966603330179},
        {-6.86476658679268, 36.4175095626911},
        {2.02441565010741, 36.6383587923643},
        {10.9135970000000, 36.8592080000000}};
    std::vector<meshkernel::Edge> edges{{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}};
    const auto mesh1d = std::make_shared<meshkernel::Mesh1D>(edges, nodes, meshkernel::Projection::Type::Cartesian);

    // Create 2d mesh
    const auto mesh2d = MakeRectangularMeshForTesting(4, 4, 10, meshkernel::Projection::Type::Cartesian, {0.0, 0.0});

    // Create contacts
    std::vector<bool> onedNodeMask(nodes.size(), true);
    meshkernel::Contacts contacts(mesh1d, mesh2d);

    // Set the polygon where to generate the contacts
    std::vector<meshkernel::Point> polygonPoints{{-30, -20}, {40, -20}, {40, 50}, {-40, 50}, {-30, -20}};
    meshkernel::Polygons polygon(polygonPoints, meshkernel::Projection::Type::Cartesian);

    // Execute
    contacts.ComputeBoundaryContacts(onedNodeMask, polygon, 200.0);

    // Assert
    ASSERT_THAT(contacts.Mesh1dIndices(), ::testing::ElementsAre(0, 2, 6, 0, 7, 7, 7, 7));
    ASSERT_THAT(contacts.Mesh2dIndices(), ::testing::ElementsAre(0, 1, 2, 3, 5, 6, 7, 8));
}
