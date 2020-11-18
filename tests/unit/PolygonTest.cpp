#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Polygons.hpp>
#include <TestUtils/MakeMeshes.hpp>
#include <gtest/gtest.h>

TEST(Polygons, MeshBoundaryToPolygon)
{
    auto mesh = ReadLegacyMeshFromFile("../../../../tests/data/SmallTriangularGrid_net.nc");

    const std::vector<meshkernel::Point> polygon;
    auto polygons = std::make_shared<meshkernel::Polygons>(polygon, meshkernel::Projections::cartesian);

    int numNodesBoundaryPolygons;
    std::vector<meshkernel::Point> meshBoundaryPolygon;
    polygons->MeshBoundaryToPolygon(*mesh, meshBoundaryPolygon, numNodesBoundaryPolygons);

    ASSERT_EQ(8, numNodesBoundaryPolygons);

    constexpr double tolerance = 1e-5;

    ASSERT_NEAR(227.00204467773401, meshBoundaryPolygon[0].x, tolerance);
    ASSERT_NEAR(259.25222778320301, meshBoundaryPolygon[1].x, tolerance);
    ASSERT_NEAR(428.00329589843801, meshBoundaryPolygon[2].x, tolerance);
    ASSERT_NEAR(536.00396728515602, meshBoundaryPolygon[3].x, tolerance);
    ASSERT_NEAR(503.75378417968801, meshBoundaryPolygon[4].x, tolerance);
    ASSERT_NEAR(350.75280761718801, meshBoundaryPolygon[5].x, tolerance);
    ASSERT_NEAR(322.25262451171898, meshBoundaryPolygon[6].x, tolerance);
    ASSERT_NEAR(227.00204467773401, meshBoundaryPolygon[7].x, tolerance);

    ASSERT_NEAR(360.37924194335898, meshBoundaryPolygon[0].y, tolerance);
    ASSERT_NEAR(241.87805175781301, meshBoundaryPolygon[1].y, tolerance);
    ASSERT_NEAR(210.37774658203099, meshBoundaryPolygon[2].y, tolerance);
    ASSERT_NEAR(310.87875366210898, meshBoundaryPolygon[3].y, tolerance);
    ASSERT_NEAR(432.37997436523398, meshBoundaryPolygon[4].y, tolerance);
    ASSERT_NEAR(458.63024902343801, meshBoundaryPolygon[5].y, tolerance);
    ASSERT_NEAR(454.88018798828102, meshBoundaryPolygon[6].y, tolerance);
    ASSERT_NEAR(360.37924194335898, meshBoundaryPolygon[7].y, tolerance);
}

TEST(Polygons, CreatePointsInPolygons)
{
    // Prepare

    std::vector<meshkernel::Point> nodes;

    nodes.push_back({302.002502, 472.130371});
    nodes.push_back({144.501526, 253.128174});
    nodes.push_back({368.752930, 112.876755});
    nodes.push_back({707.755005, 358.879242});
    nodes.push_back({301.252502, 471.380371});
    nodes.push_back({302.002502, 472.130371});

    meshkernel::Polygons polygons(nodes, meshkernel::Projections::cartesian);

    // Execute
    std::vector<std::vector<meshkernel::Point>> generatedPoints;
    polygons.CreatePointsInPolygons(generatedPoints);

    // Assert
    const double tolerance = 1e-5;

    ASSERT_NEAR(302.00250199999999, generatedPoints[0][0].x, tolerance);
    ASSERT_NEAR(472.13037100000003, generatedPoints[0][0].y, tolerance);

    ASSERT_NEAR(144.50152600000001, generatedPoints[0][1].x, tolerance);
    ASSERT_NEAR(253.12817400000000, generatedPoints[0][1].y, tolerance);

    ASSERT_NEAR(368.75292999999999, generatedPoints[0][2].x, tolerance);
    ASSERT_NEAR(112.87675500000000, generatedPoints[0][2].y, tolerance);

    ASSERT_NEAR(707.75500499999998, generatedPoints[0][3].x, tolerance);
    ASSERT_NEAR(358.87924199999998, generatedPoints[0][3].y, tolerance);

    ASSERT_NEAR(301.25250199999999, generatedPoints[0][4].x, tolerance);
    ASSERT_NEAR(471.38037100000003, generatedPoints[0][4].y, tolerance);
}

TEST(Polygons, RefinePolygon)
{
    // Prepare
    std::vector<meshkernel::Point> nodes;

    nodes.push_back({0, 0});
    nodes.push_back({3, 0});
    nodes.push_back({3, 3});
    nodes.push_back({0, 3});
    nodes.push_back({0, 0});

    meshkernel::Polygons polygons(nodes, meshkernel::Projections::cartesian);

    // Execute
    std::vector<std::vector<meshkernel::Point>> generatedPoints;
    std::vector<meshkernel::Point> refinedPolygon;
    polygons.RefinePolygonPart(0, 0, 1.0, refinedPolygon);

    ASSERT_EQ(13, refinedPolygon.size());
    const double tolerance = 1e-5;

    ASSERT_NEAR(0, refinedPolygon[0].x, tolerance);
    ASSERT_NEAR(1, refinedPolygon[1].x, tolerance);
    ASSERT_NEAR(2, refinedPolygon[2].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon[3].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon[4].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon[5].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon[6].x, tolerance);
    ASSERT_NEAR(2, refinedPolygon[7].x, tolerance);
    ASSERT_NEAR(1, refinedPolygon[8].x, tolerance);
    ASSERT_NEAR(0, refinedPolygon[9].x, tolerance);
    ASSERT_NEAR(0, refinedPolygon[10].x, tolerance);
    ASSERT_NEAR(0, refinedPolygon[11].x, tolerance);
    ASSERT_NEAR(0, refinedPolygon[12].x, tolerance);

    ASSERT_NEAR(0, refinedPolygon[0].y, tolerance);
    ASSERT_NEAR(0, refinedPolygon[1].y, tolerance);
    ASSERT_NEAR(0, refinedPolygon[2].y, tolerance);
    ASSERT_NEAR(0, refinedPolygon[3].y, tolerance);
    ASSERT_NEAR(1, refinedPolygon[4].y, tolerance);
    ASSERT_NEAR(2, refinedPolygon[5].y, tolerance);
    ASSERT_NEAR(3, refinedPolygon[6].y, tolerance);
    ASSERT_NEAR(3, refinedPolygon[7].y, tolerance);
    ASSERT_NEAR(3, refinedPolygon[8].y, tolerance);
    ASSERT_NEAR(3, refinedPolygon[9].y, tolerance);
    ASSERT_NEAR(2, refinedPolygon[10].y, tolerance);
    ASSERT_NEAR(1, refinedPolygon[11].y, tolerance);
    ASSERT_NEAR(0, refinedPolygon[12].y, tolerance);
}

TEST(Polygons, RefinePolygonOneSide)
{
    // Prepare
    std::vector<meshkernel::Point> nodes;

    nodes.push_back({0, 0});
    nodes.push_back({3, 0});
    nodes.push_back({3, 3});
    nodes.push_back({0, 3});
    nodes.push_back({0, 0});

    meshkernel::Polygons polygons(nodes, meshkernel::Projections::cartesian);

    // Execute
    std::vector<std::vector<meshkernel::Point>> generatedPoints;
    std::vector<meshkernel::Point> refinedPolygon;
    polygons.RefinePolygonPart(0, 1, 1.0, refinedPolygon);

    ASSERT_EQ(7, refinedPolygon.size());
    const double tolerance = 1e-5;

    ASSERT_NEAR(0.0, refinedPolygon[0].x, tolerance);
    ASSERT_NEAR(1.0, refinedPolygon[1].x, tolerance);
    ASSERT_NEAR(2.0, refinedPolygon[2].x, tolerance);
    ASSERT_NEAR(3.0, refinedPolygon[3].x, tolerance);
    ASSERT_NEAR(3.0, refinedPolygon[4].x, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[5].x, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[6].x, tolerance);

    ASSERT_NEAR(0.0, refinedPolygon[0].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[1].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[2].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[3].y, tolerance);
    ASSERT_NEAR(3.0, refinedPolygon[4].y, tolerance);
    ASSERT_NEAR(3.0, refinedPolygon[5].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[6].y, tolerance);
}

TEST(Polygons, RefinePolygonLongerSquare)
{
    // Prepare
    std::vector<meshkernel::Point> nodes;

    nodes.push_back({0, 0});
    nodes.push_back({3, 0});
    nodes.push_back({3, 3});
    nodes.push_back({3.5, 0});
    nodes.push_back({0, 0});

    meshkernel::Polygons polygons(nodes, meshkernel::Projections::cartesian);

    // Execute
    std::vector<std::vector<meshkernel::Point>> generatedPoints;
    std::vector<meshkernel::Point> refinedPolygon;
    polygons.RefinePolygonPart(0, 0, 1.0, refinedPolygon);

    ASSERT_EQ(15, refinedPolygon.size());
    const double tolerance = 1e-5;

    ASSERT_NEAR(0.0, refinedPolygon[0].x, tolerance);
    ASSERT_NEAR(1.0, refinedPolygon[1].x, tolerance);
    ASSERT_NEAR(2.0, refinedPolygon[2].x, tolerance);
    ASSERT_NEAR(3.0, refinedPolygon[3].x, tolerance);
    ASSERT_NEAR(3.0, refinedPolygon[4].x, tolerance);
    ASSERT_NEAR(3.0, refinedPolygon[5].x, tolerance);
    ASSERT_NEAR(3.0, refinedPolygon[6].x, tolerance);
    ASSERT_NEAR(3.1643989873053573, refinedPolygon[7].x, tolerance);
    ASSERT_NEAR(3.3287979746107146, refinedPolygon[8].x, tolerance);
    ASSERT_NEAR(3.4931969619160719, refinedPolygon[9].x, tolerance);
    ASSERT_NEAR(3.5, refinedPolygon[10].x, tolerance);
    ASSERT_NEAR(2.5, refinedPolygon[11].x, tolerance);
    ASSERT_NEAR(1.5, refinedPolygon[12].x, tolerance);
    ASSERT_NEAR(0.5, refinedPolygon[13].x, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[14].x, tolerance);

    ASSERT_NEAR(0.0, refinedPolygon[0].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[1].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[2].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[3].y, tolerance);
    ASSERT_NEAR(1.0, refinedPolygon[4].y, tolerance);
    ASSERT_NEAR(2.0, refinedPolygon[5].y, tolerance);
    ASSERT_NEAR(3.0, refinedPolygon[6].y, tolerance);
    ASSERT_NEAR(2.0136060761678563, refinedPolygon[7].y, tolerance);
    ASSERT_NEAR(1.0272121523357125, refinedPolygon[8].y, tolerance);
    ASSERT_NEAR(0.040818228503568754, refinedPolygon[9].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[10].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[11].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[12].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[13].y, tolerance);
    ASSERT_NEAR(0.0, refinedPolygon[14].y, tolerance);
}

TEST(Polygons, OffsetCopy)
{
    std::vector<meshkernel::Point> nodes;
    nodes.push_back({296.752472, 397.879639});
    nodes.push_back({294.502472, 256.128204});
    nodes.push_back({578.754211, 244.128082});
    nodes.push_back({587.754272, 400.129639});
    nodes.push_back({308.002533, 397.879639});
    nodes.push_back({296.752472, 397.879639});

    meshkernel::Polygons polygon(nodes, meshkernel::Projections::cartesian);

    meshkernel::Polygons newPolygon;
    double distance = 10.0;
    bool innerAndOuter = false;
    polygon.OffsetCopy(distance, innerAndOuter, newPolygon);

    const double tolerance = 1e-5;

    ASSERT_NEAR(newPolygon.m_nodes[0].x, 286.75373149966771, tolerance);
    ASSERT_NEAR(newPolygon.m_nodes[1].x, 284.34914611880089, tolerance);
    ASSERT_NEAR(newPolygon.m_nodes[2].x, 588.17047010011993, tolerance);
    ASSERT_NEAR(newPolygon.m_nodes[3].x, 598.35275776004642, tolerance);
    ASSERT_NEAR(newPolygon.m_nodes[4].x, 307.96231942308754, tolerance);
    ASSERT_NEAR(newPolygon.m_nodes[5].x, 296.75247200000001, tolerance);

    ASSERT_NEAR(newPolygon.m_nodes[0].y, 398.03834755999270, tolerance);
    ASSERT_NEAR(newPolygon.m_nodes[1].y, 246.54793497426144, tolerance);
    ASSERT_NEAR(newPolygon.m_nodes[2].y, 233.72165300742589, tolerance);
    ASSERT_NEAR(newPolygon.m_nodes[3].y, 410.21520441451258, tolerance);
    ASSERT_NEAR(newPolygon.m_nodes[4].y, 407.87963900000000, tolerance);
    ASSERT_NEAR(newPolygon.m_nodes[5].y, 407.87963900000000, tolerance);
}
