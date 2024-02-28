#include <gtest/gtest.h>

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/FlipEdges.hpp>
#include <MeshKernel/LandBoundaries.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Polygons.hpp>
#include <TestUtils/Definitions.hpp>
#include <TestUtils/MakeMeshes.hpp>

TEST(FlipEdges, FlipEdgesWithLandBoundary)
{
    // 1 Setup
    auto mesh = MakeRectangularMeshForTesting(3, 3, 10, meshkernel::Projection::cartesian, {0.0, 0.0});

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    // set landboundaries
    auto polygon = meshkernel::Polygons();
    std::vector<meshkernel::Point> landBoundary{{-1.369282, 21.249086},
                                                {20.885406, 21.539995},
                                                {meshkernel::constants::missing::doubleValue, meshkernel::constants::missing::doubleValue}};

    auto landBoundaries = meshkernel::LandBoundaries(landBoundary, *mesh, polygon);

    // execute flipedges
    meshkernel::FlipEdges flipEdges(*mesh, landBoundaries, true, true);

    auto undoAction = flipEdges.Compute();

    // check the values
    ASSERT_EQ(16, mesh->GetNumEdges());

    // Test the undo action has been computed correctly
    undoAction->Restore();
    // Recompute faces
    mesh->Administrate();

    ASSERT_EQ(originalNodes.size(), mesh->GetNumValidNodes());
    ASSERT_EQ(originalEdges.size(), mesh->GetNumValidEdges());

    size_t count = 0;

    for (size_t i = 0; i < mesh->Nodes().size(); ++i)
    {
        if (mesh->Node(i).IsValid())
        {
            // Check valid nodes
            EXPECT_EQ(originalNodes[count].x, mesh->Node(i).x);
            EXPECT_EQ(originalNodes[count].y, mesh->Node(i).y);
            ++count;
        }
    }

    count = 0;

    for (size_t i = 0; i < mesh->Edges().size(); ++i)
    {
        if (mesh->IsValidEdge(i))
        {
            EXPECT_EQ(originalEdges[count].first, mesh->GetEdge(i).first);
            EXPECT_EQ(originalEdges[count].second, mesh->GetEdge(i).second);
            ++count;
        }
    }
}

TEST(FlipEdges, FlipEdgesMediumTriangularMesh)
{
    // 1 Setup
    auto mesh = ReadLegacyMesh2DFromFile(TEST_FOLDER + "/data/TestOrthogonalizationMediumTriangularGrid_net.nc");

    const std::vector<meshkernel::Point> originalNodes(mesh->Nodes());
    const std::vector<meshkernel::Edge> originalEdges(mesh->Edges());

    // set landboundaries
    auto polygon = meshkernel::Polygons();

    std::vector<meshkernel::Point> landBoundary;
    auto landBoundaries = meshkernel::LandBoundaries(landBoundary, *mesh, polygon);

    // execute flipedges
    meshkernel::FlipEdges flipEdges(*mesh, landBoundaries, true, false);

    auto undoAction = flipEdges.Compute();

    // get the number of edges
    ASSERT_EQ(697, mesh->GetNumEdges());

    // check the values of flipped edges
    ASSERT_EQ(183, mesh->GetEdge(14).first);
    ASSERT_EQ(227, mesh->GetEdge(14).second);

    ASSERT_EQ(58, mesh->GetEdge(33).first);
    ASSERT_EQ(141, mesh->GetEdge(33).second);

    ASSERT_EQ(147, mesh->GetEdge(46).first);
    ASSERT_EQ(145, mesh->GetEdge(46).second);

    ASSERT_EQ(147, mesh->GetEdge(49).first);
    ASSERT_EQ(148, mesh->GetEdge(49).second);

    ASSERT_EQ(242, mesh->GetEdge(68).first);
    ASSERT_EQ(148, mesh->GetEdge(68).second);

    // Test the undo action has been computed correctly
    undoAction->Restore();
    // Recompute faces
    mesh->Administrate();

    ASSERT_EQ(originalNodes.size(), mesh->GetNumValidNodes());
    ASSERT_EQ(originalEdges.size(), mesh->GetNumValidEdges());

    size_t count = 0;

    for (size_t i = 0; i < mesh->Nodes().size(); ++i)
    {
        if (mesh->Node(i).IsValid())
        {
            // Check valid nodes
            EXPECT_EQ(originalNodes[count].x, mesh->Node(i).x);
            EXPECT_EQ(originalNodes[count].y, mesh->Node(i).y);
            ++count;
        }
    }

    count = 0;

    for (size_t i = 0; i < mesh->Edges().size(); ++i)
    {
        if (mesh->IsValidEdge(i))
        {
            EXPECT_EQ(originalEdges[count].first, mesh->GetEdge(i).first);
            EXPECT_EQ(originalEdges[count].second, mesh->GetEdge(i).second);
            ++count;
        }
    }
}
