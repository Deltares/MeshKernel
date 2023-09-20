#include <gtest/gtest.h>

#include "MeshKernel/Parameters.hpp"

using namespace meshkernel;

TEST(ParametersTests, MakeGridParameters)
{

    MakeGridParameters parameters;

    // default parameters should be valid
    EXPECT_NO_THROW(CheckMakeGridParameters(parameters));

    // reset and assign an invalid value, one struct member at a time

    parameters = MakeGridParameters();
    parameters.num_columns = -10;
    EXPECT_THROW(CheckMakeGridParameters(parameters), RangeError);

    parameters = MakeGridParameters();
    parameters.num_rows = 0;
    EXPECT_THROW(CheckMakeGridParameters(parameters), RangeError);

    parameters = MakeGridParameters();
    parameters.angle = 100.0;
    EXPECT_THROW(CheckMakeGridParameters(parameters), RangeError);

    parameters = MakeGridParameters();
    parameters.block_size_x = 0.0;
    EXPECT_THROW(CheckMakeGridParameters(parameters), RangeError);

    parameters = MakeGridParameters();
    parameters.block_size_y = -1.0;
    EXPECT_THROW(CheckMakeGridParameters(parameters), RangeError);
}

TEST(ParametersTests, CurvilinearParameters)
{
    CurvilinearParameters parameters;

    // default parameters should be valid
    EXPECT_NO_THROW(parameters);

    // reset and assign an invalid value, one struct member at a time

    parameters = CurvilinearParameters();
    parameters.m_refinement = 0;
    EXPECT_THROW(CheckCurvilinearParameters(parameters), RangeError);

    parameters = CurvilinearParameters();
    parameters.n_refinement = -1;
    EXPECT_THROW(CheckCurvilinearParameters(parameters), RangeError);

    parameters = CurvilinearParameters();
    parameters.smoothing_iterations = -1;
    EXPECT_THROW(CheckCurvilinearParameters(parameters), RangeError);

    parameters = CurvilinearParameters();
    parameters.smoothing_parameter = 1.5;
    EXPECT_THROW(CheckCurvilinearParameters(parameters), RangeError);

    parameters = CurvilinearParameters();
    parameters.attraction_parameter = -1.0;
    EXPECT_THROW(CheckCurvilinearParameters(parameters), RangeError);
}

TEST(ParametersTests, SplinesToCurvilinearParameters)
{
    SplinesToCurvilinearParameters parameters;

    // default parameters should be valid
    EXPECT_NO_THROW(parameters);

    // reset and assign an invalid value, one struct member at a time

    parameters = SplinesToCurvilinearParameters();
    parameters.aspect_ratio = -1.5;
    EXPECT_THROW(CheckSplinesToCurvilinearParameters(parameters), RangeError);

    parameters = SplinesToCurvilinearParameters();
    parameters.aspect_ratio_grow_factor = 0.0;
    EXPECT_THROW(CheckSplinesToCurvilinearParameters(parameters), RangeError);

    parameters = SplinesToCurvilinearParameters();
    parameters.average_width = -100.0;
    EXPECT_THROW(CheckSplinesToCurvilinearParameters(parameters), RangeError);

    parameters = SplinesToCurvilinearParameters();
    parameters.curvature_adapted_grid_spacing = 3;
    EXPECT_THROW(CheckSplinesToCurvilinearParameters(parameters), RangeError);

    parameters = SplinesToCurvilinearParameters();
    parameters.grow_grid_outside = -1;
    EXPECT_THROW(CheckSplinesToCurvilinearParameters(parameters), RangeError);

    parameters = SplinesToCurvilinearParameters();
    parameters.maximum_num_faces_in_uniform_part = 0;
    EXPECT_THROW(CheckSplinesToCurvilinearParameters(parameters), RangeError);

    parameters = SplinesToCurvilinearParameters();
    parameters.nodes_on_top_of_each_other_tolerance = -1.0e-6;
    EXPECT_THROW(CheckSplinesToCurvilinearParameters(parameters), RangeError);

    parameters = SplinesToCurvilinearParameters();
    parameters.min_cosine_crossing_angles = 0.0;
    EXPECT_THROW(CheckSplinesToCurvilinearParameters(parameters), RangeError);

    parameters = SplinesToCurvilinearParameters();
    parameters.check_front_collisions = 2;
    EXPECT_THROW(CheckSplinesToCurvilinearParameters(parameters), RangeError);

    parameters = SplinesToCurvilinearParameters();
    parameters.remove_skinny_triangles = 6;
    EXPECT_THROW(CheckSplinesToCurvilinearParameters(parameters), RangeError);
}

TEST(ParametersTests, MeshRefinementParameters)
{
    MeshRefinementParameters parameters;

    // default parameters should be valid
    EXPECT_NO_THROW(parameters);

    // reset and assign an invalid value, one struct member at a time

    parameters = MeshRefinementParameters();
    parameters.max_num_refinement_iterations = -3;
    EXPECT_THROW(CheckMeshRefinementParameters(parameters), RangeError);

    parameters = MeshRefinementParameters();
    parameters.refine_intersected = 2;
    EXPECT_THROW(CheckMeshRefinementParameters(parameters), RangeError);

    parameters = MeshRefinementParameters();
    parameters.use_mass_center_when_refining = 9;
    EXPECT_THROW(CheckMeshRefinementParameters(parameters), RangeError);

    parameters = MeshRefinementParameters();
    parameters.min_edge_size = -1.0e-3;
    EXPECT_THROW(CheckMeshRefinementParameters(parameters), RangeError);

    parameters = MeshRefinementParameters();
    parameters.refinement_type = 0;
    EXPECT_THROW(CheckMeshRefinementParameters(parameters), RangeError);

    parameters = MeshRefinementParameters();
    parameters.connect_hanging_nodes = 2;
    EXPECT_THROW(CheckMeshRefinementParameters(parameters), RangeError);

    parameters = MeshRefinementParameters();
    parameters.account_for_samples_outside = 3;
    EXPECT_THROW(CheckMeshRefinementParameters(parameters), RangeError);

    parameters = MeshRefinementParameters();
    parameters.smoothing_iterations = -9;
    EXPECT_THROW(CheckMeshRefinementParameters(parameters), RangeError);

    parameters = MeshRefinementParameters();
    parameters.max_courant_time = 0.0;
    EXPECT_THROW(CheckMeshRefinementParameters(parameters), RangeError);

    parameters = MeshRefinementParameters();
    parameters.directional_refinement = -1;
    EXPECT_THROW(CheckMeshRefinementParameters(parameters), RangeError);
}

TEST(ParametersTests, OrthogonalizationParameters)
{
    OrthogonalizationParameters parameters;

    // default parameters should be valid
    EXPECT_NO_THROW(parameters);

    // reset and assign an invalid value, one struct member at a time

    parameters = OrthogonalizationParameters();
    parameters.outer_iterations = 0;
    EXPECT_THROW(CheckOrthogonalizationParameters(parameters), RangeError);

    parameters = OrthogonalizationParameters();
    parameters.boundary_iterations = 0;
    EXPECT_THROW(CheckOrthogonalizationParameters(parameters), RangeError);

    parameters = OrthogonalizationParameters();
    parameters.inner_iterations = 0;
    EXPECT_THROW(CheckOrthogonalizationParameters(parameters), RangeError);

    parameters = OrthogonalizationParameters();
    parameters.orthogonalization_to_smoothing_factor = 1.1;
    EXPECT_THROW(CheckOrthogonalizationParameters(parameters), RangeError);

    parameters = OrthogonalizationParameters();
    parameters.orthogonalization_to_smoothing_factor_at_boundary = -0.2;
    EXPECT_THROW(CheckOrthogonalizationParameters(parameters), RangeError);

    parameters = OrthogonalizationParameters();
    parameters.areal_to_angle_smoothing_factor = 1.000001;
    EXPECT_THROW(CheckOrthogonalizationParameters(parameters), RangeError);
}