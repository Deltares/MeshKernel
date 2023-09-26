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

#include "MeshKernel/MeshRefinementType.hpp"
#include "MeshKernel/RangeCheck.hpp"

namespace meshkernel
{
    /// @brief This struct describes the necessary parameters to create a new curvilinear grid in a C-compatible manner
    ///
    /// @see mkernel_make_mesh
    struct MakeGridParameters
    {
        /// @brief The number of columns in x direction
        int num_columns = 3;

        /// @brief The number of columns in y direction
        int num_rows = 3;

        /// @brief The grid angle
        double angle = 0.0;

        /// @brief The x coordinate of the origin, located at the bottom left corner
        double origin_x = 0.0;

        /// @brief The y coordinate of the origin, located at the bottom left corner
        double origin_y = 0.0;

        /// @brief The grid block size in x dimension, used only for squared grids
        double block_size_x = 10.0;

        /// @brief The grid block size in y dimension, used only for squared grids
        double block_size_y = 10.0;

        /// @brief The x coordinate of the upper right corner
        double upper_right_x = 0.0;

        /// @brief The y coordinate of the upper right corner
        double upper_right_y = 0.0;
    };

    inline static void CheckMakeGridParameters(MakeGridParameters const& parameters)
    {
        range_check::CheckGreater(parameters.num_columns, 0, "Number of columns");
        range_check::CheckGreater(parameters.num_rows, 0, "Number of rows");
        range_check::CheckInClosedInterval(parameters.angle, {-90.0, 90.0}, "Grid angle");
        range_check::CheckGreater(parameters.block_size_x, 0.0, "X block size");
        range_check::CheckGreater(parameters.block_size_y, 0.0, "Y block size");
    }

    /// @brief A struct used to describe parameters for generating a curvilinear grid in a C-compatible manner
    struct CurvilinearParameters
    {
        /// @brief M-refinement factor for regular grid generation (mfacmax)
        int m_refinement = 2000;

        /// @brief N-refinement factor for regular grid generation (nfacmax)
        int n_refinement = 40;

        /// @brief Nr. of inner iterations in regular grid smoothing
        int smoothing_iterations = 10;

        /// @brief Smoothing parameter
        double smoothing_parameter = 0.5;

        /// @brief Attraction/repulsion parameter
        double attraction_parameter = 0.0;
    };

    inline static void CheckCurvilinearParameters(CurvilinearParameters const& parameters)
    {
        range_check::CheckGreater(parameters.m_refinement, 0, " M-refinement factor");
        range_check::CheckGreater(parameters.n_refinement, 0, "N-refinement factor");
        range_check::CheckGreater(parameters.smoothing_iterations, 0, "Smoothing iterations");
        range_check::CheckInClosedInterval(parameters.smoothing_parameter, {0.0, 1.0}, "Smoothing parameter"); // CHECK ME
        range_check::CheckGreaterEqual(parameters.attraction_parameter, 0.0, "Attraction parameter");          // CHECK ME
    }

    /// @brief A struct used to describe the spline to curvilinear grid parameters in a C-compatible manner
    struct SplinesToCurvilinearParameters
    {
        /// @brief Aspect ratio (mfacmax)
        double aspect_ratio = 0.1;

        /// @brief Grow factor of aspect ratio
        double aspect_ratio_grow_factor = 1.1;

        /// @brief Average mesh width on center spline
        double average_width = 500.0;

        /// @brief Curvature adapted grid spacing, 1 or not 0
        int curvature_adapted_grid_spacing = 1;

        /// @brief Grow the grid outside the prescribed grid height
        int grow_grid_outside = 0;

        /// @brief Maximum number of layers in the uniform part
        int maximum_num_faces_in_uniform_part = 5;

        /// @brief On-top-of-each-other tolerance
        double nodes_on_top_of_each_other_tolerance = 0.0001;

        /// @brief Minimum allowed absolute value of crossing-angle cosine
        double min_cosine_crossing_angles = 0.95;

        /// @brief Check for collisions with other parts of the front, 1 or not 0
        int check_front_collisions = 0;

        /// @brief Check for collisions with other parts of the front
        int remove_skinny_triangles = 1;
    };

    inline static void CheckSplinesToCurvilinearParameters(SplinesToCurvilinearParameters const& parameters)
    {
        range_check::CheckGreater(parameters.aspect_ratio, 0.0, "Aspect ratio");
        range_check::CheckGreater(parameters.aspect_ratio_grow_factor, 0.0, "Aspect ratio grow factor");
        range_check::CheckGreater(parameters.average_width, 0.0, "Average width");
        range_check::CheckOneOf(parameters.curvature_adapted_grid_spacing, {0, 1}, "Curvature adapted grid spacing");
        range_check::CheckOneOf(parameters.grow_grid_outside, {0, 1}, "Grow grid outside");
        range_check::CheckGreater(parameters.maximum_num_faces_in_uniform_part, 0, "Max number of faces in uniform part");
        range_check::CheckGreater(parameters.nodes_on_top_of_each_other_tolerance, 0.0, "nodes on top of each other tolerance");
        range_check::CheckGreater(parameters.min_cosine_crossing_angles, 0.0, "Min cosine crossing angles");
        range_check::CheckOneOf(parameters.check_front_collisions, {0, 1}, "Check front collisions");
        range_check::CheckOneOf(parameters.remove_skinny_triangles, {0, 1}, "Remove skinny triangles");
    }

    /// @brief A struct used to describe the mesh refinement parameters in a C-compatible manner
    struct MeshRefinementParameters
    {
        /// @brief Maximum number of refinement iterations, set to 1 if only one refinement is wanted
        int max_num_refinement_iterations = 10;

        /// @brief Whether to compute faces intersected by polygon (yes=1/no=0)
        int refine_intersected = 0;

        /// Whether to use the mass center when splitting a face in the refinement process (yes=1/no=0)
        int use_mass_center_when_refining = 1;

        /// @brief Minimum edge size in meters
        double min_edge_size = 0.5;

        /// @brief Refinement criterion type
        int refinement_type = 2;

        /// @brief Connect hanging nodes at the end of the iteration, 1 yes or 0 no
        int connect_hanging_nodes = 1;

        /// @brief Take samples outside face into account , 1 yes 0 no
        int account_for_samples_outside = 0;

        /// @brief The number of smoothing iterations
        int smoothing_iterations = 5;

        /// @brief Maximum courant time in seconds
        double max_courant_time = 120.0;

        /// @brief Directional refinement, cannot be used when the number of smoothing iterations is larger than 0
        int directional_refinement = 0;
    };

    inline static void CheckMeshRefinementParameters(MeshRefinementParameters const& parameters)
    {
        range_check::CheckGreater(parameters.max_num_refinement_iterations, 0, "Max num refinement iterations");
        range_check::CheckOneOf(parameters.refine_intersected, {0, 1}, "Refine intersected");
        range_check::CheckOneOf(parameters.use_mass_center_when_refining, {0, 1}, "Use mass center when refining");
        range_check::CheckGreater(parameters.min_edge_size, 0.0, "Min edge size");
        range_check::CheckOneOf(parameters.refinement_type, MeshRefinementType::ValidValues(), "Refinement type");
        range_check::CheckOneOf(parameters.connect_hanging_nodes, {0, 1}, "Connect hanging nodes");
        range_check::CheckOneOf(parameters.account_for_samples_outside, {0, 1}, "Account for samples outside");
        range_check::CheckGreaterEqual(parameters.smoothing_iterations, 0, "Smoothing iterations");
        range_check::CheckGreater(parameters.max_courant_time, 0.0, "Max courant time");
        range_check::CheckOneOf(parameters.directional_refinement, {0, 1}, "Directional refinement");
    }

    /// @brief A struct used to describe the orthogonalization parameters in a C-compatible manner
    struct OrthogonalizationParameters
    {
        /// @brief Number of outer iterations in orthogonalization. Increase this parameter for complex grids
        int outer_iterations = 2;

        /// @brief Number of boundary iterations in grid/net orthogonalization within itatp
        int boundary_iterations = 25;

        /// @brief Number of inner iterations in grid/net orthogonalization within itbnd
        int inner_iterations = 25;

        /// @brief Factor from 0 to 1. between grid smoothing and grid orthogonality
        double orthogonalization_to_smoothing_factor = 0.975;

        /// @brief Minimum ATPF on the boundary
        double orthogonalization_to_smoothing_factor_at_boundary = 1.0;

        /// @brief Factor between smoother 1d0 and area-homogenizer 0d0
        double areal_to_angle_smoothing_factor = 1.0;
    };

    inline static void CheckOrthogonalizationParameters(OrthogonalizationParameters const& parameters)
    {
        range_check::CheckGreater(parameters.outer_iterations, 0, "Outer iterations");
        range_check::CheckGreater(parameters.boundary_iterations, 0, "Boundary iterations");
        range_check::CheckGreater(parameters.inner_iterations, 0, "Inner iterations");
        range_check::CheckInClosedInterval(parameters.orthogonalization_to_smoothing_factor, {0.0, 1.0}, "Orthogonalization-to-smoothing_factor");
        range_check::CheckInClosedInterval(parameters.orthogonalization_to_smoothing_factor_at_boundary, {0.0, 1.0}, "orthogonalization-to-smoothing factor at boundary");
        range_check::CheckInClosedInterval(parameters.areal_to_angle_smoothing_factor, {0.0, 1.0}, "area to angle smoothing factor");
    }

} // namespace meshkernel
