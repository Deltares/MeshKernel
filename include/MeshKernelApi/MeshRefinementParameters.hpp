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

namespace meshkernelapi
{
    /// @brief A struct used to describe the mesh parameters in a C-compatible manner
    struct MeshRefinementParameters
    {
        /// @brief Maximum number of refinement iterations, set to 1 if only one refinement is wanted (10)
        int max_num_refinement_iterations;

        /// @brief Whether to compute faces intersected by polygon (yes=1/no=0)
        int refine_intersected;

        /// Whether to use the mass center when splitting a face in the refinement process (yes=1/no=0)
        int use_mass_center_when_refining;

        /// @brief Minimum cell size
        double min_face_size;

        /// @brief Refinement criterion type
        int refinement_type;

        /// @brief Connect hanging nodes at the end of the iteration, 1 yes or 0 no
        int connect_hanging_nodes = 1;

        /// @brief Take samples outside face into account , 1 yes 0 no
        int account_for_samples_outside;
    };
} // namespace meshkernelapi
