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

#pragma once

#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Parameters.hpp"
#include "MeshKernel/SampleInterpolator.hpp"

#include "MeshKernelApi/GeometryList.hpp"
#include "MeshKernelApi/PropertyCalculator.hpp"

namespace meshkernelapi
{

    /// @brief Calculator for the edge lengths for a mesh
    class EdgeLengthPropertyCalculator : public PropertyCalculator
    {
    public:
        /// @brief Determine is the calculator can compute the desired results correctly.
        ///
        /// This has a default of checking that the mesh2d is valid and the location is at edges
        bool IsValid(const MeshKernelState& state, const meshkernel::Location location) const override;

        /// @brief Calculate the edge-length for a mesh
        ///
        /// \note This calculator is for mesh edges only
        void Calculate(const MeshKernelState& state, const meshkernel::Location location, const GeometryList& geometryList) const override;

        /// @brief Determine the size of the edge-length vector required
        int Size(const MeshKernelState& state, const meshkernel::Location location) const override;
    };

} // namespace meshkernelapi
