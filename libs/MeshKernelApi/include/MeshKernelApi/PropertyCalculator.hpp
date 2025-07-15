//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2024.
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

namespace meshkernelapi
{

    /// @brief Forward declaration of MeshKernelState
    struct MeshKernelState;

    /// @brief Base class for calculating properties for a mesh
    class PropertyCalculator
    {
    public:
        /// @brief Destructor
        virtual ~PropertyCalculator() = default;

        /// @brief Determine is the calculator can compute the desired results correctly.
        virtual bool IsValid(const MeshKernelState& state, const meshkernel::Location location) const = 0;

        /// @brief Calculate the property
        virtual void Calculate(const MeshKernelState& state, const meshkernel::Location location, const GeometryList& geometryList) const = 0;

        /// @brief Determine the size of the vector required to store the calculated properties
        virtual int Size(const MeshKernelState& state, const meshkernel::Location location) const = 0;
    };

} // namespace meshkernelapi
