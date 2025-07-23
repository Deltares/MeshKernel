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

#include "MeshKernelApi/OrthogonalityPropertyCalculator.hpp"
#include "MeshKernelApi/PropertyCalculator.hpp"
#include "MeshKernelApi/State.hpp"

#include "MeshKernel/MeshOrthogonality.hpp"

#include <algorithm>
#include <functional>

bool meshkernelapi::OrthogonalityPropertyCalculator::IsValid(const MeshKernelState& state, const meshkernel::Location location) const
{
    return state.m_mesh2d != nullptr && state.m_mesh2d->GetNumNodes() > 0 && location == meshkernel::Location::Edges;
}

void meshkernelapi::OrthogonalityPropertyCalculator::Calculate(const MeshKernelState& state, const meshkernel::Location location, const GeometryList& geometryList) const
{

    if (geometryList.num_coordinates < static_cast<int>(state.m_mesh2d->GetNumEdges()))
    {
        throw meshkernel::ConstraintError("GeometryList with wrong dimensions, {} must be greater than or equal to {}",
                                          geometryList.num_coordinates, Size(state, location));
    }

    std::span<double> orthogonality(geometryList.values, geometryList.num_coordinates);
    meshkernel::MeshOrthogonality::Compute(*state.m_mesh2d, orthogonality);
}

int meshkernelapi::OrthogonalityPropertyCalculator::Size(const MeshKernelState& state, const meshkernel::Location location [[maybe_unused]]) const
{
    return static_cast<int>(state.m_mesh2d->GetNumEdges());
}
