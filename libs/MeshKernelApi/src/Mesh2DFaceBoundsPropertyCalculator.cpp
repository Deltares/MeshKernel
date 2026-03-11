//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2026.
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

#include "MeshKernelApi/Mesh2DFaceBoundsPropertyCalculator.hpp"
#include "MeshKernelApi/PropertyCalculator.hpp"
#include "MeshKernelApi/State.hpp"

#include "MeshKernel/Mesh2DFaceBounds.hpp"

bool meshkernelapi::Mesh2DFaceBoundsPropertyCalculator::IsValid(const MeshKernelState& state, const meshkernel::Location location) const
{
    return state.m_mesh2d != nullptr && state.m_mesh2d->GetNumNodes() > 0 && location == meshkernel::Location::Faces;
}

void meshkernelapi::Mesh2DFaceBoundsPropertyCalculator::Calculate(const MeshKernelState& state, const meshkernel::Location location, const GeometryList& geometryList) const
{
    if (geometryList.num_coordinates < static_cast<int>(meshkernel::constants::geometric::maximumNumberOfNodesPerFace * state.m_mesh2d->GetNumFaces()))
    {
        throw meshkernel::ConstraintError("GeometryList with wrong dimensions, {} must be greater than or equal to {}",
                                          geometryList.num_coordinates, Size(state, location));
    }

    std::vector<meshkernel::Point> faceBounds(meshkernel::algo::Mesh2DFaceBounds::Compute(*state.m_mesh2d));

    size_t size = static_cast<size_t>(Size(state, location));
    std::span<double> xCoord(geometryList.coordinates_x, size);
    std::span<double> yCoord(geometryList.coordinates_y, size);

    for (size_t i = 0; i < faceBounds.size(); ++i)
    {
        xCoord[i] = faceBounds[i].x;
        yCoord[i] = faceBounds[i].y;
    }
}

int meshkernelapi::Mesh2DFaceBoundsPropertyCalculator::Size(const MeshKernelState& state, const meshkernel::Location location) const
{
    int size = -1;

    if (location == meshkernel::Location::Faces)
    {
        size = meshkernel::constants::geometric::maximumNumberOfNodesPerFace * static_cast<int>(state.m_mesh2d->GetNumFaces());
    }

    return size;
}
