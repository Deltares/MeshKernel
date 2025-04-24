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

#include "MeshKernelApi/PropertyCalculator.hpp"
#include "MeshKernelApi/State.hpp"

#include "MeshKernel/MeshEdgeLength.hpp"
#include "MeshKernel/MeshOrthogonality.hpp"
#include "MeshKernel/SampleAveragingInterpolator.hpp"
#include "MeshKernel/SampleTriangulationInterpolator.hpp"

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

bool meshkernelapi::EdgeLengthPropertyCalculator::IsValid(const MeshKernelState& state, const meshkernel::Location location) const
{
    return state.m_mesh2d != nullptr && state.m_mesh2d->GetNumNodes() > 0 && location == meshkernel::Location::Edges;
}

void meshkernelapi::EdgeLengthPropertyCalculator::Calculate(const MeshKernelState& state, const meshkernel::Location location, const GeometryList& geometryList) const
{

    if (static_cast<size_t>(geometryList.num_coordinates) < state.m_mesh2d->GetNumEdges())
    {
        throw meshkernel::ConstraintError("GeometryList with wrong dimensions, {} must be greater than or equal to {}",
                                          geometryList.num_coordinates, Size(state, location));
    }

    std::span<double> edgeLengths(geometryList.values, state.m_mesh2d->GetNumEdges());
    meshkernel::algo::ComputeMeshEdgeLength(*state.m_mesh2d, edgeLengths);
}

int meshkernelapi::EdgeLengthPropertyCalculator::Size(const MeshKernelState& state, const meshkernel::Location location [[maybe_unused]]) const
{
    return static_cast<int>(state.m_mesh2d->GetNumEdges());
}

meshkernelapi::InterpolatedSamplePropertyCalculator::InterpolatedSamplePropertyCalculator(const GeometryList& sampleData,
                                                                                          const meshkernel::Projection projection,
                                                                                          const meshkernel::InterpolationParameters& interpolationParameters,
                                                                                          const int propertyId)
    : m_projection(projection),
      m_propertyId(propertyId)
{
    std::span<const double> xNodes(sampleData.coordinates_x, sampleData.num_coordinates);
    std::span<const double> yNodes(sampleData.coordinates_y, sampleData.num_coordinates);

    if (interpolationParameters.interpolation_type == 0)
    {
        m_sampleInterpolator = std::make_unique<meshkernel::SampleTriangulationInterpolator>(xNodes, yNodes, m_projection);
    }
    else if (interpolationParameters.interpolation_type == 1)
    {
        // Need to pass from api.
        m_sampleInterpolator = std::make_unique<meshkernel::SampleAveragingInterpolator>(xNodes, yNodes, m_projection, interpolationParameters);
    }

    std::span<const double> dataSamples(sampleData.values, sampleData.num_coordinates);
    m_sampleInterpolator->SetData(m_propertyId, dataSamples);
}

bool meshkernelapi::InterpolatedSamplePropertyCalculator::IsValid(const MeshKernelState& state, const meshkernel::Location location [[maybe_unused]]) const
{
    return state.m_mesh2d != nullptr &&
           state.m_mesh2d->GetNumNodes() > 0 &&
           m_sampleInterpolator->Contains(m_propertyId) &&
           m_projection == state.m_projection;
}

void meshkernelapi::InterpolatedSamplePropertyCalculator::Calculate(const MeshKernelState& state, const meshkernel::Location location, const GeometryList& geometryList) const
{
    std::span<double> interpolatedSampleData(geometryList.values, geometryList.num_coordinates);
    m_sampleInterpolator->Interpolate(m_propertyId, *state.m_mesh2d, location, interpolatedSampleData);
}

int meshkernelapi::InterpolatedSamplePropertyCalculator::Size(const MeshKernelState& state, const meshkernel::Location location) const
{
    using enum meshkernel::Location;

    switch (location)
    {
    case Nodes:
        return static_cast<int>(state.m_mesh2d->GetNumNodes());
    case Edges:
        return static_cast<int>(state.m_mesh2d->GetNumEdges());
    case Faces:
        return static_cast<int>(state.m_mesh2d->GetNumFaces());
    default:
        return -1;
    }
}
