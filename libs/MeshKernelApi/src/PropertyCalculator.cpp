#include "MeshKernelApi/PropertyCalculator.hpp"

#include <algorithm>

bool meshkernelapi::PropertyCalculator::IsValid(const MeshKernelState& state, const int propertyId [[maybe_unused]]) const
{
    return state.m_mesh2d != nullptr && state.m_mesh2d->GetNumNodes() > 0;
}

void meshkernelapi::OrthogonalityPropertyCalculator::Calculate(const MeshKernelState& state, const int propertyId, const GeometryList& geometryList) const
{

    if (propertyId != static_cast<int>(meshkernel::Mesh2D::Property::Orthogonality))
    {
        throw meshkernel::ConstraintError("Incorrect property id, expected = {}, actual = {}", static_cast<int>(meshkernel::Mesh2D::Property::Orthogonality), propertyId);
    }

    std::vector<double> values = state.m_mesh2d->GetOrthogonality();

    if (static_cast<size_t>(geometryList.num_coordinates) < values.size())
    {
        throw meshkernel::ConstraintError("GeometryList with wrong dimensions, {} must be greater than or equal to {}",
                                          geometryList.num_coordinates, Size(state));
    }

    std::ranges::copy(values, geometryList.values);
}

int meshkernelapi::OrthogonalityPropertyCalculator::Size(const MeshKernelState& state) const
{
    return static_cast<int>(state.m_mesh2d->GetNumEdges());
}

void meshkernelapi::EdgeLengthPropertyCalculator::Calculate(const MeshKernelState& state, const int propertyId, const GeometryList& geometryList) const
{

    if (propertyId != static_cast<int>(meshkernel::Mesh2D::Property::EdgeLength))
    {
        throw meshkernel::ConstraintError("Incorrect property id, expected = {}, actual = {}", static_cast<int>(meshkernel::Mesh2D::Property::EdgeLength), propertyId);
    }

    state.m_mesh2d->ComputeEdgesLengths();
    std::vector<double> values = state.m_mesh2d->m_edgeLengths;

    if (static_cast<size_t>(geometryList.num_coordinates) < values.size())
    {
        throw meshkernel::ConstraintError("GeometryList with wrong dimensions, {} must be greater than or equal to {}",
                                          geometryList.num_coordinates, Size(state));
    }

    std::ranges::copy(values, geometryList.values);
}

int meshkernelapi::EdgeLengthPropertyCalculator::Size(const MeshKernelState& state) const
{
    return static_cast<int>(state.m_mesh2d->GetNumEdges());
}

bool meshkernelapi::DepthSamplePropertyCalculator::IsValid(const MeshKernelState& state, const int propertyId) const
{
    return state.m_mesh2d != nullptr && state.m_mesh2d->GetNumNodes() > 0 && state.m_sampleInterpolator != nullptr && state.m_sampleInterpolator->Contains(propertyId);
}

void meshkernelapi::DepthSamplePropertyCalculator::Calculate(const MeshKernelState& state, const int propertyId, const GeometryList& geometryList) const
{

    if (geometryList.num_coordinates < Size(state))
    {
        throw meshkernel::ConstraintError("GeometryList with wrong dimensions, {} must be greater than or equal to {}",
                                          geometryList.num_coordinates, Size(state));
    }

    std::span<double> interpolatedSampleData(geometryList.values, geometryList.num_coordinates);
    state.m_sampleInterpolator->Interpolate(propertyId, state.m_mesh2d->Nodes(), interpolatedSampleData);
}

int meshkernelapi::DepthSamplePropertyCalculator::Size(const MeshKernelState& state) const
{
    return static_cast<int>(state.m_mesh2d->GetNumNodes());
}
