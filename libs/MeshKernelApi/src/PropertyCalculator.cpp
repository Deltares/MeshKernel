#include "MeshKernelApi/PropertyCalculator.hpp"

#include <algorithm>

bool meshkernelapi::PropertyCalculator::IsValid(const MeshKernelState& state [[maybe_unused]]) const
{
    return true;
}

void meshkernelapi::OrthogonalityPropertyCalculator::Calculate(const MeshKernelState& state, const GeometryList& geometryList) const
{
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

void meshkernelapi::EdgeLengthPropertyCalculator::Calculate(const MeshKernelState& state, const GeometryList& geometryList) const
{
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

const std::string meshkernelapi::DepthSamplePropertyCalculator::SampleName = "depth";

bool meshkernelapi::DepthSamplePropertyCalculator::IsValid(const MeshKernelState& state) const
{
    return state.m_sampleInterpolator != nullptr;
}

void meshkernelapi::DepthSamplePropertyCalculator::Calculate(const MeshKernelState& state, const GeometryList& geometryList) const
{

    if (geometryList.num_coordinates < Size(state))
    {
        throw meshkernel::ConstraintError("GeometryList with wrong dimensions, {} must be greater than or equal to {}",
                                          geometryList.num_coordinates, Size(state));
    }

    std::span<double> interpolatedSampleData(geometryList.values, geometryList.num_coordinates);
    state.m_sampleInterpolator->Interpolate(SampleName, state.m_mesh2d->Nodes(), interpolatedSampleData);
}

int meshkernelapi::DepthSamplePropertyCalculator::Size(const MeshKernelState& state) const
{
    return static_cast<int>(state.m_mesh2d->GetNumNodes());
}
