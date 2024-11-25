#include "MeshKernelApi/PropertyCalculator.hpp"

#include <algorithm>

void meshkernelapi::OrthogonalityPropertyCalculator::Calculate(const MeshKernelState& state, const GeometryList& geometryList) const
{
    std::vector<double> values = state.m_mesh2d->GetOrthogonality();

    if (static_cast<size_t>(geometryList.num_coordinates) < values.size())
    {
        throw meshkernel::MeshKernelError("GeometryList with wrong dimensions");
    }

    std::copy(values.begin(), values.end(), geometryList.values);
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
        throw meshkernel::MeshKernelError("GeometryList with wrong dimensions");
    }

    std::copy(values.begin(), values.end(), geometryList.values);
}

int meshkernelapi::EdgeLengthPropertyCalculator::Size(const MeshKernelState& state) const
{
    return static_cast<int>(state.m_mesh2d->GetNumEdges());
}
