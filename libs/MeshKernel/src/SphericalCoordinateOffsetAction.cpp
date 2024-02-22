#include "MeshKernel/SphericalCoordinateOffsetAction.hpp"

/// \brief Compute the approximate amount of memory being used, in bytes.
std::uint64_t meshkernel::SphericalCoordinateOffsetAction::MemorySize() const
{
    std::uint64_t size = sizeof(*this) + m_offsetNodes.size() * sizeof(UInt);
    return size;
}

/// @brief Print the add node action to the stream
void meshkernel::SphericalCoordinateOffsetAction::Print(std::ostream& out) const
{
    out << "SphericalCoordinateOffsetActionPrint: state " << to_string(State())
        << ", x-minimum = " << m_xMin
        << ", x-maximum = " << m_xMax
        << std::endl;

    out << "  nodes with offset:";

    for (UInt node : m_offsetNodes)
    {
        out << " " << node;
    }

    out << std::endl;
}
