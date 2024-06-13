#include "MeshKernel/UndoActions/SphericalCoordinatesOffsetAction.hpp"
#include "MeshKernel/Mesh2D.hpp"

std::unique_ptr<meshkernel::SphericalCoordinatesOffsetAction>
meshkernel::SphericalCoordinatesOffsetAction::Create(Mesh2D& mesh, const double minx, const double maxx)
{
    return std::make_unique<meshkernel::SphericalCoordinatesOffsetAction>(mesh, minx, maxx);
}

meshkernel::SphericalCoordinatesOffsetAction::SphericalCoordinatesOffsetAction(Mesh2D& mesh, const double minx, const double maxx)
    : BaseMeshUndoAction<SphericalCoordinatesOffsetAction, Mesh2D>(mesh), m_xMin(minx), m_xMax(maxx) {}

double meshkernel::SphericalCoordinatesOffsetAction::MinX() const
{
    return m_xMin;
}

double meshkernel::SphericalCoordinatesOffsetAction::MaxX() const
{
    return m_xMax;
}

void meshkernel::SphericalCoordinatesOffsetAction::AddDecrease(const UInt nodeId)
{
    m_offsetNodesDecrease.push_back(nodeId);
}

void meshkernel::SphericalCoordinatesOffsetAction::AddIncrease(const UInt nodeId)
{
    m_offsetNodesIncrease.push_back(nodeId);
}

void meshkernel::SphericalCoordinatesOffsetAction::ApplyOffset(std::vector<Point>& nodes) const
{
    for (UInt index : m_offsetNodesDecrease)
    {
        nodes[index].x -= 360.0;
    }

    for (UInt index : m_offsetNodesIncrease)
    {
        nodes[index].x += 360.0;
    }
}

void meshkernel::SphericalCoordinatesOffsetAction::UndoOffset(std::vector<Point>& nodes) const
{
    for (UInt index : m_offsetNodesDecrease)
    {
        nodes[index].x += 360.0;
    }

    for (UInt index : m_offsetNodesIncrease)
    {
        nodes[index].x -= 360.0;
    }
}

std::uint64_t meshkernel::SphericalCoordinatesOffsetAction::MemorySize() const
{
    std::uint64_t size = sizeof(*this) + static_cast<std::uint64_t>(m_offsetNodesDecrease.capacity() + m_offsetNodesIncrease.capacity()) * sizeof(UInt);
    return size;
}
