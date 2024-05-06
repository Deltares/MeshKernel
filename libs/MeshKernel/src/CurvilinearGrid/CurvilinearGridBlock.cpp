#include "MeshKernel/CurvilinearGrid/CurvilinearGridBlock.hpp"
#include "MeshKernel/Constants.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"

meshkernel::CurvilinearGridBlock::CurvilinearGridBlock(const CurvilinearGridNodeIndices& bottomLeft, const CurvilinearGridNodeIndices& topRight)
    : m_bottomLeft(bottomLeft), m_topRight(topRight)
{
    const UInt rows = m_topRight.m_n - m_bottomLeft.m_n;
    const UInt cols = m_topRight.m_m - m_bottomLeft.m_m;

    lin_alg::ResizeAndFillMatrix(m_gridNodes, rows, cols, false, {constants::missing::doubleValue, constants::missing::doubleValue});
}

void meshkernel::CurvilinearGridBlock::CopyFrom(const CurvilinearGrid& grid)
{
    const UInt rows = m_topRight.m_n - m_bottomLeft.m_n;
    const UInt cols = m_topRight.m_m - m_bottomLeft.m_m;

    for (UInt r = 0; r < rows; ++r)
    {
        for (UInt c = 0; c < cols; ++c)
        {
            m_gridNodes(r, c) = grid.GetNode(r + m_bottomLeft.m_n, c + m_bottomLeft.m_m);
        }
    }
}

void meshkernel::CurvilinearGridBlock::Swap(CurvilinearGrid& grid)
{
    const UInt rows = m_topRight.m_n - m_bottomLeft.m_n;
    const UInt cols = m_topRight.m_m - m_bottomLeft.m_m;

    for (UInt r = 0; r < rows; ++r)
    {
        for (UInt c = 0; c < cols; ++c)
        {
            std::swap(m_gridNodes(r, c), grid.GetNode(r + m_bottomLeft.m_n, c + m_bottomLeft.m_m));
        }
    }
}

std::uint64_t meshkernel::CurvilinearGridBlock::MemorySize() const
{
    std::uint64_t result = 0;

    result += sizeof(*this);
    result += static_cast<std::uint64_t>(m_gridNodes.rows() * m_gridNodes.cols() * sizeof(Point));

    return result;
}
