#include "MeshKernel/SimpleAveragingStrategy.hpp"

namespace meshkernel::averaging
{
    void SimpleAveragingStrategy::Add(Point const& samplePoint, double const sampleValue)
    {
        m_result += sampleValue;
        m_nAdds += 1;
    }

    double SimpleAveragingStrategy::Calculate() const
    {
        return m_nAdds > 0 ? m_result / static_cast<double>(m_nAdds) : doubleMissingValue;
    }
} // namespace meshkernel::averaging
