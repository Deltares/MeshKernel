#include "MeshKernel/MaxAveragingStrategy.hpp"

namespace meshkernel::averaging
{
    MaxAveragingStrategy::MaxAveragingStrategy(double const missingValue) : m_missingValue(missingValue) {}

    void MaxAveragingStrategy::Add(Point const& samplePoint, double const sampleValue)
    {
        m_result = std::max(m_result, sampleValue);
    }

    double MaxAveragingStrategy::Calculate() const
    {
        return m_result != std::numeric_limits<double>::lowest() ? m_result : m_missingValue;
    }
} // namespace meshkernel::averaging
