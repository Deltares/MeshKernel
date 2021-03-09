#include "MeshKernel/MinAveragingStrategy.hpp"

namespace meshkernel::averaging
{
    MinAveragingStrategy::MinAveragingStrategy(double const missingValue) : m_missingValue(missingValue) {}

    void MinAveragingStrategy::Add(Point const& samplePoint, double const sampleValue)
    {
        m_result = std::min(m_result, sampleValue);
    }

    double MinAveragingStrategy::Calculate() const
    {
        return m_result != std::numeric_limits<double>::max() ? m_result : m_missingValue;
    }
} // namespace meshkernel::averaging
