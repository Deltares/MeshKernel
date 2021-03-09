#include "MeshKernel/MinAbsAveragingStrategy.hpp"

namespace meshkernel::averaging
{
    MinAbsAveragingStrategy::MinAbsAveragingStrategy(double const missingValue) : m_missingValue(missingValue) {}

    void MinAbsAveragingStrategy::Add(Point const& samplePoint, double const sampleValue)
    {
        m_result = std::min(m_result, std::abs(sampleValue));
    }

    double MinAbsAveragingStrategy::Calculate() const
    {
        return m_result != std::numeric_limits<double>::max() ? m_result : m_missingValue;
    }
} // namespace meshkernel::averaging
