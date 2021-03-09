#include "MeshKernel/InverseWeightedAveragingStrategy.hpp"
#include "MeshKernel/Operations.hpp"

namespace meshkernel::averaging
{
    InverseWeightedAveragingStrategy::InverseWeightedAveragingStrategy(double const missingValue,
                                                                       Point const& interpolation_point,
                                                                       Projection const projection) : m_missingValue(missingValue),
                                                                                                      m_interpolationPoint(interpolation_point),
                                                                                                      m_projection(projection) {}

    void InverseWeightedAveragingStrategy::Add(Point const& samplePoint, double const sampleValue)
    {
        double const distance = std::max(0.01, ComputeDistance(m_interpolationPoint, samplePoint, m_projection));
        double const weight = 1.0 / distance;
        m_wall += weight;
        m_result += weight * sampleValue;
    }

    double InverseWeightedAveragingStrategy::Calculate() const
    {
        return m_wall > 0.0 ? m_result / m_wall : m_missingValue;
    }
} // namespace meshkernel::averaging
