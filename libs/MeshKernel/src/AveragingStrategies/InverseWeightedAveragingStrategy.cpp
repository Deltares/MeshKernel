#include <MeshKernel/AveragingStrategies/InverseWeightedAveragingStrategy.hpp>
#include <MeshKernel/Operations.hpp>

namespace meshkernel::averaging
{
    InverseWeightedAveragingStrategy::InverseWeightedAveragingStrategy(size_t minNumSamples,
                                                                       Projection const projection) : m_minNumSamples(minNumSamples),
                                                                                                      m_projection(projection) {}

    void InverseWeightedAveragingStrategy::Reset(const Point& interpolationPoint)
    {
        m_interpolationPoint = interpolationPoint;
        m_result = 0.0;
        m_wall = 0.0;
    }

    void InverseWeightedAveragingStrategy::Add(Point const& samplePoint, double const sampleValue)
    {
        double const distance = std::max(0.01, ComputeDistance(m_interpolationPoint, samplePoint, m_projection));
        double const weight = 1.0 / distance;
        m_wall += weight;
        m_result += weight * sampleValue;
    }

    double InverseWeightedAveragingStrategy::Calculate() const
    {
        return m_wall >= m_minNumSamples ? m_result / m_wall : constants::missing::doubleValue;
    }
} // namespace meshkernel::averaging
