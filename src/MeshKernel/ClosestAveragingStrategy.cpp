#include "MeshKernel/ClosestAveragingStrategy.hpp"
#include "MeshKernel/Operations.hpp"

namespace meshkernel::averaging
{
    ClosestAveragingStrategy::ClosestAveragingStrategy(double const missingValue,
                                                       Point const& interpolationPoint,
                                                       Projection const projection) : m_result(missingValue),
                                                                                      m_interpolationPoint(interpolationPoint),
                                                                                      m_projection(projection) {}

    void ClosestAveragingStrategy::Add(Point const& samplePoint, double const sampleValue)
    {
        if (const auto squaredDistance = ComputeSquaredDistance(m_interpolationPoint, samplePoint, m_projection);
            squaredDistance < m_closestSquaredValue)
        {
            m_closestSquaredValue = squaredDistance;
            m_result = sampleValue;
        }
    }

    double ClosestAveragingStrategy::Calculate() const
    {
        return m_result;
    }
} // namespace meshkernel::averaging
