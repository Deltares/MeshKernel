#include "MeshKernel/ClosestAveragingStrategy.hpp"
#include "MeshKernel/Operations.hpp"

namespace meshkernel::averaging
{
    ClosestAveragingStrategy::ClosestAveragingStrategy(double const missing_value,
                                                       Point const& interpolation_point,
                                                       Projection const projection) : result_(missing_value),
                                                                                      interpolationPoint_(interpolation_point),
                                                                                      proj_(projection) {}

    void ClosestAveragingStrategy::Add(Point const& samplePoint, double const sampleValue)
    {
        if (const auto squaredDistance = ComputeSquaredDistance(interpolationPoint_, samplePoint, proj_);
            squaredDistance < closestSquaredValue_)
        {
            closestSquaredValue_ = squaredDistance;
            result_ = sampleValue;
        }
    }

    double ClosestAveragingStrategy::Calculate() const
    {
        return result_;
    }
} // namespace meshkernel::averaging
