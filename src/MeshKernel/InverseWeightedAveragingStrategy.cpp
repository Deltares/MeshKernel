#include "MeshKernel/InverseWeightedAveragingStrategy.hpp"
#include "MeshKernel/Operations.hpp"

namespace meshkernel::averaging
{
    InverseWeightedAveragingStrategy::InverseWeightedAveragingStrategy(double const missingValue,
                                                                       Point const& interpolation_point,
                                                                       Projection const projection) : missingValue_(missingValue),
                                                                                                      interpolationPoint_(interpolation_point),
                                                                                                      proj_(projection) {}

    void InverseWeightedAveragingStrategy::Add(Point const& samplePoint, double const sampleValue)
    {
        double const distance = std::max(0.01, ComputeDistance(interpolationPoint_, samplePoint, proj_));
        double const weight = 1.0 / distance;
        wall_ += weight;
        result_ += weight * sampleValue;
    }

    double InverseWeightedAveragingStrategy::Calculate() const
    {
        return wall_ > 0.0 ? result_ / wall_ : missingValue_;
    }
} // namespace meshkernel::averaging
