#include "MeshKernel/MaxAveragingStrategy.hpp"

namespace meshkernel::averaging
{
    MaxAveragingStrategy::MaxAveragingStrategy(double const missingValue) : missingValue_(missingValue) {}

    void MaxAveragingStrategy::Add(Point const& samplePoint, double const sampleValue)
    {
        result_ = std::max(result_, sampleValue);
    }

    double MaxAveragingStrategy::Calculate() const
    {
        return result_ != std::numeric_limits<double>::lowest() ? result_ : missingValue_;
    }
} // namespace meshkernel::averaging
