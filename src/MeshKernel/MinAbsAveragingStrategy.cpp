#include "MeshKernel/MinAbsAveragingStrategy.hpp"

namespace meshkernel::averaging
{
	MinAbsAveragingStrategy::MinAbsAveragingStrategy(double const missingValue) : missingValue_(missingValue) { }

    void MinAbsAveragingStrategy::Add(Point const& samplePoint, double const sampleValue)
    {
	    result_ = std::min(result_, std::abs(sampleValue));
    }

    double MinAbsAveragingStrategy::Calculate() const
    {
        return result_ != std::numeric_limits<double>::max() ? result_ : missingValue_;
    }
}
