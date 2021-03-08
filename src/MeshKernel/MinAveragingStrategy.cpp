#include "MeshKernel/MinAveragingStrategy.hpp"

namespace meshkernel::averaging
{
	MinAveragingStrategy::MinAveragingStrategy(double const missingValue) : missingValue_(missingValue) { }

    void MinAveragingStrategy::Add(Point const& samplePoint, double const sampleValue)
    {
	    result_ = std::min(result_, sampleValue);
    }

    double MinAveragingStrategy::Calculate() const
    {
        return result_ != std::numeric_limits<double>::max() ? result_ : missingValue_;
    }
}
