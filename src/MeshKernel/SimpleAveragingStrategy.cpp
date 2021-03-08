#include "MeshKernel/SimpleAveragingStrategy.hpp"

namespace meshkernel::averaging
{
    SimpleAveragingStrategy::SimpleAveragingStrategy(double const missingValue) : missingValue_(missingValue) { }


    void SimpleAveragingStrategy::Add(Point const& samplePoint, double const sampleValue)
    {
        result_ += sampleValue;
        nAdds_ += 1;
    }

    double SimpleAveragingStrategy::Calculate() const
    {
        return nAdds_ > 0 ? result_ / static_cast<double>(nAdds_) : missingValue_;
    }
}
