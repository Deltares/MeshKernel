#pragma once
#include "AveragingStrategy.hpp"

namespace meshkernel::averaging
{
    /// @brief SimpleAveragingStrategy implements the averaging strategy which takes the average of the values of all added points.
    class SimpleAveragingStrategy final : public AveragingStrategy
    {
    public:        
        /// @brief Construct a new SimpleAveragingStrategy implements the simple averaging strategy.
        explicit SimpleAveragingStrategy(double missingValue);

        void Add(Point const& samplePoint, double sampleValue) override;
        [[nodiscard]] double Calculate() const override;
    private:
        double result_ = 0.0;
        size_t nAdds_ = 0;

        double const missingValue_;
    };
} // namespace meshkernel