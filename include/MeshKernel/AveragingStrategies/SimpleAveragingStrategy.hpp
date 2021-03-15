#pragma once
#include <MeshKernel/AveragingStrategies/AveragingStrategy.hpp>
#include <MeshKernel/AveragingStrategies/SimpleAveragingStrategy.hpp>
#include <MeshKernel/Entities.hpp>

namespace meshkernel::averaging
{
    /// @brief SimpleAveragingStrategy implements the averaging strategy which takes the average of the values of all added points.
    class SimpleAveragingStrategy final : public AveragingStrategy
    {
    public:
        void Add(Point const& samplePoint, double sampleValue) override;
        [[nodiscard]] double Calculate() const override;

    private:
        /// @brief The current result from which Calculate calculates the final value.
        double m_result = 0.0;

        /// @brief The number of times a value has been added to this strategy.
        size_t m_nAdds = 0;
    };
} // namespace meshkernel::averaging