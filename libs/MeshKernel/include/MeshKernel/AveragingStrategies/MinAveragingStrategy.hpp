//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2023.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// contact: delft3d.support@deltares.nl
// Stichting Deltares
// P.O. Box 177
// 2600 MH Delft, The Netherlands
//
// All indications and logos of, and references to, "Delft3D" and "Deltares"
// are registered trademarks of Stichting Deltares, and remain the property of
// Stichting Deltares. All rights reserved.
//
//------------------------------------------------------------------------------

#pragma once
#include <MeshKernel/AveragingStrategies/AveragingStrategy.hpp>
#include <MeshKernel/Entities.hpp>

namespace meshkernel::averaging
{
    /// @brief MinAbsAveragingStrategy implements the AveragingStrategy which takes the minimum value of the added points.
    class MinAveragingStrategy final : public AveragingStrategy
    {
    public:
        /// @brief Reset the state of the minimum value averaging-strategy.
        void Reset(const Point& interpolationPoint) override;

        void Add(Point const& samplePoint, double sampleValue) override;

        [[nodiscard]] double Calculate() const override;

        /// @brief Calculates the average value based on the sample values.
        /// @param[in] interpolationPoint The point for which the average should be calculated.
        /// @param[in] samplePoints The sample points to used by this strategy.
        /// @param[in] sampleValues The sample values  associated with each sample point.
        /// @return The calculated average
        double Calculate(const Point& interpolationPoint,
                         const std::vector<Point>& samplePoints,
                         const std::vector<double>& sampleValues) const override;

    private:
        /// @brief The current result returned in Calculate
        double m_result = std::numeric_limits<double>::max();
    };
} // namespace meshkernel::averaging
