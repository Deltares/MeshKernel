//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2021.
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

#include <MeshKernel/AveragingStrategies/MaxAveragingStrategy.hpp>

namespace meshkernel::averaging
{

    void MaxAveragingStrategy::Reset(const Point& interpolationPoint [[maybe_unused]])
    {
        m_result = std::numeric_limits<double>::lowest();
    }

    void MaxAveragingStrategy::Add(Point const& /*samplePoint*/, double const sampleValue)
    {
        m_result = std::max(m_result, sampleValue);
    }

    double MaxAveragingStrategy::Calculate() const
    {
        return m_result != std::numeric_limits<double>::lowest() ? m_result : constants::missing::doubleValue;
    }

    double MaxAveragingStrategy::Calculate(const Point& interpolationPoint [[maybe_unused]],
                                           const std::vector<Point>& samplePoints,
                                           const std::vector<double>& sampleValues) const
    {
        double result = std::numeric_limits<double>::lowest();

        for (UInt i = 0; i < samplePoints.size (); ++i)
        {
            result = std::max(result, sampleValues[i]);
        }

        return result != std::numeric_limits<double>::lowest() ? result : constants::missing::doubleValue;
    }

} // namespace meshkernel::averaging
