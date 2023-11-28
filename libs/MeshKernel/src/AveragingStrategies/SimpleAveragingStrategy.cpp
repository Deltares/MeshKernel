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

#include <MeshKernel/AveragingStrategies/SimpleAveragingStrategy.hpp>

namespace meshkernel::averaging
{

    SimpleAveragingStrategy::SimpleAveragingStrategy(size_t minNumSamples) : m_minNumPoints(minNumSamples) {}

    void SimpleAveragingStrategy::Reset(const Point& interpolationPoint [[maybe_unused]])
    {
        m_result = 0.0;
        m_nAdds = 0;
    }

    void SimpleAveragingStrategy::Add(Point const& /*samplePoint*/, double const sampleValue)
    {
        m_result += sampleValue;
        m_nAdds += 1;
    }

    double SimpleAveragingStrategy::Calculate() const
    {
        return m_nAdds >= m_minNumPoints ? m_result / static_cast<double>(m_nAdds) : constants::missing::doubleValue;
    }

    double SimpleAveragingStrategy::Calculate(const Point& interpolationPoint [[maybe_unused]],
                                              const std::vector<Point>& samplePoints [[maybe_unused]],
                                              const std::vector<double>& sampleValues) const
    {
        double result = 0.0;

        for (UInt i = 0; i < sampleValues.size (); ++i)
        {
            result += sampleValues[i];
        }

        return sampleValues.size () >= m_minNumPoints ? result / static_cast<double>(sampleValues.size ()) : constants::missing::doubleValue;
    }

} // namespace meshkernel::averaging
