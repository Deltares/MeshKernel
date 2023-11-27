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

#include <MeshKernel/AveragingStrategies/ClosestAveragingStrategy.hpp>
#include <MeshKernel/Operations.hpp>

namespace meshkernel::averaging
{
    ClosestAveragingStrategy::ClosestAveragingStrategy(Projection const projection) : m_projection(projection) {}

    void ClosestAveragingStrategy::Reset(const Point& interpolationPoint)
    {
        m_interpolationPoint = interpolationPoint;
        m_result = constants::missing::doubleValue;
        m_closestSquaredValue = std::numeric_limits<double>::max();
    }

    void ClosestAveragingStrategy::Add(Point const& samplePoint, double const sampleValue)
    {
        if (const auto squaredDistance = ComputeSquaredDistance(m_interpolationPoint, samplePoint, m_projection);
            squaredDistance < m_closestSquaredValue)
        {
            m_closestSquaredValue = squaredDistance;
            m_result = sampleValue;
        }
    }

    double ClosestAveragingStrategy::Calculate() const
    {
        return m_result;
    }
} // namespace meshkernel::averaging
