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

#include <MeshKernel/AveragingStrategies/InverseWeightedAveragingStrategy.hpp>
#include <MeshKernel/Operations.hpp>

namespace meshkernel::averaging
{
    InverseWeightedAveragingStrategy::InverseWeightedAveragingStrategy(size_t minNumSamples,
                                                                       Projection const projection) : m_minNumSamples(minNumSamples),
                                                                                                      m_projection(projection) {}

    double InverseWeightedAveragingStrategy::Calculate(const Point& interpolationPoint,
                                                       const std::vector<Point>& samplePoints,
                                                       const std::vector<double>& sampleValues) const
    {
        double result = 0.0;
        double wall = 0.0;

        for (UInt i = 0; i < samplePoints.size(); ++i)
        {
            double weight = 1.0 / std::max(0.01, ComputeDistance(interpolationPoint, samplePoints[i], m_projection));
            wall += weight;
            result += weight * sampleValues[i];
        }

        return wall >= static_cast<double>(m_minNumSamples) ? result / wall : constants::missing::doubleValue;
    }

} // namespace meshkernel::averaging
