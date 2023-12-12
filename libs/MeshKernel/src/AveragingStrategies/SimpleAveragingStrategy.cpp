﻿//---- GPL ---------------------------------------------------------------------
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

#include <algorithm>

#include <MeshKernel/AveragingStrategies/SimpleAveragingStrategy.hpp>

namespace meshkernel::averaging
{

    SimpleAveragingStrategy::SimpleAveragingStrategy(size_t minNumSamples) : m_minNumPoints(minNumSamples) {}

    double SimpleAveragingStrategy::Calculate(const Point& interpolationPoint [[maybe_unused]],
                                              const std::vector<Sample>& samples) const
    {
        auto sumSampleValue = [](double value, const Sample& sample)
        { return value + sample.value; };
        double result = std::accumulate(samples.begin(), samples.end(), 0.0, sumSampleValue);
        return samples.size() >= m_minNumPoints ? result / static_cast<double>(samples.size()) : constants::missing::doubleValue;
    }

} // namespace meshkernel::averaging
