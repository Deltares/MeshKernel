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

#include <MeshKernel/AveragingStrategies/MaxAveragingStrategy.hpp>

namespace meshkernel::averaging
{

    double MaxAveragingStrategy::Calculate(const Point& interpolationPoint [[maybe_unused]],
                                           const std::vector<Sample>& samples) const
    {
        double result = std::numeric_limits<double>::lowest();

        for (UInt i = 0; i < samples.size(); ++i)
        {
            result = std::max(result, samples[i].value);
        }

        return result != std::numeric_limits<double>::lowest() ? result : constants::missing::doubleValue;
    }

} // namespace meshkernel::averaging
