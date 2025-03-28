//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2024.
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

#include "MeshKernel/SampleInterpolator.hpp"
#include "MeshKernel/Exceptions.hpp"

void meshkernel::SampleInterpolator::SetData(const int propertyId, const std::span<const double> sampleData)
{
    if (Size() != sampleData.size())
    {
        throw ConstraintError("The sample data array does not have the same number of values as the sample point set: {} /= {}",
                              Size(), sampleData.size());
    }

    m_sampleData[propertyId].assign(sampleData.begin(), sampleData.end());
}

const std::vector<double>& meshkernel::SampleInterpolator::GetSampleData(const int propertyId) const
{

    if (!Contains(propertyId))
    {
        throw ConstraintError("The sample data set for property id {}, has not been defined", propertyId);
    }

    return m_sampleData.at(propertyId);
}
