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

#pragma once

#include <map>
#include <string>
#include <vector>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/MeshTriangulation.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/Point.hpp"
#include "MeshKernel/Utilities/RTreeFactory.hpp"

namespace meshkernel
{

    class SampleInterpolator
    {
    public:
        template <class VectorType>
        SampleInterpolator(const VectorType& xNodes, const VectorType& yNodes) : m_triangulation(xNodes, yNodes) {}

        /// @brief Set sample data
        template <class VectorType>
        void SetData(const std::string& name, const VectorType& sampleData);

        /// @brief Interpolate the sample data set at the interpolation nodes.
        template <class PointVectorType, class ScalarVectorType>
        void Interpolate(const std::string& name, const PointVectorType& iterpolationNodes, ScalarVectorType& result) const;

        /// @brief Determine if the SampleInterpolator already has this sample set.
        bool Contains(const std::string& name) const;

    private:
        /// @brief Interpolate the sample data on the element at the interpolation point.
        double InterpolateOnElement(const UInt elementId, const Point& interpolationPoint, const std::vector<double>& sampleValues) const;

        /// @brief Triangulation of the sample points
        MeshTriangulation m_triangulation;

        /// @brief Map from sample name to sample data.
        std::map<std::string, std::vector<double>> m_sampleData;
    };

} // namespace meshkernel

inline bool meshkernel::SampleInterpolator::Contains(const std::string& name) const
{
    return m_sampleData.contains(name);
}

template <class VectorType>
void meshkernel::SampleInterpolator::SetData(const std::string& name, const VectorType& sampleData)
{
    if (m_triangulation.NumberOfNodes() != sampleData.size())
    {
        throw ConstraintError("The sample data array does not have the same number of elements as the number of nodes in the triangulation: {} /= {}",
                              m_triangulation.NumberOfNodes(), sampleData.size());
    }

    if (!Contains(name))
    {
        m_sampleData.emplace(name, std::vector<double>());
    }

    m_sampleData[name].assign(sampleData.begin(), sampleData.end());
}

template <class PointVectorType, class ScalarVectorType>
void meshkernel::SampleInterpolator::Interpolate(const std::string& name, const PointVectorType& iterpolationNodes, ScalarVectorType& result) const
{
    if (!Contains(name))
    {
        throw ConstraintError("Sample interpolator does not contain the name: {}.", name);
    }

    if (iterpolationNodes.size() != result.size())
    {
        throw ConstraintError("The arrays for interpolation nodes and the results are different sizes: {} /= {}",
                              iterpolationNodes.size(), result.size());
    }

    const std::vector<double>& propertyValues = m_sampleData.at(name);

    for (size_t i = 0; i < iterpolationNodes.size(); ++i)
    {
        result[i] = constants::missing::doubleValue;

        if (!iterpolationNodes[i].IsValid())
        {
            continue;
        }

        const UInt elementId = m_triangulation.FindFace(iterpolationNodes[i]);

        if (elementId == constants::missing::uintValue)
        {
            continue;
        }

        if (m_triangulation.PointIsInElement(iterpolationNodes[i], elementId))
        {
            result[i] = InterpolateOnElement(elementId, iterpolationNodes[i], propertyValues);
        }
    }
}
