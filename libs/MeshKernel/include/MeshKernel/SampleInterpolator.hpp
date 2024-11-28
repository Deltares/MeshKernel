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

    /// @brief Interpolator for sample data on a triangulated grid.
    ///
    /// The triangulation does not have to match any mesh.
    class SampleInterpolator
    {
    public:
        /// @brief Constructor.
        ///
        /// The VectorType can be any array type of double precision values, e.g. std::vector, std::span.
        template <ValidConstDoubleArray VectorType>
        SampleInterpolator(const VectorType& xNodes,
                           const VectorType& yNodes,
                           const Projection projection)
            : m_triangulation(xNodes, yNodes, projection)
        {
        }

        /// @brief Get the number of nodes of size of the sample data.
        UInt Size() const;

        /// @brief Set sample data
        template <ValidConstDoubleArray VectorType>
        void SetData(const int sampleId, const VectorType& sampleData);

        /// @brief Interpolate the sample data set at the interpolation nodes.
        template <ValidConstPointArray PointVectorType, ValidConstDoubleArray ScalarVectorType>
        void Interpolate(const int sampleId, const PointVectorType& iterpolationNodes, ScalarVectorType& result) const;

        /// @brief Determine if the SampleInterpolator already has this sample set.
        bool Contains(const int sampleId) const;

    private:
        /// @brief Interpolate the sample data on the element at the interpolation point.
        double InterpolateOnElement(const UInt elementId, const Point& interpolationPoint, const std::vector<double>& sampleValues) const;

        /// @brief Triangulation of the sample points
        MeshTriangulation m_triangulation;

        /// @brief Map from sample id (int) to sample data.
        std::map<int, std::vector<double>> m_sampleData;
    };

} // namespace meshkernel

inline meshkernel::UInt meshkernel::SampleInterpolator::Size() const
{
    return m_triangulation.NumberOfNodes();
}

inline bool meshkernel::SampleInterpolator::Contains(const int sampleId) const
{
    return m_sampleData.contains(sampleId);
}

template <meshkernel::ValidConstDoubleArray VectorType>
void meshkernel::SampleInterpolator::SetData(const int sampleId, const VectorType& sampleData)
{
    if (m_triangulation.NumberOfNodes() != sampleData.size())
    {
        throw ConstraintError("The sample data array does not have the same number of elements as the number of nodes in the triangulation: {} /= {}",
                              m_triangulation.NumberOfNodes(), sampleData.size());
    }

    m_sampleData[sampleId].assign(sampleData.begin(), sampleData.end());
}

template <meshkernel::ValidConstPointArray PointVectorType, meshkernel::ValidConstDoubleArray ScalarVectorType>
void meshkernel::SampleInterpolator::Interpolate(const int sampleId, const PointVectorType& iterpolationNodes, ScalarVectorType& result) const
{
    if (!Contains(sampleId))
    {
        throw ConstraintError("Sample interpolator does not contain the id: {}.", sampleId);
    }

    if (iterpolationNodes.size() != result.size())
    {
        throw ConstraintError("The arrays for interpolation nodes and the results are different sizes: {} /= {}",
                              iterpolationNodes.size(), result.size());
    }

    const std::vector<double>& propertyValues = m_sampleData.at(sampleId);

    for (size_t i = 0; i < iterpolationNodes.size(); ++i)
    {
        result[i] = constants::missing::doubleValue;

        if (!iterpolationNodes[i].IsValid())
        {
            continue;
        }

        const UInt elementId = m_triangulation.FindNearestFace(iterpolationNodes[i]);

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
