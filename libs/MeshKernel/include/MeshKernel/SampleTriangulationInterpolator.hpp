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
#include <span>
#include <vector>

#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/MeshTriangulation.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/Point.hpp"
#include "MeshKernel/SampleInterpolator.hpp"
#include "MeshKernel/Utilities/RTreeFactory.hpp"

namespace meshkernel
{

    /// @brief Interpolator for sample data on a triangulated grid.
    ///
    /// The triangulation does not have to match any mesh.
    class SampleTriangulationInterpolator : public SampleInterpolator
    {
    public:
        /// @brief Constructor.
        ///
        /// The VectorType can be any array type of double precision values, e.g. std::vector, std::span.
        template <meshkernel::ValidConstDoubleArray VectorType>
        SampleTriangulationInterpolator(const VectorType& xNodes,
                                        const VectorType& yNodes,
                                        const Projection projection)
            : m_triangulation(xNodes, yNodes, projection) {}

        /// @brief Constructor.
        ///
        /// The VectorType can be any array type of double precision values, e.g. std::vector, std::span.
        template <meshkernel::ValidConstPointArray PointVector>
        SampleTriangulationInterpolator(const PointVector& nodes,
                                        const Projection projection)
            : m_triangulation(nodes, projection) {}

        /// @brief Get the number of nodes of size of the sample data.
        UInt Size() const override;

        /// @brief Set sample data
        // template <meshkernel::ValidConstDoubleArray VectorType>
        // void SetData(const int propertyId, const VectorType& sampleData) override;

        /// @brief Interpolate the sample data set at a single interpolation point.
        ///
        /// If interpolation at multiple points is required then better performance
        /// can be obtained using the Interpolate function above.
        double InterpolateValue(const int propertyId, const Point& evaluationPoint) const override;

        /// @brief Determine if the SampleInterpolator already has this sample set.
        bool Contains(const int propertyId) const override;

    private:
        /// @brief Set sample data from std::span object
        void SetDataSpan(const int propertyId, const std::span<const double>& sampleData) override;

        // DOes nothing at the moment.
        void InterpolateSpan(const int propertyId [[maybe_unused]], const Mesh2D& mesh [[maybe_unused]], const Location location [[maybe_unused]], std::span<double>& result [[maybe_unused]]) const override
        {
        }

        /// @brief Interpolate the sample data set at the interpolation nodes.
        void InterpolateSpan(const int propertyId, const std::span<const Point>& iterpolationNodes, std::span<double>& result) const override;

        /// @brief Interpolate the sample data on the element at the interpolation point.
        double InterpolateOnElement(const UInt elementId, const Point& interpolationPoint, const std::vector<double>& sampleValues) const;

        /// @brief Triangulation of the sample points
        MeshTriangulation m_triangulation;

        /// @brief Map from sample id (int) to sample data.
        std::map<int, std::vector<double>> m_sampleData;
    };

} // namespace meshkernel

inline meshkernel::UInt meshkernel::SampleTriangulationInterpolator::Size() const
{
    return m_triangulation.NumberOfNodes();
}

inline bool meshkernel::SampleTriangulationInterpolator::Contains(const int propertyId) const
{
    return m_sampleData.contains(propertyId);
}
