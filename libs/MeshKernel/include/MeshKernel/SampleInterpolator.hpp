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
#include <string>
#include <vector>

#include "MeshKernel/AveragingInterpolation.hpp" //  Only for the enum Method. should move to Definitions.hpp
#include "MeshKernel/AveragingStrategies/AveragingStrategy.hpp"
#include "MeshKernel/AveragingStrategies/AveragingStrategyFactory.hpp"
#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/MeshTriangulation.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/Point.hpp"
#include "MeshKernel/Utilities/RTreeFactory.hpp"

namespace meshkernel
{

    /// @brief Interface for sample interpolation
    class SampleInterpolator
    {
    public:
        virtual ~SampleInterpolator() = default;

        /// @brief Get the number of sample points
        virtual UInt Size() const = 0;

        /// @brief Set sample data
        template <meshkernel::ValidConstDoubleArray VectorType>
        void SetData(const int propertyId, const VectorType& sampleData)
        {
            const std::span<const double> spanSampleData(sampleData.data(), sampleData.size());
            SetDataSpan(propertyId, sampleData);
        }

        /// @brief Interpolate the sample data set at the interpolation nodes.
        template <meshkernel::ValidConstPointArray PointVectorType, meshkernel::ValidConstDoubleArray ScalarVectorType>
        void Interpolate(const int propertyId, const PointVectorType& iterpolationNodes, ScalarVectorType& result) const
        {
            const std::span<const Point> spanInterpolationNodes(iterpolationNodes.data(), iterpolationNodes.size());
            std::span<double> spanResult(result.data(), result.size());
            InterpolateSpan(propertyId, spanInterpolationNodes, spanResult);
        }

        /// @brief Interpolate the sample data set at the interpolation point from a mesh.
        template <meshkernel::ValidConstDoubleArray ScalarVectorType>
        void Interpolate(const int propertyId, const Mesh2D& mesh, const Location location, ScalarVectorType& result) const
        {
            std::span<double> spanResult(result.data(), result.size());
            InterpolateSpan(propertyId, mesh, location, spanResult);
        }

        /// @brief Interpolate the sample data set at a single interpolation point.
        ///
        /// If interpolation at multiple points is required then better performance
        /// can be obtained using the Interpolate function above.
        virtual double InterpolateValue(const int propertyId, const Point& evaluationPoint) const = 0;

        /// @brief Determine if the SampleInterpolator already has this sample set.
        virtual bool Contains(const int propertyId) const = 0;

    private:
        /// @brief Set sample data contained in an std::span object
        virtual void SetDataSpan(const int propertyId, const std::span<const double>& sampleData) = 0;

        /// @brief Interpolate the sample data set at the interpolation nodes.
        virtual void InterpolateSpan(const int propertyId, const std::span<const Point>& iterpolationNodes, std::span<double>& result) const = 0;

        /// @brief Interpolate the sample data.
        // TODO may need a polygon too: const Polygons& polygon
        virtual void InterpolateSpan(const int propertyId, const Mesh2D& mesh, const Location location, std::span<double>& result) const = 0;

        // TODO the map -> sample-data can be here
        // with protected members to access the data.
    };

} // namespace meshkernel
