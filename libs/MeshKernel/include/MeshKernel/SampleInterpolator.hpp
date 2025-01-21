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

#include "MeshKernel/AveragingInterpolation.hpp"
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

        /// @brief Destructor
        virtual ~SampleInterpolator() = default;

        /// @brief Set sample data
        void SetData(const int propertyId, const std::span<const double> sampleData);

        /// @brief Get the number of sample points
        virtual UInt Size() const = 0;

        /// @brief Determine if the SampleInterpolator already has this sample set.
        bool Contains(const int propertyId) const;

        /// @brief Interpolate the sample data set at the interpolation nodes.
        virtual void Interpolate(const int propertyId, const std::span<const Point> iterpolationNodes, std::span<double> result) const = 0;

        /// @brief Interpolate the sample data.
        virtual void Interpolate(const int propertyId, const Mesh2D& mesh, const Location location, std::span<double> result) const = 0;

        /// @brief Interpolate the sample data set at a single interpolation point.
        ///
        /// If interpolation at multiple points is required then better performance
        /// can be obtained using the Interpolate function above.
        virtual double InterpolateValue(const int propertyId, const Point& evaluationPoint) const = 0;

    protected:
        /// @brief Get the sample property data for the id.
        const std::vector<double>& GetSampleData(const int propertyId) const;

    private:
        /// @brief Map from sample id (int) to sample data.
        std::map<int, std::vector<double>> m_sampleData;
    };

} // namespace meshkernel

inline bool meshkernel::SampleInterpolator::Contains(const int propertyId) const
{
    return m_sampleData.contains(propertyId);
}
