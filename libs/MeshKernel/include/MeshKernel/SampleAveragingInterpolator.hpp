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
#include "MeshKernel/SampleInterpolator.hpp"
#include "MeshKernel/Utilities/RTreeFactory.hpp"

namespace meshkernel
{

    /// @brief Parameters used by the averaging interpolation
    struct InterpolationParameters
    {
        /// @brief Which averaging method should be used
        AveragingInterpolation::Method m_method = AveragingInterpolation::Method::SimpleAveraging;

        /// @brief The relative search radius
        double m_relativeSearchRadius = 1.0;

        /// @brief If no point is found in polygon then just used the closest point
        bool m_useClosestIfNoneFound = true;

        /// @brief The minimum number of samples for several averaging methods.
        UInt m_minimumNumberOfSamples = 10;
    };

    /// @brief Interpolator for sample data using an averaging scheme
    class SampleAveragingInterpolator : public SampleInterpolator
    {
    public:
        /// @brief Constructor.
        SampleAveragingInterpolator(const std::span<const double> xNodes,
                                    const std::span<const double> yNodes,
                                    const Projection projection,
                                    const InterpolationParameters& interpolationParameters);

        /// @brief Constructor.
        SampleAveragingInterpolator(const std::span<const Point> nodes,
                                    const Projection projection,
                                    const InterpolationParameters& interpolationParameters);

        /// @brief Get the number of nodes of size of the sample data.
        UInt Size() const override;

        /// @brief Interpolate the sample data set at the interpolation nodes.
        void Interpolate(const int propertyId, const std::span<const Point> iterpolationNodes, std::span<double> result) const override;

        /// @brief Interpolate the sample data set at the locationd defined.
        void Interpolate(const int propertyId, const Mesh2D& mesh, const Location location, std::span<double> result) const override;

        /// @brief Interpolate the sample data set at a single interpolation point.
        ///
        /// If interpolation at multiple points is required then better performance
        /// can be obtained using the Interpolate function above.
        double InterpolateValue(const int propertyId, const Point& evaluationPoint) const override;

    private:
        static constexpr UInt MaximumNumberOfEdgesPerNode = 16; ///< Maximum number of edges per node

        /// @brief Combine x- and y-coordinate arrays to a point array
        static std::vector<Point> CombineCoordinates(const std::span<const double> xNodes, const std::span<const double> yNodes);

        /// @brief Interpolate the sample data on the element at the interpolation point.
        double InterpolateOnElement(const UInt elementId, const Point& interpolationPoint, const std::vector<double>& sampleValues) const;

        /// @brief Interpolate at the points from points found within a polygon
        double ComputeOnPolygon(const int propertyId,
                                std::vector<Point>& polygon,
                                const Point& interpolationPoint,
                                const Projection projection,
                                std::vector<Sample>& sampleCache) const;

        /// @brief Compute the search radius
        double GetSearchRadiusSquared(const std::vector<Point>& searchPolygon,
                                      const Point& interpolationPoint,
                                      const Projection projection) const;

        /// @brief Compute the search polygon
        void GenerateSearchPolygon(const double relativeSearchRadius,
                                   const Point& interpolationPoint,
                                   std::vector<Point>& polygon,
                                   const Projection projection) const;

        /// @brief Get the sample value at the index from the rtree result cache
        double GetSampleValueFromRTree(const int propertyId, const UInt index) const;

        /// @brief Compute the average from the neighbours.
        double ComputeInterpolationResultFromNeighbors(const int propertyId,
                                                       const Point& interpolationPoint,
                                                       const std::vector<Point>& searchPolygon,
                                                       const Projection projection,
                                                       std::vector<Sample>& sampleCache) const;

        /// @brief Interpolate at the mesh nodes
        void InterpolateAtNodes(const int propertyId, const Mesh2D& mesh, std::span<double>& result) const;

        /// @brief Interpolate at edge centres by averaging the node values
        void InterpolateAtEdgeCentres(const Mesh2D& mesh,
                                      const std::span<double>& nodeResult,
                                      std::span<double>& result) const;

        /// @brief Interpolate at the face centres
        void InterpolateAtFaces(const int propertyId, const Mesh2D& mesh, std::span<double>& result) const;

        /// @brief The sample points
        std::vector<Point> m_samplePoints;

        // Should use the m_projection from the mesh in the interpolate function or
        // needs to be passed to the interpolate function in the no mesh version
        /// @brief The projection used
        Projection m_projection = Projection::cartesian;

        /// @brief Interpolation parameter
        InterpolationParameters m_interpolationParameters;

        /// @brief Averaging strategy
        std::unique_ptr<averaging::AveragingStrategy> m_strategy;

        /// @brief Map from sample id (int) to sample data.
        std::map<int, std::vector<double>> m_sampleData;

        /// @brief RTree of mesh nodes
        std::unique_ptr<RTreeBase> m_nodeRTree;
    };

} // namespace meshkernel

inline meshkernel::UInt meshkernel::SampleAveragingInterpolator::Size() const
{
    return static_cast<UInt>(m_samplePoints.size());
}
