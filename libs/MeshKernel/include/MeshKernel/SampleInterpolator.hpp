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

    class SampleInterpolator
    {
    public:
        virtual ~SampleInterpolator() = default;

        virtual UInt Size() const = 0;

        /// @brief Set sample data
        template <meshkernel::ValidConstDoubleArray VectorType>
        void SetData(const int propertyId, const VectorType& sampleData)
        {
            const std::span<const double> spanSampleData(sampleData.data(), sampleData.size());
            SetDataSpan(propertyId, sampleData);
        }

        /// @brief Set sample data contained in an std::span object
        virtual void SetDataSpan(const int propertyId, const std::span<const double>& sampleData) = 0;

        /// @brief Interpolate the sample data set at the interpolation nodes.
        template <meshkernel::ValidConstPointArray PointVectorType, meshkernel::ValidConstDoubleArray ScalarVectorType>
        void Interpolate(const int propertyId, const PointVectorType& iterpolationNodes, ScalarVectorType& result) const
        {
            const std::span<const Point> spanInterpolationNodes(iterpolationNodes.data(), iterpolationNodes.size());
            std::span<double> spanResult(result.data(), result.size());
            InterpolateSpan(propertyId, spanInterpolationNodes, spanResult);
        }

        /// @brief Interpolate the sample data set at the interpolation nodes.
        virtual void InterpolateSpan(const int propertyId, const std::span<const Point>& iterpolationNodes, std::span<double>& result) const = 0;

        /// @brief Interpolate the sample data set at the interpolation nodes.
        template <meshkernel::ValidConstDoubleArray ScalarVectorType>
        void Interpolate(const int propertyId, const Mesh2D& mesh, const Location location, ScalarVectorType& result) const
        {
            std::span<double> spanResult(result.data(), result.size());
            InterpolateSpan(propertyId, mesh, location, spanResult);
        }

        /// @brief Interpolate the sample data.
        // TODO may need a polygon too: const Polygons& polygon
        virtual void InterpolateSpan(const int propertyId, const Mesh2D& mesh, const Location location, std::span<double>& result) const = 0;

        /// @brief Interpolate the sample data set at a single interpolation point.
        ///
        /// If interpolation at multiple points is required then better performance
        /// can be obtained using the Interpolate function above.
        virtual double InterpolateValue(const int propertyId, const Point& evaluationPoint) const = 0;

        /// @brief Determine if the SampleInterpolator already has this sample set.
        virtual bool Contains(const int propertyId) const = 0;

    private:
        // TODO the map -> sample-data can be here
        // with protected members to access the data.
    };

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

        void SetDataSpan(const int propertyId, const std::span<const double>& sampleData) override;

        /// @brief Set sample data
        // template <meshkernel::ValidConstDoubleArray VectorType>
        // void SetData(const int propertyId, const VectorType& sampleData) override;

        /// @brief Interpolate the sample data set at the interpolation nodes.
        // template <meshkernel::ValidConstPointArray PointVectorType, meshkernel::ValidConstDoubleArray ScalarVectorType>
        // void Interpolate(const int propertyId, const PointVectorType& iterpolationNodes, ScalarVectorType& result) const override;
        void InterpolateSpan(const int propertyId, const std::span<const Point>& iterpolationNodes, std::span<double>& result) const override;

        // DOes nothing at the moment.
        void InterpolateSpan(const int propertyId [[maybe_unused]], const Mesh2D& mesh [[maybe_unused]], const Location location [[maybe_unused]], std::span<double>& result [[maybe_unused]]) const override
        {
        }

        /// @brief Interpolate the sample data set at a single interpolation point.
        ///
        /// If interpolation at multiple points is required then better performance
        /// can be obtained using the Interpolate function above.
        double InterpolateValue(const int propertyId, const Point& evaluationPoint) const override;

        /// @brief Determine if the SampleInterpolator already has this sample set.
        bool Contains(const int propertyId) const override;

    private:
        /// @brief Interpolate the sample data on the element at the interpolation point.
        double InterpolateOnElement(const UInt elementId, const Point& interpolationPoint, const std::vector<double>& sampleValues) const;

        /// @brief Triangulation of the sample points
        MeshTriangulation m_triangulation;

        /// @brief Map from sample id (int) to sample data.
        std::map<int, std::vector<double>> m_sampleData;
    };

    struct InterpolationParameters
    {
        AveragingInterpolation::Method m_method = AveragingInterpolation::Method::SimpleAveraging;
        double m_relativeSearchRadius = 0.0;
        bool m_useClosestIfNoneFound = true;
        UInt m_minimumNumberOfSamples = 10;
    };

    /// @brief Interpolator for sample data
    class SampleAveragingInterpolator : public SampleInterpolator
    {
    public:
        /// @brief Constructor.
        ///
        /// The VectorType can be any array type of double precision values, e.g. std::vector, std::span.
        // TODO need InterpolationParameters
        // containing e.g.
        //   - averaging method
        //   - relative search radius
        //   - use closest if non found
        //   - min number of samples
        template <meshkernel::ValidConstDoubleArray VectorType>
        SampleAveragingInterpolator(const VectorType& xNodes,
                                    const VectorType& yNodes,
                                    const Projection projection,
                                    const InterpolationParameters& interpolationParameters)
            : m_samplePoints(CombineCoordinates(xNodes, yNodes)),
              m_projection(projection),
              m_interpolationParameters(interpolationParameters),
              m_strategy(averaging::AveragingStrategyFactory::GetAveragingStrategy(interpolationParameters.m_method,
                                                                                   interpolationParameters.m_minimumNumberOfSamples,
                                                                                   projection)),
              m_nodeRTree(RTreeFactory::Create(projection))
        {
            m_nodeRTree->BuildTree(m_samplePoints);
        }

        /// @brief Constructor.
        ///
        /// The VectorType can be any array type of double precision values, e.g. std::vector, std::span.
        template <meshkernel::ValidConstPointArray PointVector>
        SampleAveragingInterpolator(const PointVector& nodes,
                                    const Projection projection,
                                    const InterpolationParameters& interpolationParameters)
            : m_samplePoints(nodes.begin(), nodes.end()),
              m_projection(projection),
              m_interpolationParameters(interpolationParameters),
              m_strategy(averaging::AveragingStrategyFactory::GetAveragingStrategy(interpolationParameters.m_method,
                                                                                   interpolationParameters.m_minimumNumberOfSamples,
                                                                                   projection)),
              m_nodeRTree(RTreeFactory::Create(projection))
        {
            m_nodeRTree->BuildTree(m_samplePoints);
        }

        /// @brief Get the number of nodes of size of the sample data.
        UInt Size() const override;

        void SetDataSpan(const int propertyId, const std::span<const double>& sampleData) override;

        /// @brief Interpolate the sample data set at the interpolation nodes.
        void InterpolateSpan(const int propertyId, const std::span<const Point>& iterpolationNodes, std::span<double>& result) const override;

        void InterpolateSpan(const int propertyId, const Mesh2D& mesh, const Location location, std::span<double>& result) const override;

        /// @brief Interpolate the sample data set at a single interpolation point.
        ///
        /// If interpolation at multiple points is required then better performance
        /// can be obtained using the Interpolate function above.
        double InterpolateValue(const int propertyId, const Point& evaluationPoint) const override;

        /// @brief Determine if the SampleInterpolator already has this sample set.
        bool Contains(const int propertyId) const override;

    private:
        static constexpr UInt MaximumNumberOfEdgesPerNode = 16; ///< Maximum number of edges per node

        template <meshkernel::ValidConstDoubleArray VectorType>
        static std::vector<Point> CombineCoordinates(const VectorType& xNodes, const VectorType& yNodes)
        {
            std::vector<Point> result(xNodes.size());

            for (size_t i = 0; i < xNodes.size(); ++i)
            {
                result[i].x = xNodes[i];
                result[i].y = yNodes[i];
            }

            return result;
        }

        /// @brief Interpolate the sample data on the element at the interpolation point.
        double InterpolateOnElement(const UInt elementId, const Point& interpolationPoint, const std::vector<double>& sampleValues) const;

        double ComputeOnPolygon(const int propertyId,
                                std::vector<Point>& polygon,
                                const Point& interpolationPoint,
                                const double relativeSearchRadius,
                                const bool useClosestSampleIfNoneAvailable,
                                const Projection projection,
                                std::vector<Sample>& sampleCache) const;

        double GetSearchRadiusSquared(const std::vector<Point>& searchPolygon,
                                      const Point& interpolationPoint,
                                      const Projection projection) const;

        void GenerateSearchPolygon(const double relativeSearchRadius,
                                   const Point& interpolationPoint,
                                   std::vector<Point>& polygon,
                                   const Projection projection) const;

        double GetSampleValueFromRTree(const int propertyId, const UInt index) const;

        double ComputeInterpolationResultFromNeighbors(const int propertyId,
                                                       const Point& interpolationPoint,
                                                       const std::vector<Point>& searchPolygon,
                                                       const Projection projection,
                                                       std::vector<Sample>& sampleCache) const;

        std::vector<Point> m_samplePoints;

        // SHould use the m_projection from the mesh in the interpolate function or
        // needs to be passed to the interpolate function in the no mesh version
        Projection m_projection = Projection::cartesian; ///< The projection used

        InterpolationParameters m_interpolationParameters;

        std::unique_ptr<averaging::AveragingStrategy> m_strategy; ///< Averaging strategy

        /// @brief Map from sample id (int) to sample data.
        std::map<int, std::vector<double>> m_sampleData;

        ///< RTree of mesh nodes
        std::unique_ptr<RTreeBase> m_nodeRTree;
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

inline meshkernel::UInt meshkernel::SampleAveragingInterpolator::Size() const
{
    return static_cast<UInt>(m_samplePoints.size());
}

inline bool meshkernel::SampleAveragingInterpolator::Contains(const int propertyId) const
{
    return m_sampleData.contains(propertyId);
}
