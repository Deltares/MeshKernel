//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2021.
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

#include "MeshInterpolationInterface.hpp"

#include <MeshKernel/AveragingStrategies/AveragingStrategy.hpp>
#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/RTree.hpp>

namespace meshkernel
{
    // Forward declarations
    class Mesh;
    struct Sample;

    /// @brief The class used to interpolate based on averaging
    ///
    /// The averaging interpolation operates on three specific \ref Mesh::Location - Faces
    /// (m_facesMassCenters), Nodes, and Edges(m_edgesCenters). The idea is to
    /// collect all samples close to the locations and perform a mathematical
    /// operation on their values. The \ref Method enum describes available operations.
    ///
    /// The algorithm operates as follow:
    ///
    /// -   The samples are ordered in an RTree for a fast search
    ///
    /// -   The search area around the location is constructed as follow:
    ///
    ///     1.  For face locations, the sample areas correspond to the faces,
    ///         increased/decreased by the relativeSearchRadius parameter
    ///         (relativeSearchRadius > 1 increased, relativeSearchRadius < 1
    ///         decreased).
    ///
    ///     2.  For \ref Mesh::Location Nodes and \ref Mesh::Location Edges locations, the dual face around the node is
    ///         constructed by connecting the mid-points of all edges connected
    ///         to the node. As above, the resulting polygon can be
    ///         increased/decreased by the relativeSearchRadius parameter.
    ///
    /// -   The search radius is computed from the constructed polygon nodes and
    ///     the locations, the maximum value is taken.
    ///
    /// -   The sample RTree is inquired to retrieve the indices of the samples
    ///     within the search radius.
    ///
    /// -   The operations described above are executed on the found samples.
    ///
    /// -   For the \ref Mesh::Location Edges location, the interpolated values at the node are
    ///     averaged.
    class AveragingInterpolation : public MeshInterpolationInterface
    {
    public:
        /// @brief Averaging methods
        enum class Method
        {
            SimpleAveraging = 1,         ///< Computes a simple mean
            Closest = 2,                 ///< Takes the value of the closest sample to the interpolation location
            Max = 3,                     ///< Takes the maximum sample value
            Min = 4,                     ///< Takes the minimum sample value
            InverseWeightedDistance = 5, ///< Computes the inverse weighted sample mean
            MinAbsValue = 6              ///< Computes the minimum absolute value
        };

        /// @brief Interpolation based on averaging
        /// @param[in] mesh                            The input mesh
        /// @param[in] samples                         The samples with x,y locations and values
        /// @param[in] method                          The averaging method to use
        /// @param[in] locationType                    The location type (faces, edges, nodes).
        /// @param[in] relativeSearchRadius            The relative search radius, used to enlarge the search area when looking for samples.
        /// @param[in] useClosestSampleIfNoneAvailable If no samples are found, use the closest one.
        /// @param[in] subtractSampleValues            For some algorithms (e.g. refinement based on levels) we need to subtract 1 to the sample value.
        /// @param[in] minNumSamples                   The minimum a of samples used for certain interpolation algorithms
        AveragingInterpolation(std::shared_ptr<Mesh2D> mesh,
                               std::vector<Sample>& samples,
                               Method method,
                               Mesh::Location locationType,
                               double relativeSearchRadius,
                               bool useClosestSampleIfNoneAvailable,
                               bool subtractSampleValues,
                               size_t minNumSamples);

        /// @brief Compute interpolation
        void Compute() override;

        [[nodiscard]] double GetNodeResult(size_t node) const override { return m_nodeResults[node]; }
        [[nodiscard]] double GetEdgeResult(size_t edge) const override { return m_edgeResults[edge]; }
        [[nodiscard]] double GetFaceResult(size_t face) const override { return m_faceResults[face]; }

        [[nodiscard]] const std::vector<double>& GetNodeResults() const override { return m_nodeResults; }
        [[nodiscard]] const std::vector<double>& GetEdgeResults() const override { return m_edgeResults; }
        [[nodiscard]] const std::vector<double>& GetFaceResults() const override { return m_faceResults; }

    private:
        /// @brief Compute the averaging results in polygon
        /// @param[in]  polygon            The bounding polygon where the samples are included
        /// @param[in]  interpolationPoint The interpolation point
        /// @returns The resulting value
        double ComputeOnPolygon(const std::vector<Point>& polygon,
                                Point interpolationPoint);

        /// @brief Decreases the values of samples
        void DecreaseValueOfSamples();

        /// @brief Generate the search polygon from an input polygon
        /// @param[in]  polygon            The input polygon
        /// @param[in]  interpolationPoint The interpolation point
        /// @return The search polygon
        [[nodiscard]] std::vector<Point> GetSearchPolygon(std::vector<Point> const& polygon, Point const& interpolationPoint) const;

        /// @brief Computes the average value from the neighbors using a strategy
        /// param[in]  strategy            The input strategy
        /// param[in] searchPolygon        The bounding polygon
        /// @return The interpolated result
        [[nodiscard]] double ComputeInterpolationResultFromNeighbors(std::unique_ptr<averaging::AveragingStrategy> strategy, std::vector<Point> const& searchPolygon);

        /// @brief Gets the sample value from an r-tree query
        /// param[in] index            The query index
        /// @return The sample value
        [[nodiscard]] double GetSampleValueFromRTree(size_t index);

        /// @brief Compute a search radius from a point and a polygon
        /// @param searchPolygon The input polygon
        /// @param interpolationPoint The input point
        /// @return the search radius including the polygon from the input point
        [[nodiscard]] double GetSearchRadiusSquared(std::vector<Point> const& searchPolygon,
                                                    Point const& interpolationPoint) const;

        const std::shared_ptr<Mesh2D> m_mesh;           ///< Pointer to the mesh
        std::vector<Sample>& m_samples;                 ///< The samples
        Method m_method;                                ///< The method to use for the interpolation
        Mesh::Location m_interpolationLocation;         ///< Interpolation location
        double m_relativeSearchRadius;                  ///< Relative search radius
        bool m_useClosestSampleIfNoneAvailable = false; ///< Whether to use the closest sample if there is none available
        bool m_transformSamples = false;                ///< Wheher to transform samples
        size_t m_minNumSamples = 1;                     ///< The minimum amount of samples for a valid interpolation. Used in some interpolation algorithms.

        std::vector<double> m_nodeResults; ///< The interpolation results at nodes
        std::vector<double> m_edgeResults; ///< The interpolation results at edges
        std::vector<double> m_faceResults; ///< The interpolation results at faces

        std::vector<bool> m_visitedSamples; ///< The visited samples

        RTree m_samplesRtree; ///< The samples tree
    };
} // namespace meshkernel
