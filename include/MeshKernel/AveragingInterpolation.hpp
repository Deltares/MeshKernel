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
    /// The averaging interpolation operates on three specific \ref MeshLocations - Faces
    /// (m face mass centers), Nodes, and Edges(m edge centers). The idea is to
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
    ///     2.  For \ref MeshLocations Nodes and \ref MeshLocations Edges locations, the dual face around the node is
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
    /// -   For the \ref MeshLocations Edges location, the interpolated values at the node are
    ///     averaged.
    class AveragingInterpolation
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
        /// @param[in] mesh The input mesh
        /// @param[in] samples The samples with x,y locations and values
        /// @param[in] method The averaging method to use
        /// @param[in] locationType The location type (faces, edges, nodes).
        /// @param[in] relativeSearchRadius The relative search radius, used to enlarge the search area when looking for samples.
        /// @param[in] useClosestSampleIfNoneAvailable If no samples are found use the closest one.
        /// @param[in] subtractSampleValues For some algorithms (e.g. refinement based on levels) we need to subtract 1 to the sample value.
        explicit AveragingInterpolation(std::shared_ptr<Mesh2D> mesh,
                                        std::vector<Sample>& samples,
                                        Method method,
                                        MeshLocations locationType,
                                        double relativeSearchRadius,
                                        bool useClosestSampleIfNoneAvailable,
                                        bool subtractSampleValues);

        /// @brief Compute interpolation
        void Compute();

        /// @brief Get the result values
        /// @return the results
        [[nodiscard]] const auto& GetResults() const
        {
            return m_results;
        }

    private:
        /// @brief Compute the averaging results in polygon
        /// @param[in] polygon The bounding polygon where the samples are included
        /// @param[in] interpolationPoint The interpolation point
        /// @param[out] result The resulting value
        void ComputeOnPolygon(const std::vector<Point>& polygon,
                              Point interpolationPoint,
                              double& result);

        /// @brief Compute the interpolated results on designed location
        /// @return the interpolated results
        [[nodiscard]] std::vector<double> ComputeOnLocations();

        /// @brief Compute the interpolated results on faces
        /// @return the interpolated results
        [[nodiscard]] std::vector<double> ComputeOnFaces();

        /// @brief Compute the interpolated results on nodes or edges
        /// @return the interpolated results
        [[nodiscard]] std::vector<double> ComputeOnNodesOrEdges();

        /// @brief Decreases the values of samples
        void DecreaseValueOfSamples();

        const std::shared_ptr<Mesh2D> m_mesh;           ///< Pointer to the mesh
        std::vector<Sample>& m_samples;                 ///< The samples
        Method m_method;                                ///< The method to use for the interpolation
        MeshLocations m_interpolationLocation;          ///< Interpolation location
        double m_relativeSearchRadius;                  ///< Relative search radius
        bool m_useClosestSampleIfNoneAvailable = false; ///< Whether to use the closest sample if there is none available
        bool m_transformSamples = false;                ///< Wheher to transform samples

        RTree m_samplesRtree;               ///< The samples tree
        std::vector<double> m_results;      ///< The results
        std::vector<bool> m_visitedSamples; ///< The visited samples
    };
} // namespace meshkernel
