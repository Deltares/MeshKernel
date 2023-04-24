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

#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/MeshInterpolationInterface.hpp>

namespace meshkernel
{
    /// @brief A class for performing bilinear interpolation on gridded samples
    class BilinearInterpolationOnGriddedSamples : public MeshInterpolationInterface
    {
    public:
        /// @brief Bilinear interpolation with constant cell size (faster because no linear search is performed for each mesh node)
        /// @param[in] mesh The input mesh
        /// @param[in] numColumns The number of grid columns
        /// @param[in] numRows The number of grid rows
        /// @param[in] xOrigin The x coordinate of the grid origin
        /// @param[in] yOrigin The y coordinate of the grid origin
        /// @param[in] cellSize The grid cell size
        /// @param[in] values The values of the gridded samples
        BilinearInterpolationOnGriddedSamples(std::shared_ptr<Mesh2D> mesh,
                                              size_t numColumns,
                                              size_t numRows,
                                              double xOrigin,
                                              double yOrigin,
                                              double cellSize,
                                              const std::vector<double>& values);

        /// @brief Bilinear interpolation with non constant cell size (slower because linear search is performed for each mesh node)
        /// @param[in] mesh The input mesh
        /// @param[in] xCoordinates The x coordinates of the grid
        /// @param[in] yCoordinates The y coordinates of the grid
        /// @param[in] values The values of the gridded samples
        BilinearInterpolationOnGriddedSamples(std::shared_ptr<Mesh2D> mesh,
                                              const std::vector<double>& xCoordinates,
                                              const std::vector<double>& yCoordinates,
                                              const std::vector<double>& values);

        /// @brief Compute interpolation
        void Compute() override;

        /// @brief Gets the interpolation value at a specific node
        /// @param[in] node The node index
        /// @return The interpolated value
        [[nodiscard]] double GetNodeResult(size_t node) const override { return m_nodeResults[node]; }

        /// @brief Gets the interpolation value at a specific edge
        /// @param[in] edge The edge index
        /// @return The interpolated value
        [[nodiscard]] double GetEdgeResult(size_t edge) const override { return m_edgeResults[edge]; }

        /// @brief Gets the interpolation value at a specific face
        /// @param[in] face The face index
        /// @return  The interpolated value
        [[nodiscard]] double GetFaceResult(size_t face) const override { return m_faceResults[face]; }

        /// @brief Gets all interpolated values at nodes
        /// @return The interpolated values
        [[nodiscard]] const std::vector<double>& GetNodeResults() const override { return m_nodeResults; }

        /// @brief Gets all interpolated values at edges
        /// @return The interpolated values
        [[nodiscard]] const std::vector<double>& GetEdgeResults() const override { return m_edgeResults; }

        /// @brief Gets all interpolated values at faces
        /// @return The interpolated values
        [[nodiscard]] const std::vector<double>& GetFaceResults() const override { return m_faceResults; }

    private:
        /// @brief Performs bilinear interpolation
        /// @param[in] point The input point
        /// @return The result of bilinear interpolation at the point
        double bilinearInterpolation(const Point& point) const;

        /// @brief For a specific point, gets the fractional number of columns
        /// @param[in] point The input point
        /// @return The fractional column index
        [[nodiscard]] inline double GetFractionalNumberOfColumns(const Point& point) const;

        /// @brief For a specific point, gets the fractional number of columns
        /// @param[in] point The input point
        /// @return The fractional row index
        [[nodiscard]] inline double GetFractionalNumberOfRows(const Point& point) const;

        /// @brief Gets the sample value at specific row and column
        /// @return The sample value
        [[nodiscard]] double getGriddedValue(size_t columnIndex, size_t rowIndex) const
        {
            const auto index = rowIndex * (m_numColumns + 1) + columnIndex;
            return m_values[index];
        }

        const std::shared_ptr<Mesh2D> m_mesh; ///< Pointer to the mesh

        size_t m_numColumns; ///< The number of grid columns
        size_t m_numRows;    ///< The number of grid rows
        double m_xOrigin;    ///< The x coordinate of the grid origin
        double m_yOrigin;    ///< The y coordinate of the grid origin
        double m_cellSize;   ///< The grid cell size

        std::vector<double> m_xCoordinates; ///< The x coordinates of the grid
        std::vector<double> m_yCoordinates; ///< The y coordinates of the grid
        std::vector<double> m_values;       ///< The gridded sample values
        bool m_isCellSizeConstant;          ///< If the grid coordinates are specified using vectors of coordinates

        std::vector<double> m_nodeResults; ///< The interpolation results at nodes
        std::vector<double> m_edgeResults; ///< The interpolation results at edges
        std::vector<double> m_faceResults; ///< The interpolation results at faces
    };
} // namespace meshkernel
