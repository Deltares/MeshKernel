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
#include <MeshKernel/MeshInterpolation.hpp>

namespace meshkernel
{
    /// @brief A class for performing bilinear interpolation on gridded samples
    class BilinearInterpolationOnGriddedSamples : public MeshInterpolation
    {
    public:
        /// @brief Bilinear interpolation with constant cell size (faster because no linear search is performed for each mesh node)
        /// @param[in] mesh The input mesh
        /// @param[in] numXCoord The number of x coordinates of the gridded data
        /// @param[in] numYCoord The number of y coordinates of the gridded data
        /// @param[in] origin The coordinate of the grid origin
        /// @param[in] cellSize The grid cell size
        /// @param[in] values The values of the gridded samples
        BilinearInterpolationOnGriddedSamples(const Mesh2D& mesh,
                                              UInt numXCoord,
                                              UInt numYCoord,
                                              const Point& origin,
                                              double cellSize,
                                              const std::vector<double>& values);

        /// @brief Bilinear interpolation with non constant cell size (slower because linear search is performed for each mesh node)
        /// @param[in] mesh The input mesh
        /// @param[in] xCoordinates The x coordinates of the grid
        /// @param[in] yCoordinates The y coordinates of the grid
        /// @param[in] values The values of the gridded samples
        BilinearInterpolationOnGriddedSamples(const Mesh2D& mesh,
                                              const std::vector<double>& xCoordinates,
                                              const std::vector<double>& yCoordinates,
                                              const std::vector<double>& values);

        /// @brief Compute interpolation
        void Compute() override;

    private:
        /// @brief Performs bilinear interpolation
        /// @param[in] point The input point
        /// @return The result of bilinear interpolation at the point
        double Interpolation(const Point& point) const;

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
        [[nodiscard]] double getGriddedValue(UInt columnIndex, UInt rowIndex) const
        {
            const auto index = rowIndex * m_numXCoord + columnIndex;
            return m_values[index];
        }

        const Mesh2D& m_mesh; ///< Pointer to the mesh
        UInt m_numXCoord;     ///< The number of x coordinates of the gridded data
        UInt m_numYCoord;     ///< The number of y coordinates of the gridded data
        Point m_origin;       ///< The coordinate of the origin
        double m_cellSize;    ///< The grid cell size

        std::vector<double> m_xCoordinates; ///< The x coordinates of the grid
        std::vector<double> m_yCoordinates; ///< The y coordinates of the grid
        std::vector<double> m_values;       ///< The gridded sample values
        bool m_isCellSizeConstant;          ///< If the grid coordinates are specified using vectors of coordinates
    };
} // namespace meshkernel
