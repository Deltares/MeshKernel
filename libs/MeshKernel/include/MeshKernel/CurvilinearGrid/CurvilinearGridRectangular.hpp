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

#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Parameters.hpp>
#include <MeshKernel/Utilities/LinearAlgebra.hpp>

#include <memory>

namespace meshkernel
{
    class Polygons;
    class CurvilinearGrid;

    /// @brief A class implementing the generation of a rectangular curvilinear grid
    class CurvilinearGridRectangular
    {
    public:
        /// @brief Class constructor
        ///
        /// @param[in] projection The projection to use
        CurvilinearGridRectangular(Projection projection);

        /// @brief Compute a rectangular curvilinear grid, given the origin and the number of rows, columns and block size
        /// @param[in] numColumns The number of columns in x direction
        /// @param[in] numRows The number of columns in y direction
        /// @param[in] originX The x coordinate of the origin, located at the bottom left corner
        /// @param[in] originY The y coordinate of the origin, located at the bottom left corner
        /// @param[in] angle The grid angle
        /// @param[in] blockSizeX The grid block size in x dimension
        /// @param[in] blockSizeY The grid block size in y dimension
        /// @returns[in] A curvilinear grid
        std::unique_ptr<CurvilinearGrid> Compute(const int numColumns,
                                                 const int numRows,
                                                 const double originX,
                                                 const double originY,
                                                 const double angle,
                                                 const double blockSizeX,
                                                 const double blockSizeY) const;

        /// @brief Compute a rectangular curvilinear grid in one polygon, given an angle and the block sizes
        /// @param[in] angle The grid angle
        /// @param[in] blockSizeX The grid block size in x dimension
        /// @param[in] blockSizeY The grid block size in y dimension
        /// @returns[in] A curvilinear grid
        std::unique_ptr<CurvilinearGrid> Compute(const double angle,
                                                 const double blockSizeX,
                                                 const double blockSizeY,
                                                 std::shared_ptr<Polygons> polygons,
                                                 UInt polygonIndex) const;

        /// @brief Compute a rectangular curvilinear grid in one polygon, given the block size and the extension. The grid angle is 0
        /// @param[in] originX The x coordinate of the origin, located at the bottom left corner
        /// @param[in] originY The y coordinate of the origin, located at the bottom left corner
        /// @param[in] blockSizeX The grid block size in x dimension
        /// @param[in] blockSizeY The grid block size in y dimension
        /// @param[in] upperRightX The x coordinate of the upper right corner
        /// @param[in] upperRightY The y coordinate of the upper right corner
        /// @returns[in] A curvilinear grid
        std::unique_ptr<CurvilinearGrid> Compute(const double originX,
                                                 const double originY,
                                                 const double blockSizeX,
                                                 const double blockSizeY,
                                                 const double upperRightX,
                                                 const double upperRightY) const;

    private:
        /// @brief Compute a rectangular curvilinear grid on cartesian coordinates.
        /// @param[in] numColumns The number of columns in x direction
        /// @param[in] numRows The number of columns in y direction
        /// @param[in] originX The x coordinate of the origin, located at the bottom left corner
        /// @param[in] originY The y coordinate of the origin, located at the bottom left corner
        /// @param[in] angle The grid angle
        /// @param[in] blockSizeX The grid block size in x dimension
        /// @param[in] blockSizeY The grid block size in y dimension
        /// @returns[in] The coordinates of the grid point
        static lin_alg::Matrix<Point> ComputeCartesian(const int numColumns,
                                                       const int numRows,
                                                       const double originX,
                                                       const double originY,
                                                       const double angle,
                                                       const double blockSizeX,
                                                       const double blockSizeY);

        /// @brief Compute a rectangular curvilinear grid on spherical coordinates.
        /// A correction to the longitudinal discretization is applied to preserve an aspect ratio ds/dy = 1 on real distances.
        /// For preventing the creation of small edges, the correction stops when the distance is less than 2000 meters and the grid is generated around the poles.
        /// This is a customized fix for GTSM models.
        /// @param[in] numColumns The number of columns in x direction
        /// @param[in] numRows The number of columns in y direction
        /// @param[in] originX The x coordinate of the origin, located at the bottom left corner
        /// @param[in] originY The y coordinate of the origin, located at the bottom left corner
        /// @param[in] angle The grid angle
        /// @param[in] blockSizeX The grid block size in x dimension
        /// @param[in] blockSizeY The grid block size in y dimension
        /// @returns[in] The coordinates of the grid point
        static lin_alg::Matrix<Point> ComputeSpherical(const int numColumns,
                                                       const int numRows,
                                                       const double originX,
                                                       const double originY,
                                                       const double angle,
                                                       const double blockSizeX,
                                                       const double blockSizeY);

        /// @brief Compute the adjusted latitude for keeping an aspect ratio of 1, considering the spherical coordinates
        /// @param[in] blockSize The grid block size in y dimension
        /// @param[in] latitude The current latitude
        /// @returns[in] The adjusted latitude
        static double ComputeLatitudeIncrementWithAdjustment(double blockSize, double latitude);

        /// @brief Compute the number of rows required to generate a grid from minY to maxY
        /// @param[in] minY The min latitude
        /// @param[in] maxY The max latitude
        /// @param[in] projection The projection to use
        /// @returns[in] The number of rows
        static int ComputeNumRows(double minY,
                                  double maxY,
                                  double blockSizeY,
                                  Projection projection);

        Projection m_projection; ///< The projection to use
    };
} // namespace meshkernel
