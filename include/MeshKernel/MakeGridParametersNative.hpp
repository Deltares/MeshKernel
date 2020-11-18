//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2020.
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

namespace meshkernelapi
{
    struct MakeGridParametersNative
    {
        /// <summary>
        /// The type of grid to create : square = 0, wieber = 1, hexagonal type 1 = 2,  hexagonal type 2 = 3, triangular = 4 (0)
        /// </summary>
        int GridType;

        /// <summary>
        /// The number of columns in x direction
        /// </summary>
        int NumberOfColumns;

        /// <summary>
        /// The number of columns in y direction
        /// </summary>
        int NumberOfRows;

        /// <summary>
        /// The grid angle
        /// </summary>
        double GridAngle;

        /// <summary>
        /// The grid block size, used in x and y direction
        /// </summary>
        double GridBlockSize;

        /// <summary>
        /// The line thickness in mm (interactor setting)
        /// </summary>
        double LineThickness;

        /// <summary>
        /// (interactor setting)
        /// </summary>
        double hSize;

        /// <summary>
        /// The x coordinate of the origin, located at the bottom left corner
        /// </summary>
        double OriginXCoordinate;

        /// <summary>
        /// The y coordinate of the origin, located at the bottom left corner
        /// </summary>
        double OriginYCoordinate;

        /// <summary>
        /// The z coordinate of the origin, located at the bottom left corner
        /// </summary>
        double OriginZCoordinate;

        /// <summary>
        /// The grid block size in x dimension, used only for squared grids
        /// </summary>
        double XGridBlockSize;

        /// <summary>
        /// The grid block size in y dimension, used only for squared grids
        /// </summary>
        double YGridBlockSize;
    };
} // namespace meshkernelapi
