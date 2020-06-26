#pragma once

namespace GridGeomApi
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
}


