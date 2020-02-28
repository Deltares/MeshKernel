#pragma once

namespace GridGeomApi
{
    struct MakeGridParametersNative
    {
        /// <summary>
        /// * The type of grid to create : square = 0, wieber = 1, hexagonal type 1 = 2,  hexagonal type 2 = 3, triangular = 4 (0) 
        /// </summary>  
        int GridType;

        /// <summary>   
        /// * The number of columns in x direction (3)
        /// </summary>
        int NumberOfColumns;

        /// <summary>
        /// * The number of columns in y direction (3)
        /// </summary>
        int NumberOfRows;

        /// <summary>
        /// * The grid angle (0.0)
        /// </summary>
        double GridAngle;

        /// <summary>
        /// * The grid block size, used in x and y direction (50.0)
        /// </summary>
        double GridBlockSize;

        /// <summary>
        /// Interactor setting: line thickness in mm (8.0)
        /// </summary>
        double LineThickness;

        /// <summary>
        /// Interactor setting: line thickness in mm (8.0)
        /// </summary>
        double hSize;

        /// <summary>
        /// * The x coordinate of the origin, located at the bottom left corner (0.0)	 
        /// </summary>
        double OriginXCoordinate;

        /// <summary>
        /// * The y coordinate of the origin, located at the bottom left corner (0.0)	 
        /// </summary>
        double OriginYCoordinate;

        /// <summary>
        /// * The z coordinate of the origin, located at the bottom left corner (0.0)	 
        /// </summary>
        double OriginZCoordinate;

        /// <summary>
        /// * The grid block size in x dimension, used only for squared grids (10.0) 
        /// </summary>
        double XGridBlockSize;

        /// <summary>
        /// * The grid block size in y dimension, used only for squared grids (10.0) 
        /// </summary>
        double YGridBlockSize;
    };
}


