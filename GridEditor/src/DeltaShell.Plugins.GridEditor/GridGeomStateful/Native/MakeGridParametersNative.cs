using System.Runtime.InteropServices;

namespace DeltaShell.Plugins.GridEditor.GridGeomStateful.Native
{
    [StructLayout(LayoutKind.Sequential)]
    public struct MakeGridParametersNative
    {
        /// <summary>
        /// * The type of grid to create : square = 0, wieber = 1, hexagonal type 1 = 2,  hexagonal type 2 = 3, triangular = 4 (0) 
        /// </summary>  
        public int GridType;

        /// <summary>   
        /// * The number of columns in x direction (3)
        /// </summary>
        public int NumberOfColumns;

        /// <summary>
        /// * The number of columns in y direction (3)
        /// </summary>
        public int NumberOfRows;

        /// <summary>
        /// * The grid angle (0.0)
        /// </summary>
        public double GridAngle;

        /// <summary>
        /// * The grid block size, used in x and y direction (50.0)
        /// </summary>
        public double GridBlockSize;

        /// <summary>
        /// Interactor setting: line thickness in mm (8.0)
        /// </summary>
        public double LineThickness;

        /// <summary>
        /// Interactor setting: line thickness in mm (8.0)
        /// </summary>
        public double hSize;

        /// <summary>
        /// * The x coordinate of the origin, located at the bottom left corner (0.0)	 
        /// </summary>
        public double OriginXCoordinate;

        /// <summary>
        /// * The y coordinate of the origin, located at the bottom left corner (0.0)	 
        /// </summary>
        public double OriginYCoordinate;

        /// <summary>
        /// * The z coordinate of the origin, located at the bottom left corner (0.0)	 
        /// </summary>
        public double OriginZCoordinate;

        /// <summary>
        /// * The grid block size in x dimension, used only for squared grids (10.0) 
        /// </summary>
        public double XGridBlockSize;

        /// <summary>
        /// * The grid block size in y dimension, used only for squared grids (10.0) 
        /// </summary>
        public double YGridBlockSize;
    }
}
