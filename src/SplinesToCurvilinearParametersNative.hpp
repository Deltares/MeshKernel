#pragma once

namespace GridGeomApi
{
    struct SplinesToCurvilinearParametersNative
    {
        /// <summary>
        /// * Aspect ratio (0.1)
        /// </summary>
        double AspectRatio;

        /// <summary>
        /// * Grow factor of aspect ratio (1.1)
        /// </summary>
        double AspectRatioGrowFactor;

        /// <summary>
        /// * Average mesh width on center spline (0.005)
        /// </summary>
        double AverageWidth;

        /// <summary>
        /// * Curvature adapted grid spacing, 1 or not 0 (1)
        /// </summary>
        int CurvatureAdapetedGridSpacing;

        /// <summary>
        /// * Grow the grid outside the prescribed grid height (1)
        /// </summary>
        int GrowGridOutside;

        /// <summary>
        /// * Maximum number of layers in the uniform part (5)
        /// </summary>
        int MaximumNumberOfGridCellsInTheUniformPart;

        /// <summary>
        /// * On-top-of-each-other tolerance (0.0001)
        /// </summary>
        double GridsOnTopOfEachOtherTolerance;

        /// <summary>
        /// * Minimum allowed absolute value of crossing-angle cosine (0.95)
        /// </summary>
        double MinimumCosineOfCrossingAngles;

        /// <summary>
        /// * Check for collisions with other parts of the front, 1 or not 0 (0)
        /// </summary>
        int CheckFrontCollisions;

        /// <summary>
        /// * Uniform grid size, netboundary to grid only (0.0)
        /// </summary>
        double UniformGridSize;

        /// <summary>
        /// * Remove skinny triangles (1)
        /// </summary>
        int RemoveSkinnyTriangles;
    };
}


