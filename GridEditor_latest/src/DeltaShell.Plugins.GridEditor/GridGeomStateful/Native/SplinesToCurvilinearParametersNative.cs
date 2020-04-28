namespace DeltaShell.Plugins.GridEditor.GridGeomStateful.Native
{
    public struct SplinesToCurvilinearParametersNative
    {
        /// <summary>
        /// * Aspect ratio (0.1)
        /// </summary>
        public double AspectRatio { get; set; }

        /// <summary>
        /// * Grow factor of aspect ratio (1.1)
        /// </summary>
        public double AspectRatioGrowFactor { get; set; }

        /// <summary>
        /// * Average mesh width on center spline (0.005)
        /// </summary>
        public double AverageWidth { get; set; }

        /// <summary>
        /// * Curvature adapted grid spacing, 1 or not 0 (1)
        /// </summary>
        public int CurvatureAdapetedGridSpacing { get; set; }

        /// <summary>
        /// * Grow the grid outside the prescribed grid height (1)
        /// </summary>
        public int GrowGridOutside { get; set; }

        /// <summary>
        /// * Maximum number of layers in the uniform part (5)
        /// </summary>
        public int MaximumNumberOfGridCellsInTheUniformPart { get; set; }

        /// <summary>
        /// * On-top-of-each-other tolerance (0.0001)
        /// </summary>
        public double GridsOnTopOfEachOtherTolerance { get; set; }

        /// <summary>
        /// * Minimum allowed absolute value of crossing-angle cosine (0.95)
        /// </summary>
        public double MinimumCosineOfCrossingAngles { get; set; }

        /// <summary>
        /// * Check for collisions with other parts of the front, 1 or not 0 (0)
        /// </summary>
        public int CheckFrontCollisions { get; set; }

        /// <summary>
        /// * Uniform grid size, netboundary to grid only (0.0)
        /// </summary>
        public double UniformGridSize { get; set; }

        /// <summary>
        /// * Remove skinny triangles (1)
        /// </summary>
        public int RemoveSkinnyTriangles { get; set; }
    }
}
