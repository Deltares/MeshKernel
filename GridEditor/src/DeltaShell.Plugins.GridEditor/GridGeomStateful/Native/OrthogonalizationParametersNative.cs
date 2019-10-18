namespace DeltaShell.Plugins.GridEditor.GridGeomStateful.Native
{
    /// <summary>
    /// All orthogonalization options
    /// </summary>
    public struct OrthogonalizationParametersNative
    {
        /// <summary>
        /// *Number of outer iterations in orthogonalization. Increase this parameter for complex grids (2) 
        /// </summary>
        public int OuterIterations { get; set; }

        /// <summary>
        /// *Number of boundary iterations in grid/net orthogonalization within itatp (25) 
        /// </summary>
        public int BoundaryIterations { get; set; }

        /// <summary>
        /// *Number of inner    iterations in grid/net orthogonalization within itbnd (25)
        /// </summary>
        public int InnerIterations { get; set; }

        /// <summary>
        /// *Factor from 0 to 1. between grid smoothing and grid orthogonality (0.975)	 
        /// </summary>
        public double OrthogonalizationToSmoothingFactor { get; set; }

        /// <summary>
        /// Minimum ATPF on the boundary (1.0) 
        /// </summary>
        public double AtpfB { get; set; }

        /// <summary>
        /// Factor to weight between circumcentres 1.0 and masscentre 0.0 (1.0)
        /// </summary>
        public double Circumormasscenter { get; set; }

        /// <summary>
        /// Factor between smoother 1d0 and area-homogenizer 0d0 (1.0)
        /// </summary>
        public double Smoothorarea { get; set; }

        /// <summary>
        /// Mesh-adaptation method; 0: Winslow, 1: arc-length, 2: harmonic map (1) 
        /// </summary>
        public int AdaptMethod { get; set; }

        /// <summary>
        /// Mesh-refinement factor; between 0d0 and 1d0 (0.0)
        /// </summary>
        public double AdaptBeta { get; set; }

        /// <summary>
        /// Number of smoothing iterations of `solution` u in adaptation (0) 
        /// </summary>
        public int AdaptNiterU { get; set; }

        /// <summary>
        /// Number of smoothing iterations of monitor matrix G in adaptation (4) 
        /// </summary>
        public int AdaptNiterG { get; set; }

        /// <summary>
        /// Curvi-linear-like 0d0 or pure 1d0 orthogonalisation (0.5)
        /// </summary>
        public double OrthoPure { get; set; }
    }
}