using System.Runtime.InteropServices;

namespace DeltaShell.Plugins.GridEditor.GridGeomStateful.Native
{
    [StructLayout(LayoutKind.Sequential)]
    public struct SampleRefineParametersNative
    {
        /// <summary>
        /// Sample vector dimension (1) 
        /// </summary>
        public int SampleVectorDimension;

        /// <summary>
        /// *Maximum number of refinement iterations, set to 1 if only one refinement is wanted (10) 
        /// </summary>
        public int MaxNumberOfRefinementIterations;

        /// <summary>
        /// *Minimum cell size (50000.0)
        /// </summary>
        public double MinimumCellSize;

        /// <summary>
        /// *Directional refinement, 1 yes 0 no (0)
        /// </summary>
        public int DirectionalRefinement;

        /// <summary>
        /// *Refinement criterion type (2)
        /// </summary>
        public int RefinementType;

        /// <summary>
        /// *Connect hanging nodes at the end of the iteration, 1 yes 0 no (1) 
        /// </summary>
        public int ConnectHangingNodes;

        /// <summary>
        /// *Maximum time-step in courant grid (0.0)
        /// </summary>
        public double MaximumTimeStepInCourantGrid;

        /// <summary>
        /// *Take samples outside cell into account , 1 yes 0 no (1)
        /// </summary>
        public int AccountForSamplesOutside;
    }
}
