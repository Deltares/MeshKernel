using System.Runtime.InteropServices;

namespace DeltaShell.Plugins.GridEditor.GridGeomStateful.Native
{
    [StructLayout(LayoutKind.Sequential)]
    public struct CurvilinearParametersNative
    {
        /// <summary>
        /// *M-refinement factor for regular grid generation (2000) 
        /// </summary>
        public int MRefinement;

        /// <summary>
        /// *N-refinement factor for regular grid generation (40) 
        /// </summary>
        public int NRefinement;

        /// <summary>
        /// Nr. of inner iterations in regular grid smoothing (10).
        /// </summary>
        public int SmoothingIterations;

        /// <summary>
        /// * Smoothing parameter (0.5).
        /// </summary>
        public double SmoothingParameter;

        /// <summary>
        /// *Attraction/repulsion parameter (0).
        /// </summary>
        public double AttractionParameter;
    }
}
