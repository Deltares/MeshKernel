using System.Runtime.InteropServices;

namespace DeltaShell.Plugins.GridEditor.GridGeomStateful.Native
{
    [StructLayout(LayoutKind.Sequential)]
    public struct InterpolationParametersNative
    {
        /// <summary>
        /// *Actual interpolation type (1)
        /// </summary>
        public int InterpolationType;

        /// <summary>
        /// Variable related to interactor behaviour (0) 
        /// </summary>
        public int DisplayInterpolationProcess;

        /// <summary>
        /// *Maximum number of refinement iterations, set to 1 if only one refinement is wanted (10) 
        /// </summary>
        public int MaxNumberOfRefinementIterations;

        /// <summary>
        /// *Averaging method : 1 = simple averaging, 2 = closest point, 3 = max, 4 = min, 5 = inverse weighted distance, 6 = minabs, 7 = kdtree (1)
        /// </summary>
        public int AveragingMethod;

        /// <summary>
        /// *Minimum number of points needed inside cell to handle the cell (1)
        /// </summary>
        public int MinimumNumberOfPoints;

        /// <summary>
        /// *Relative search cell size, default 1= actual cell size, 2= twice as large, search radius can be larger than cell so more sample are included. (1.01)
        /// </summary>
        public double RelativeSearchRadius;

        /// <summary>
        /// *Interpolation settings, 1=bathy, 2=zk, 3=s1, 4=Zc (2)
        /// </summary>
        public int InterpolateTo;
    }
}
