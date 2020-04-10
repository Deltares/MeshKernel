#pragma once

namespace GridGeomApi
{
    struct InterpolationParametersNative
    {
        /// *Actual interpolation type (1)
        int InterpolationType;

        /// Variable related to interactor behaviour (0) 
        int DisplayInterpolationProcess;

        /// *Maximum number of refinement iterations, set to 1 if only one refinement is wanted (10) 
        int MaxNumberOfRefinementIterations;

        /// *Averaging method : 1 = simple averaging, 2 = closest point, 3 = max, 4 = min, 5 = inverse weighted distance, 6 = minabs, 7 = kdtree (1)
        int AveragingMethod;

        /// *Minimum number of points needed inside cell to handle the cell (1)
        int MinimumNumberOfPoints;

        /// *Relative search cell size, default 1= actual cell size, 2= twice as large, search radius can be larger than cell so more sample are included. (1.01)
        double RelativeSearchRadius;

        /// *Interpolation settings, 1=bathy, 2=zk, 3=s1, 4=Zc (2)
        int InterpolateTo;
    };
}


