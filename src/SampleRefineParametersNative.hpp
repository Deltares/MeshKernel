#pragma once

namespace GridGeomApi
{
    struct SampleRefineParametersNative
    {
        /// Sample vector dimension (1) 
        int SampleVectorDimension;

        /// *Maximum number of refinement iterations, set to 1 if only one refinement is wanted (10) 
        int MaxNumberOfRefinementIterations;

        /// *Minimum cell size (50000.0)
        double MinimumCellSize;

        /// *Directional refinement, 1 yes 0 no (0)
        int DirectionalRefinement;

        /// *Refinement criterion type (2)
        int RefinementType;

        /// *Connect hanging nodes at the end of the iteration, 1 yes 0 no (1) 
        int ConnectHangingNodes;

        /// *Maximum time-step in courant grid (0.0)
        double MaximumTimeStepInCourantGrid;

        /// *Take samples outside cell into account , 1 yes 0 no (1)
        int AccountForSamplesOutside;
    };
}


