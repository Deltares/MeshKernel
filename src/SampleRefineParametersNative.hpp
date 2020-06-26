#pragma once

namespace GridGeomApi
{
    struct SampleRefineParametersNative
    {
        /// Sample vector dimension
        int SampleVectorDimension;

        /// Maximum number of refinement iterations, set to 1 if only one refinement is wanted
        int MaxNumberOfRefinementIterations;

        /// Minimum cell size
        double MinimumCellSize;

        /// Directional refinement, 1 yes 0 no
        int DirectionalRefinement;

        /// Refinement criterion type
        int RefinementType;

        /// Connect hanging nodes at the end of the iteration, 1 yes 0 no
        int ConnectHangingNodes;

        /// Maximum time-step in courant grid
        double MaximumTimeStepInCourantGrid;

        /// Take samples outside face into account , 1 yes 0 no
        int AccountForSamplesOutside;

    };
}


