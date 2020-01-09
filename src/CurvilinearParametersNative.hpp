#pragma once

namespace GridGeomApi
{
    struct CurvilinearParametersNative
    {
        /// <summary>
        /// *M-refinement factor for regular grid generation (2000) 
        /// </summary>
        int MRefinement;

        /// <summary>
        /// *N-refinement factor for regular grid generation (40) 
        /// </summary>
        int NRefinement;

        /// <summary>
        /// Nr. of inner iterations in regular grid smoothing (10).
        /// </summary>
        int SmoothingIterations;

        /// <summary>
        /// * Smoothing parameter (0.5).
        /// </summary>
        double SmoothingParameter;

        /// <summary>
        /// *Attraction/repulsion parameter (0).
        /// </summary>
        double AttractionParameter;
    };
}


