#pragma once

namespace GridGeomApi
{
    struct OrthogonalizationParametersNative
    {
        int OuterIterations;
        int BoundaryIterations;
        int InnerIterations;
        double OrthogonalizationToSmoothingFactor;
        double OrthogonalizationToSmoothingFactorBoundary;
        double Circumormasscenter;
        double Smoothorarea;
        int AdaptMethod;
        double AdaptBeta;
        int AdaptNiterU;
        int AdaptNiterG;
        double OrthoPure;
    };
}


