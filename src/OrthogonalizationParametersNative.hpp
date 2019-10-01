#pragma once

namespace GridGeomApi
{
struct OrthogonalizationParametersNative
{
    int OuterIterations;
    int BoundaryIterations;
    int InnerIterations;
    double OrthogonalizationToSmoothingFactor;
    double AtpfB;
    double Circumormasscenter;
    double Smoothorarea;
    int AdaptMethod;
    double AdaptBeta;
    int AdaptNiterU;
    int AdaptNiterG;
    double OrthoPure;
};
}


