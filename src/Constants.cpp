#pragma once

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif
#include <math.h> 

namespace GridGeom
{
    // missing value
    static constexpr double doubleMissingValue = -999.0;
    static constexpr int    intMissingValue    = -999;
    static constexpr double squareRootOfThree  = 1.73205080756887729352;
    static constexpr double oneThird           = 0.33333333333333333333;

    //geometric constants
    static constexpr double degrad_hp = M_PI / 180.0; // conversion factor from degrees to radians(pi / 180)
    static constexpr double raddeg_hp = 180.0 / M_PI; // conversion factor from radians to degrees(180 / pi)
    static constexpr double earth_radius = 6378137.0; // earth radius(m)
    static constexpr double dtol_pole = 0.0001;       // pole tolerance in degrees

    //mesh constants
    static constexpr double minimumDeltaCoordinate = 1e-14;
    static constexpr int maximumNumberOfEdgesPerNode = 10;
    static constexpr int maximumNumberOfEdgesPerFace = 6;
    static constexpr int maximumNumberOfNodesPerFace = 8;
    static constexpr int maximumNumberOfConnectedNodes = maximumNumberOfEdgesPerNode * 4;
    static constexpr double minimumCellArea = 1e-12;

    //orthogonalization 
    static constexpr double minimumEdgeLength = 1e-4;
    static constexpr double curvilinearToOrthogonalRatio= 0.5; //curvi - linear - like(0.0) or pure(1.0) orthogonalisation
    static constexpr double orthogonalizationToSmoothingFactor = 0.975; //Factor between grid smoothing and grid ortho resp (0.<=ATPF<=1.)
    static constexpr double relaxationFactorOrthogonalizationUpdate = 0.75;

    // merging distance
    static constexpr double mergingDistance = 0.001;

    // physical constants
    static constexpr double gravity = 9.81;

    // Operations averaging methods
    enum AveragingMethod 
    {
        SimpleAveraging = 1,
        ClosestPoint= 2,
        Max = 3,
        Min = 4,
        InverseWeightDistance = 5,
        MinAbs = 6,
        KdTree = 7
    };

}