#ifndef CONSTANTS_CPP
#define CONSTANTS_CPP

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

namespace GridGeom
{
    // missing value
    static constexpr double doubleMissingValue = -999.0;
    static constexpr int    intMissingValue  = -999;

    //geometric constants
    static constexpr double degrad_hp = M_PI / 180.0; // conversion factor from degrees to radians(pi / 180)
    static constexpr double raddeg_hp = 180.0 / M_PI; // conversion factor from radians to degrees(180 / pi)
    static constexpr double earth_radius = 6378137.0; // earth radius(m)
    static constexpr double dtol_pole = 0.0001;       // pole tolerance in degrees

    //mesh constants
    static constexpr double minimumDeltaCoordinate = 1e-14;
    static constexpr size_t maximumNumberOfEdgesPerNode = 10;
    static constexpr size_t maximumNumberOfNodesPerFace = 8;
    static constexpr size_t maximumNumberOfConnectedNodes = maximumNumberOfEdgesPerNode * 4;
    static constexpr double minimumCellArea = 1e-12;

    //orthogonalization 
    static constexpr double minimumEdgeLength = 1e-4;
    static constexpr double curvilinearToOrthogonalRatio= 0.5; //curvi - linear - like(0.0) or pure(1.0) orthogonalisation
    static constexpr double orthogonalizationToSmoothingFactor = 0.975; //Factor between grid smoothing and grid ortho resp (0.<=ATPF<=1.)
    static constexpr int orthogonalizationOuterIterations = 2; // Nr.of outer    iterations in grid / net orthogonalisation.
    static constexpr int orthogonalizationBoundaryIterations = 25; // Nr.of boundary iterations in grid / net orthogonalisation. (within ITBND)
    static constexpr int orthogonalizationInnerIterations = 25;    // Nr.of inner    iterations in grid / net orthogonalisation. (within ITIN)
    static constexpr double relaxationFactorOrthogonalizationUpdate = 0.75;
}

#endif