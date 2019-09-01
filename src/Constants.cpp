#ifndef CONSTANTS_CPP
#define CONSTANTS_CPP

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

namespace GridGeom
{
    // missing value
    static constexpr double doubleMissingValue = -999.0;

    //geometric constants
    static constexpr double degrad_hp = M_PI / 180.0; // conversion factor from degrees to radians(pi / 180)
    static constexpr double earth_radius = 6378137.0; // earth radius(m)
    static constexpr double dtol_pole = 0.0001;       // pole tolerance in degrees
    static constexpr double raddeg_hp = 180.0 / M_PI; // conversion factor from radians to degrees(180 / pi)

    //mesh constants
    static constexpr double minimumDeltaCoordinate = 1e-14;
    static constexpr size_t maximumNumberOfEdgesPerNode = 8;
    static constexpr size_t maximumNumberOfNodesPerFace = 10;

    //orthogonalization 
    static constexpr double minimumEdgeLength = 1e-4;
    static constexpr double curvilinearToOrthogonalRatio= 0.5; //curvi - linear - like(0.0) or pure(1.0) orthogonalisation

}

#endif