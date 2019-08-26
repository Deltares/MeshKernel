#ifndef ENTITIES_CPP
#define ENTITIES_CPP

#include <utility>
#include <vector>
#define _USE_MATH_DEFINES

namespace GridGeom
{
    struct Point
    {
        double x;
        double y;
    };

    struct Points
    {
        std::vector<double> x;
        std::vector<double> y;
    };

    struct Point3D
    {
        double x;
        double y;
        double z;
    };

    typedef std::pair<size_t, size_t> Edge;
}

#endif
