#ifndef ENTITIES_CPP
#define ENTITIES_CPP

#include <utility>
#define _USE_MATH_DEFINES

namespace GridGeom
{
    struct Point
    {
        double x;
        double y;
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
