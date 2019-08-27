#ifndef ENTITIES_CPP
#define ENTITIES_CPP

#include <utility>
#include <vector>
#include <Eigen/Dense>
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

    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorED;
    typedef Eigen::Matrix<int, Eigen::Dynamic, 1> VectorEI;
}

#endif
