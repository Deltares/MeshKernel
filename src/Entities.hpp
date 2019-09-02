#ifndef ENTITIES_HPP
#define ENTITIES_HPP

#include <utility>
#include <vector>
#include <Eigen/Dense>

#define _USE_MATH_DEFINES

namespace GridGeom
{
    struct cartesianPoint
    {
        double x;
        double y;
    };

    struct sphericalPoint
    {
        double x;
        double y;
    };

    struct cartesian3DPoint
    {
        double x;
        double y;
        double z;
    };

    struct Nodes
    {
        std::vector<double> x;
        std::vector<double> y;
    };

    typedef std::pair<size_t, size_t> Edge;

    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorED;
    typedef Eigen::Matrix<int, Eigen::Dynamic, 1> VectorEI;

}

#endif
