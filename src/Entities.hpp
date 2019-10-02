#pragma once

#include <utility>
#include <vector>
#include <Eigen/Dense>

#define _USE_MATH_DEFINES

namespace GridGeom
{
    enum OperationTypes 
    {
        cartesianOperations,
        sphericalOperations
    };

    // functions that depends of the point type
    template<OperationTypes operationType>
    struct Operations;

    struct Point
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

    struct MeshBase
    {
        //getters..
    };

    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorED;
    typedef Eigen::Matrix<int, Eigen::Dynamic, 1> VectorEI;

}