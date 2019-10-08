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

    // contains a store of static functions
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
        virtual bool setMesh(const std::vector<GridGeom::Edge>& edges, const std::vector<Point>& nodes) = 0;
        virtual bool setState() = 0;
        virtual int getNumFaces() = 0;
        //Used for internal state
        std::vector<double> m_nodex;
        std::vector<double> m_nodey;
        std::vector<double> m_nodez;
        std::vector<int>    m_edgeNodes;
    };

    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorED;
    typedef Eigen::Matrix<int, Eigen::Dynamic, 1> VectorEI;

}