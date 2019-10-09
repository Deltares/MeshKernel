#pragma once

#include "Entities.hpp"

namespace GridGeom
{
    struct IOperations
    {
        //normalout, Creates the relative unit normal vector to edge 1->2
        virtual void normalVector(const Point& firstPoint, const Point& secondPoint, const Point& insidePoint, Point& result) = 0;

        //normaloutchk
        virtual void normalVectorInside(const Point& firstPoint, const Point& secondPoint, const Point& insidePoint, Point& result, bool& flippedNormal) =0;

        virtual double getDx(const Point& firstPoint, const Point& secondPoint) = 0;

        virtual double getDy(const Point& firstPoint, const Point& secondPoint) = 0;

        virtual void add(Point& point, const Point& normal, const double increment) = 0;

        virtual void referencePoint(std::vector<Point>& polygon, const int numPoints, double& minX, double& minY) = 0;

        //dbdistance
        virtual double distance(const Point& firstPoint, const Point& secondPoint) =0;

        //dLINEDIS3
        virtual double distanceFromLine(const Point& p3, const Point& p1, const Point& p2, Point& normalPoint, double& rlout) = 0;

        //out product of two segments
        virtual double outerProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment, const Point& firstPointSecondSegment, const Point& secondPointSecondSegment)=0;

        //inner product of two segments
        virtual double innerProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment, const Point& firstPointSecondSegment, const Point& secondPointSecondSegment)=0;

        virtual bool orthogonalizationComputeLocalCoordinates(const std::vector<size_t>& m_nodesNumEdges, const std::vector<size_t>& numConnectedNodes, std::vector<int>& localCoordinates)=0;

        virtual bool orthogonalizationComputeJacobian(const int currentNode, const std::vector<double>& Jxi, const std::vector<double>& Jeta, const std::vector<size_t>& connectedNodes, const int numNodes, const std::vector<Point>& nodes, std::vector<double>& J) = 0;

        virtual bool orthogonalizationComputeDeltas(int firstNode, int secondNode, double wwx, double wwy, const std::vector<Point>& nodes, double& dx0, double& dy0, std::vector<double>& increments) = 0;

        virtual bool orthogonalizationComputeCoordinates(double dx0, double dy0, const Point& point, Point& updatedPoint) = 0;

        virtual bool circumcenterOfTriangle(const Point& p1, const Point& p2, const Point& p3, Point& circumcenter) = 0;

    };
}