#pragma once

#include <cmath>
#include "Entities.hpp"
#include "Constants.cpp"
#include "IOperations.hpp"
#include <algorithm>

namespace GridGeom
{
    // cartesian points
    struct OperationsCartesian: IOperations
    {
        OperationsCartesian() = default;

        //normalout, Creates the relative unit normal vector to edge 1->2
        void normalVector(const Point& firstPoint, const Point& secondPoint, const Point& insidePoint, Point& result) override
        {
            double dx = getDx(firstPoint, secondPoint);
            double dy = getDy(firstPoint, secondPoint);
            const double squaredDistance = dx * dx + dy * dy;
            if (squaredDistance != 0.0)
            {
                const double distance = sqrt(squaredDistance);
                result.x = dy / distance;
                result.y = -dx / distance;
            }
        }

        //normaloutchk
        void normalVectorInside(const Point& firstPoint, const Point& secondPoint, const Point& insidePoint, Point& result, bool& flippedNormal) override
        {
            normalVector(firstPoint, secondPoint, insidePoint, result);
            flippedNormal = false;
            Point thirdPoint{ firstPoint.x + result.x, firstPoint.y + result.y };

            if (outerProductTwoSegments(firstPoint, thirdPoint, firstPoint, secondPoint) * outerProductTwoSegments(firstPoint, insidePoint, firstPoint, secondPoint) > 0.0)
            {
                result.x = -result.x;
                result.y = -result.y;
                flippedNormal = true;
            }
            else
            {
                flippedNormal = false;
            }
        }

        double getDx(const Point& firstPoint, const Point& secondPoint) override
        {
            return secondPoint.x - firstPoint.x;
        }

        double getDy(const Point& firstPoint, const Point& secondPoint) override
        {
            return secondPoint.y - firstPoint.y;
        }

        void add(Point& point, const Point& normal, const double increment) override
        {
            point.x = point.x + normal.x * increment;
            point.y = point.y + normal.y * increment;
        }

        void referencePoint(std::vector<Point>& polygon, const int numPoints, double& minX, double& minY) override
        {
            minX = std::numeric_limits<double>::max();
            minY = std::numeric_limits<double>::max();
            for (int i=0; i< numPoints;i++)
            {
                if (polygon[i].x < minX)
                {
                    minX = polygon[i].x;
                }
                if (abs(polygon[i].y) < abs(minY))
                {
                    minY = polygon[i].y;
                }
            }
        }

        //dbdistance
        double distance(const Point& firstPoint, const Point& secondPoint) override
        {
            double dx = getDx(firstPoint, secondPoint);
            double dy = getDy(firstPoint, secondPoint);
            const double squaredDistance = dx * dx + dy * dy;
            double distance = 0.0;
            if (squaredDistance != 0.0)
            {
                distance = sqrt(squaredDistance);
            }
            return distance;
        }

        //dLINEDIS3
        double distanceFromLine(const Point& p3, const Point& p1, const Point& p2, Point& normalPoint, double& rlout) override
        {
            double dis = 0.0;
            double r2 = distance(p2, p1);
            if (r2 != 0.0)
            {
                double rl = (getDx(p1, p3) * getDx(p1, p2) + getDy(p1, p3) * getDy(p1, p2)) / (r2 * r2);
                rlout = std::max(std::min(1.0, rl), 0.0);
                normalPoint.x = p1.x + rlout * (p2.x - p1.x);
                normalPoint.y = p1.y + rlout * (p2.y - p1.y);
                dis = distance(p3, normalPoint);
            }
            return dis;
        }

        //out product of two segments
        double outerProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment, const Point& firstPointSecondSegment, const Point& secondPointSecondSegment) override
        {
            double dx1 = getDx(firstPointFirstSegment, secondPointFirstSegment);
            double dx2 = getDx(firstPointSecondSegment, secondPointSecondSegment);

            double dy1 = getDy(firstPointFirstSegment, secondPointFirstSegment);
            double dy2 = getDy(firstPointSecondSegment, secondPointSecondSegment);
         
            return dx1 * dy2 - dy1 * dx2;
        }

        //inner product of two segments
        double innerProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment, const Point& firstPointSecondSegment, const Point& secondPointSecondSegment) override
        {
            double dx1 = getDx(firstPointFirstSegment, secondPointFirstSegment);
            double dx2 = getDx(firstPointSecondSegment, secondPointSecondSegment);

            double dy1 = getDy(firstPointFirstSegment, secondPointFirstSegment);
            double dy2 = getDy(firstPointSecondSegment, secondPointSecondSegment);

            return dx1 * dx2 + dy1 * dy2;
        }


        bool orthogonalizationComputeLocalCoordinates(const std::vector<size_t>& m_nodesNumEdges, const std::vector<size_t>& numConnectedNodes, std::vector<int>& localCoordinates) override
        {
            //do nothing
            return true;
        }


        bool orthogonalizationComputeJacobian(const int currentNode, const std::vector<double>& Jxi, const std::vector<double>& Jeta, const std::vector<size_t>& connectedNodes, const int numNodes, const std::vector<Point>& nodes, std::vector<double>& J) override
        {
            J[0] = 0.0; 
            J[1] = 0.0;
            J[2] = 0.0;
            J[3] = 0.0;
            for (int i = 0; i < numNodes; i++)
            {
                J[0] += Jxi[i] * nodes[connectedNodes[i]].x;
                J[1] += Jxi[i] * nodes[connectedNodes[i]].y;
                J[2] += Jeta[i] * nodes[connectedNodes[i]].x;
                J[3] += Jeta[i] * nodes[connectedNodes[i]].y;
            }
            return true;
        }

        bool orthogonalizationComputeDeltas(int firstNode, int secondNode, double wwx, double wwy, const std::vector<Point>& nodes, double& dx0, double& dy0, std::vector<double>& increments) override
        {

            increments[0] += wwx;
            increments[1] += wwy;

            dx0 = dx0 + wwx * (nodes[firstNode].x - nodes[secondNode].x);
            dy0 = dy0 + wwy * (nodes[firstNode].y - nodes[secondNode].y);
            return true;
        }

        bool orthogonalizationComputeCoordinates(double dx0, double dy0, const Point& point, Point& updatedPoint) override
        {
            double x0 = point.x + dx0;
            double y0 = point.y + dy0;
            static constexpr double relaxationFactorCoordinates = 1.0 - relaxationFactorOrthogonalizationUpdate;

            updatedPoint.x = relaxationFactorOrthogonalizationUpdate * x0 + relaxationFactorCoordinates * point.x;
            updatedPoint.y = relaxationFactorOrthogonalizationUpdate * y0 + relaxationFactorCoordinates * point.y;

            return true;
        }

        bool circumcenterOfTriangle(const Point& p1, const Point& p2, const Point& p3, Point& circumcenter) override
        {
            double dx2 = getDx(p1, p2);
            double dy2 = getDy(p1, p2);

            double dx3 = getDx(p1, p3);
            double dy3 = getDy(p1, p3);

            double den = dy2 * dx3 - dy3 * dx2;
            double z = 0.0;
            if (den != 0.0)
            {
                z = (dx2 * (dx2 - dx3) + dy2 * (dy2 - dy3)) / den;
            }

            circumcenter.x = p1.x + 0.5 * (dx3 - z * dy3);
            circumcenter.y = p1.y + 0.5 * (dy3 + z * dx3);
            return true;
        }

        bool linesCrossing(const Point& firstSegmentFistPoint, const Point& firstSegmentSecondPoint, const Point& secondSegmentFistPoint, const Point& secondSegmentSecondPoint, bool adimensional, Point& intersection, double& crossProduct) override
        {
            bool isCrossing = false;

            double x21 = getDx(firstSegmentFistPoint, firstSegmentSecondPoint);
            double y21 = getDy(firstSegmentFistPoint, firstSegmentSecondPoint);

            double x43 = getDx(secondSegmentFistPoint, secondSegmentSecondPoint);
            double y43 = getDy(secondSegmentFistPoint, secondSegmentSecondPoint);

            double x31 = getDx(firstSegmentFistPoint, secondSegmentFistPoint);
            double y31 = getDy(firstSegmentFistPoint, secondSegmentFistPoint);

            double det = x43 * y21 - y43 * x21;

            std::vector<double> values{ x21, y21, x43, y43 };
            double eps = std::max(0.00001 * (*std::max_element(values.begin(), values.end())), std::numeric_limits<double>::denorm_min());

            if (det < eps)
            {
                return isCrossing;
            }

            double sm = (y31 * y21 - x31 * y21) / det;
            double sl = (y31 * x43 - x31 * y43) / det;
            if (sm >= 0.0 && sm <= 1.0 && sl >= 0.0 && sl <= 1.0)
            {
                isCrossing = true;
            }
            intersection.x = firstSegmentFistPoint.x + sl * (firstSegmentSecondPoint.x - firstSegmentFistPoint.x);
            intersection.y = firstSegmentFistPoint.y + sl * (firstSegmentSecondPoint.x - firstSegmentFistPoint.y);
            crossProduct = -det;
            if (adimensional)
            {
                crossProduct = -det / (std::sqrt(x21 * x21 + y21 * y21) * std::sqrt(x43 * x43 + y43 * y43) + 1e-8);
            }
            return isCrossing;
        }

    };
}
