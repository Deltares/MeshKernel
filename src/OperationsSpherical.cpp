#pragma once

#include <cmath>
#include "Entities.hpp"
#include "Constants.cpp"
#include "Operations.cpp"
#include "IOperations.hpp"

namespace GridGeom
{
    // spherical point
    struct OperationsSpherical : IOperations
    {
        void normalVector(const Point& firstPoint, const Point& secondPoint, const Point& insidePoint, Point& result) override
        {
            cartesian3DPoint firstPointCartesianCoordinates;
            cartesian3DPoint secondPointCartesianCoordinates;
            sphericalToCartesian(firstPoint, firstPointCartesianCoordinates);
            sphericalToCartesian(secondPoint, secondPointCartesianCoordinates);

            double lambda = insidePoint.x * degrad_hp;
            double phi = insidePoint.y * degrad_hp;
            double elambda[3] = { -sin(lambda), cos(lambda), 0.0 };
            double ephi[3] = { -sin(lambda), cos(lambda), 0.0 };

            double dx = (secondPointCartesianCoordinates.x - firstPointCartesianCoordinates.x) * elambda[0] +
                (secondPointCartesianCoordinates.y - firstPointCartesianCoordinates.y) * elambda[1] +
                (secondPointCartesianCoordinates.z - firstPointCartesianCoordinates.z) * elambda[2];

            double dy = (secondPointCartesianCoordinates.x - firstPointCartesianCoordinates.x) * ephi[0] +
                (secondPointCartesianCoordinates.y - firstPointCartesianCoordinates.y) * ephi[1] +
                (secondPointCartesianCoordinates.z - firstPointCartesianCoordinates.z) * ephi[2];

            double squaredDistance = dx * dx + dy * dy;
            if (squaredDistance != 0.0)
            {
                const double distance = sqrt(squaredDistance);
                result.x = dx / distance;
                result.y = dy / distance;
            }
        }

        void normalVectorInside(const Point& firstPoint, const Point& secondPoint, const Point& insidePoint, Point& result, bool& flippedNormal) override
        {

            //if (JSFERIC.eq.1 . and .jasfer3D.eq.0) xn = xn * cos(dg2rd * 0.5d0 * (y0 + y1)) !normal vector needs to be in Cartesian coordinates
        }

        double getDx(const Point& firstPoint, const Point& secondPoint) override
        {
            double  firstPointYDiff = abs(abs(firstPoint.y) - 90.0);
            double  secondPointYDiff = abs(abs(secondPoint.y) - 90.0);
            if (firstPointYDiff <= dtol_pole && secondPointYDiff > dtol_pole || firstPointYDiff > dtol_pole && secondPointYDiff <= dtol_pole)
            {
                return 0.0;
            }
            double firstPointX = firstPoint.x;
            double secondPointX = secondPoint.x;
            if (firstPointX - secondPointX > 180.0)
            {
                firstPointX -= 360.0;
            }
            else if (firstPointX - secondPointX < -180.0)
            {
                firstPointX += 360.0;
            }

            firstPointX = firstPointX * degrad_hp;
            secondPointX = secondPointX * degrad_hp;
            double firstPointY = firstPoint.y * degrad_hp;
            double secondPointY = secondPoint.y * degrad_hp;
            double cosPhi = cos(0.5 * (firstPointY + secondPointY));
            double dx = earth_radius * cosPhi * (secondPointX - firstPointX);
            return dx;
        }

        double getDy(const Point& firstPoint, const Point& secondPoint) override
        {
            double firstPointY = firstPoint.y * degrad_hp;
            double secondPointY = secondPoint.y * degrad_hp;
            double dy = earth_radius * (secondPointY - firstPointY);
            return dy;
        }

        void add(Point& point, const Point& normal, const double increment) override
        {
            double convertedIncrement = raddeg_hp * increment / earth_radius;
            double xf = 1.0 / cos(degrad_hp * point.y);
            point.x = point.x + normal.x * convertedIncrement * xf;
            point.y = point.y + normal.y * convertedIncrement;
        }

        void referencePoint(std::vector<Point>& polygon, const int numPoints, double& minX, double& minY) override
        {
            minX = std::numeric_limits<double>::max();
            minY = std::numeric_limits<double>::max();
            double xmax = std::numeric_limits<double>::min();
            for (int i = 0; i < numPoints; i++)
            {
                if (polygon[i].x < minX)
                {
                    minX = polygon[i].x;
                }
                if (abs(polygon[i].y) < abs(minY))
                {
                    minY = polygon[i].y;
                }
                if (polygon[i].x > xmax)
                {
                    xmax = polygon[i].x;
                }
            }

            if (xmax - minX > 180.0)
            {
                double deltaX = xmax - 180.0;
                for (int i = 0; i < numPoints; i++)
                {
                    if (polygon[i].x < deltaX)
                    {
                        polygon[i].x = polygon[i].x + 360.0;
                    }
                }
                minX = minX + 360.0;
            }
            //TODO: check result
            minX = std::min_element(polygon.begin(), polygon.end(), [](const Point& p1, const Point& p2) { return p1.x < p2.x; })->x;
        }

        double distance(const Point& firstPoint, const Point& secondPoint) override
        {
            return -1.0;
        }

        //dLINEDIS3
        double distanceFromLine(const Point& p3, const Point& p1, const Point& p2, Point& normalPoint, double& rlout) override
        {
            //TODO: implement me
            return -1.0;
        }


        //out product of two segments
        double outerProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment, const Point& firstPointSecondSegment, const Point& secondPointSecondSegment) override
        {
            //TODO: IMPLEMENTATION IS MISSING

            return 0.0;
        }

        double innerProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment, const Point& firstPointSecondSegment, const Point& secondPointSecondSegment) override
        {
            cartesian3DPoint firstPointFirstSegment3D;
            cartesian3DPoint secondPointFirstSegment3D;
            cartesian3DPoint firstPointSecondSegment3D;
            cartesian3DPoint secondPointSecondSegment3D;

            sphericalToCartesian(firstPointFirstSegment, firstPointFirstSegment3D);
            sphericalToCartesian(secondPointFirstSegment, secondPointFirstSegment3D);
            sphericalToCartesian(firstPointSecondSegment, firstPointSecondSegment3D);
            sphericalToCartesian(secondPointSecondSegment, secondPointSecondSegment3D);

            double dx1 = secondPointFirstSegment3D.x - firstPointFirstSegment3D.x;
            double dy1 = secondPointFirstSegment3D.y - firstPointFirstSegment3D.y;
            double dz1 = secondPointFirstSegment3D.z - firstPointFirstSegment3D.z;

            double dx2 = secondPointSecondSegment3D.x - firstPointSecondSegment3D.x;
            double dy2 = secondPointSecondSegment3D.y - firstPointSecondSegment3D.y;
            double dz2 = secondPointSecondSegment3D.z - firstPointSecondSegment3D.z;

            return dotProduct(dx1, dx2, dy1, dy2, dz1, dz2);
        }

        //TODO:: comp_local_coords
        bool orthogonalizationComputeLocalCoordinates(const std::vector<size_t>& m_nodesNumEdges, const std::vector<size_t>& numConnectedNodes, std::vector<int>& localCoordinates) override
        {
            localCoordinates.resize(m_nodesNumEdges.size(), 0);
            localCoordinates[0] = 1;
            for (int i = 0; i < m_nodesNumEdges.size(); i++)
            {
                localCoordinates[i + 1] = localCoordinates[i] + std::max(m_nodesNumEdges[i] + 1, numConnectedNodes[i]);
            }
            return true;
        }


        bool orthogonalizationComputeJacobian(const int currentNode, const std::vector<double>& Jxi, const std::vector<double>& Jeta, const std::vector<size_t>& connectedNodes, const int numNodes, const std::vector<Point>& nodes, std::vector<double>& J) override
        {
            double factor = std::cos(nodes[currentNode].y) * degrad_hp;
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
            double wwxTransformed = wwx * earth_radius * degrad_hp;
            double wwyTransformed = wwy * earth_radius * degrad_hp;

            increments[0] += wwxTransformed;
            increments[1] += wwyTransformed;

            dx0 = dx0 + wwxTransformed * (nodes[firstNode].x - nodes[secondNode].x);
            dy0 = dy0 + wwxTransformed * (nodes[firstNode].y - nodes[secondNode].y);

            return true;
        }

        bool orthogonalizationComputeCoordinates(double dx0, double dy0, const Point& point, Point& updatedPoint) override
        {
            //TODO: implement
            //if (jsferic.eq.1 . and .jasfer3D.eq.1) then
            //    dumx(1) = relaxin * Dx0
            //    dumy(1) = relaxin * Dy0
            //    call loc2spher(xk(k), yk(k), 1, dumx, dumy, xk1(k), yk1(k))
            //else

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
            if (den > 1e-16)
            {
                z = (dx2 * (dx2 - dx3) + dy2 * (dy2 - dy3)) / den;
            }

            //TODO circumcenter3 FINISH
            //phi = (y(1) + y(2) + y(3)) / 3d0
            //    xf = 1d0 / dcos(degrad_hp * phi)
            //    xz = x(1) + xf * 0.5d0 * (dx3 - z * dy3) * raddeg_hp / earth_radius
            //    yz = y(1) + 0.5d0 * (dy3 + z * dx3) * raddeg_hp / earth_radius
            return true;
        }
    };
}