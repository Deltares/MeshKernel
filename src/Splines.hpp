#pragma once

#include <vector>
#include "Operations.cpp"
#include "Entities.hpp"

namespace GridGeom
{
    class Splines
    {
    public:

        Splines() : m_numAllocatedSplines(0), m_numSplines(0)
        {
            AllocateVector(m_numAllocatedSplines,
                m_splines,
                std::vector<Point>(10, { doubleMissingValue, doubleMissingValue }),
                5);

            m_numAllocatedSplineNodes.resize(m_numAllocatedSplines, 10);
            m_numSplineNodes.resize(m_numAllocatedSplines, 0);
        }

        /// add a new spline
        bool Set(const std::vector<Point>& splines)
        {
            AllocateVector(m_numAllocatedSplines, m_splines, std::vector<Point>(10, { doubleMissingValue, doubleMissingValue }),5);

            m_numAllocatedSplineNodes.resize(m_numAllocatedSplines, 10);
            m_numSplineNodes.resize(m_numAllocatedSplines, 0);

            m_splines[m_numSplines] = splines;
            m_numSplineNodes[m_numSplines] = splines.size();
            m_numSplines++;
            return true;
        }

        /// add a new spline point in an existing spline
        bool Set(int splineIndex, const Point& point)
        {
            if (splineIndex >= m_numSplines)
            {
                return false;
            }
            AllocateVector(m_numAllocatedSplineNodes[splineIndex], m_splines[splineIndex],{ doubleMissingValue, doubleMissingValue },10);

            m_splines[splineIndex][m_numSplineNodes[splineIndex]] = point;
            m_numSplineNodes[splineIndex]++;
            return true;
        }

        /// splint
        static bool Interpolate(const std::vector<Point>& coordinates, const std::vector<Point>& coordinatesDerivatives, double pointAdimensionalCoordinate, Point& pointCoordinate)
        {

            const double eps = 1e-5;
            const double splFac = 1.0;
            int intCoordinate = std::floor(pointAdimensionalCoordinate);
            if (pointAdimensionalCoordinate - intCoordinate < eps)
            {
                pointCoordinate = coordinates[intCoordinate];
                return true;
            }

            int low = intCoordinate;
            int high = low + 1;
            double a = high - pointAdimensionalCoordinate;
            double b = pointAdimensionalCoordinate - low;
            
            pointCoordinate = coordinates[low] * a  + coordinates[high] * b +
                (coordinatesDerivatives[low] * (pow(a, 3) - a) + coordinatesDerivatives[high] * (pow(b, 3) - b)) / 6.0 * splFac ;

            return true;
        }

        /// SPLINE
        static bool Derivative(const std::vector<Point>& coordinates, std::vector<Point>& coordinatesDerivatives)
        {
            std::vector<Point> u(coordinates.size());
            u[0] = { 0.0, 0.0 };
            coordinatesDerivatives[0] = { 0.0, 0.0 };

            for (int i = 1; i < coordinates.size() - 1; i++)
            {
                const Point p =  coordinatesDerivatives[i - 1] * 0.5 + 2.0;
                coordinatesDerivatives[i].x = -0.5 / p.x;
                coordinatesDerivatives[i].y = -0.5 / p.y;

                const Point delta = coordinates[i + 1] - coordinates[i] - (coordinates[i] - coordinates[i - 1]);
                u[i] = (delta *6.0 / 2.0 - u[i - 1] * 0.5 ) / p;                
            }

            coordinatesDerivatives[coordinates.size() - 1] = { 0.0, 0.0 };
            for (int i = coordinates.size() - 2; i >= 0; i--)
            {
                coordinatesDerivatives[i] = coordinatesDerivatives[i] * coordinatesDerivatives[i + 1] + u[i];
            }

            return true;
        }

        int m_numSplines;
        int m_numAllocatedSplines;
        std::vector<int> m_numAllocatedSplineNodes;
        std::vector<int> m_numSplineNodes;
        std::vector<std::vector<Point>> m_splines;
    };

}