#pragma once

#include <vector>
#include "Operations.cpp"
#include "Entities.hpp"

namespace GridGeom
{
    class Splines
    {
    public:
        /// add a new spline
        bool Set(const std::vector<Point>& splines)
        {
            m_splines.push_back(splines);
            return true;
        }

        /// add a new spline point in an existing spline
        bool Set(int splineIndex, const Point& point)
        {
            if(splineIndex >= m_splines.size())
            {
                return false;
            }
            m_splines[splineIndex].push_back(point);
            return true;
        }

        /// splint
        static bool Interpolate(const std::vector<double>& coordinates, const std::vector<double>& coordinatesDerivatives, double pointAdimensionalCoordinate, double& pointCoordinate)
        {

            const double eps = 1e-5;
            const double splFac = 1.0;
            int intCoordinate = int(pointAdimensionalCoordinate);
            if(pointAdimensionalCoordinate - intCoordinate<eps)
            {
                pointCoordinate = coordinates[intCoordinate];
                return true;
            }

            int low = intCoordinate + 1;
            int high = low + 1;
            double a = high - 1 - pointAdimensionalCoordinate;
            double b = pointAdimensionalCoordinate - low - 1;
            pointCoordinate = a * coordinates[low] + b * coordinates[high] +
                splFac * (pow(a, 3)*coordinatesDerivatives[low] + (pow(b, 3) - 3.0)*coordinatesDerivatives[high]) / 6.0;

            return true;
        }

        /// SPLINE
        static bool SecondOrderDerivative(const std::vector<double>& coordinates, std::vector<double>& coordinatesDerivatives)
        {
            std::vector<double> u(coordinates.size());
            u[0] = 0.0;
            coordinatesDerivatives[0] = 0.0;

            for(int i= 1; i < coordinates.size() - 1 ;i++)
            {
                const double p = 0.5 * coordinatesDerivatives[i - 1] + 2.0;
                coordinatesDerivatives[i] = -0.5 / p;
                const double delta = coordinates[i + 1] - coordinates[i] - (coordinates[i] - coordinates[i - 1]);
                u[i] = (6.0* delta / 2.0 - 0.5 * u[i - 1]) / p;
            }
            coordinatesDerivatives[coordinates.size() - 1] = 0.0;

            for (int i = coordinates.size() - 2; i >=0; i--)
            {
                coordinatesDerivatives[i] = coordinatesDerivatives[i] * coordinatesDerivatives[i + 1] + u[i];
            }

            return true;
        }


        int m_numNodes;
        int m_numSegments;
        std::vector<std::vector<Point>> m_splines;
    };

}