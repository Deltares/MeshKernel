#pragma once

#include <vector>
#include "Entities.hpp"
#include "Polygons.hpp"
#include "CurvilinearParametersNative.hpp"
#include "SplinesToCurvilinearParametersNative.hpp"

namespace GridGeom
{
    class CurvilinearGrid;
    
    class Splines
    {

    public:

        /// <summary>
        /// Ctor
        /// </summary>
        /// <returns></returns>
        Splines();

        /// <summary>
        /// Ctor and set projection
        /// </summary>
        /// <param name="projection">The map projection</param>
        /// <returns></returns>
        Splines(Projections projection);

        /// <summary>
        /// Adds a new spline to m_splineCornerPoints 
        /// </summary>
        /// <param name="splines">A vector containing the spline corner points</param>
        /// <param name="start">The index in splines of the starting node to add</param>
        /// <param name="size">The index in splines of the last node to add</param>
        /// <returns></returns>
        bool AddSpline(const std::vector<Point>& splines, int start, int size);

        /// <summary>
        /// Second order derivative at spline corner points
        /// </summary>
        /// <param name="coordinates"></param>
        /// <param name="numNodes"></param>
        /// <param name="coordinatesDerivatives"></param>
        /// <returns></returns>
        static bool SecondOrderDerivative(const std::vector<Point>& coordinates,
                                          int numNodes,
                                          std::vector<Point>& coordinatesDerivatives);

        /// <summary>
        /// Second order derivative at spline corner points (result in a flat array)
        /// </summary>
        /// <param name="coordinates"></param>
        /// <param name="numNodes"></param>
        /// <param name="coordinatesDerivatives"></param>
        /// <returns></returns>
        static bool SecondOrderDerivative(const std::vector<double>& coordinates,
                                          int numNodes,
                                          std::vector<double>& coordinatesDerivatives);

        /// <summary>
        /// Computes the intersection of two splines (sect3r)
        /// </summary>
        /// <param name="first"></param>
        /// <param name="second"></param>
        /// <param name="projection"></param>
        /// <param name="crossProductIntersection"></param>
        /// <param name="intersectionPoint"></param>
        /// <param name="firstSplineRatio"></param>
        /// <param name="secondSplineRatio"></param>
        /// <returns></returns>
        bool GetSplinesIntersection(int first,
                                    int second,
                                    const Projections& projection,
                                    double& crossProductIntersection,
                                    Point& intersectionPoint,
                                    double& firstSplineRatio,
                                    double& secondSplineRatio);

        /// <summary>
        /// Computes the spline length in s coordinates (splinelength)
        /// </summary>
        /// <param name="index"></param>
        /// <param name="beginFactor"></param>
        /// <param name="endFactor"></param>
        /// <param name="numSamples"></param>
        /// <param name="accountForCurvature"></param>
        /// <param name="height"></param>
        /// <param name="assignedDelta"></param>
        /// <returns></returns>
        double GetSplineLength(int index,
                               double beginFactor,
                               double endFactor,
                               int numSamples = 100,
                               bool accountForCurvature = false,
                               double height = 1.0,
                               double assignedDelta = -1);

        std::vector<std::vector<Point>> m_splineCornerPoints;    // The spline corner points
        std::vector<std::vector<Point>> m_splineDerivatives;     // The spline derivatives at the corner points  
        std::vector<int> m_numSplineNodes;                       // Number of spline nodes in each spline
        std::vector<int> m_numAllocatedSplineNodes;              // Number of allocated node in each spline
        std::vector<double> m_splinesLength;                     // The length of each spline
        Projections m_projection;                                // The map projection  
        int m_numSplines = 0;                                    // Current number of splines

        int m_numAllocatedSplines = 0;                           // Total number of allocated splines
        int m_allocationSize = 5;                                // allocation cache size
        
    private:

        /// <summary>
        /// Adds a new corner point in an existing spline
        /// </summary>
        /// <param name="splineIndex">The spline index</param>
        /// <param name="point">The point to add</param>
        /// <returns></returns>
        bool AddPointInExistingSpline(int splineIndex, const Point& point);
        
        /// <summary>
        /// Computes curvature in a point on a spline (comp_curv)
        /// </summary>
        /// <param name="splineIndex"></param>
        /// <param name="adimensionalPointCoordinate"></param>
        /// <param name="curvatureFactor"></param>
        /// <param name="normalVector"></param>
        /// <param name="tangentialVector"></param>
        /// <returns></returns>
        bool ComputeCurvatureOnSplinePoint(int splineIndex,
                                           double adimensionalPointCoordinate,
                                           double& curvatureFactor,
                                           Point& normalVector,
                                           Point& tangentialVector);

        /// <summary>
        /// Delete a spline
        /// </summary>
        /// <param name="splineIndex">The spline index to delete</param>
        /// <returns></returns>
        bool DeleteSpline(int splineIndex);

        /// <summary>
        /// Allocate spline properties arrays 
        /// </summary>
        /// <returns></returns>
        bool AllocateSplinesProperties();

    };
}


