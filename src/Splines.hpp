#pragma once

/// TODO: CHECK FOR FRONT COLLISION

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
        Splines();

        Splines(Projections projection, Polygons& polygon);

        /// add a new spline, return the index
        bool AddSpline(const std::vector<Point>& splines, int start, int size);


        //Run-time settings
        bool SetParameters(const GridGeomApi::CurvilinearParametersNative& curvilinearParametersNative,
            const GridGeomApi::SplinesToCurvilinearParametersNative& splinesToCurvilinearParametersNative);


        // set the reference to polygon
        bool SetPolygon(const Polygons& polygon);


        bool DeleteSpline(int splineIndex);

        /// to be called after all splines have been stored
        bool AllocateSplinesProperties();

        /// add a new spline point in an existing spline
        bool AddPointInExistingSpline(int splineIndex, const Point& point);

        /// spline2curvi
        /// 1. Eliminate spline that are not in polygon
        /// 2. Compute properties (crossings)
        /// 3. Make all grid lines of the central spline
        /// 4. Add artificial splines
        /// 5. Compute properties with artificial spline added
        /// 6. Compute edge velocities
        /// 7. Grow layers
        bool OrthogonalCurvilinearGridFromSplines(CurvilinearGrid& curvilinearGrid);

        bool OrthogonalCurvilinearGridFromSplinesInitialize();

        bool OrthogonalCurvilinearGridFromSplinesIteration(int layer);

        bool OrthogonalCurvilinearGridFromSplinesRefreshMesh(CurvilinearGrid& curvilinearGrid);

        bool ConvertSplineMeshToCurvilinearMesh(const std::vector<std::vector<Point>>& gridPoints, CurvilinearGrid& curvilinearGrid);

        ///get_isub
        bool GetSubIntervalAndGridLayer(int layer, int& gridLayer, int& subLayerIndex);

        /// growlayer
        /// m_edgeVelocities, m_validFrontNodes, m_gridPoints, m_timeStep

        bool GrowLayer(int layerIndex);

        bool ComputeMaximumGridLayerGrowTimeOtherFront();

        ///comp_tmax_self
        bool ComputeMaximumGridLayerGrowTime(const std::vector<Point>& coordinates,
            const std::vector<Point>& velocities,
            std::vector<double>& maximumGridLayerGrowTime);

        /// copy growth velocities to the front, and add points in the front at corners
        /// copy_vel_to_front
        bool CopyVelocitiesToFront(
            const int layerIndex,
            const std::vector<Point>& previousVelocities,
            int& numFrontPoints,
            std::vector<std::vector<int>>& gridPointsIndexses,
            std::vector<Point>& frontGridPoints,
            std::vector<Point>& velocities);

        ///find the frontline of the old (static) grid
        ///findfront
        bool FindFront(
            std::vector<std::vector<int>>& gridPointsIndexses,
            std::vector<Point>& frontGridPoints,
            int& numFrontPoints);

        //comp_vel
        bool ComputeVelocitiesAtGridPoints(int layerIndex, std::vector<Point>& velocityVector);

        /// get_LR
        bool GetNeighbours(const std::vector<Point>& gridPoints, const int index, int& currentLeftIndex, int& currentRightIndex);

        ///comp_edgevel: TODO: can this be splitted in compute heights and computeGrowFactors
        bool ComputeEdgeVelocities(
            std::vector<double>& edgeVelocities, // edgevel
            std::vector<std::vector<double>>& growFactorOnSubintervalAndEdge, //dgrow1
            std::vector<std::vector<int>>& numPerpendicularFacesOnSubintervalAndEdge //nfac1
        );

        ///comp_dgrow: this is another root finding algorithm, could go in the general part
        bool ComputeGrowFactor(
            double totalGridHeight,
            double firstGridLayerHeight,
            double numberOfGridLayers,
            double& result);

        double ComputeTotalExponentialHeight(double aspectRatioGrowFactor, double firstGridLayerHeights, int numberOfGridLayers);

        ///comp_nfac
        ///compute the number of grid layers for a given grow factor, first grid layer height and total grid height
        int ComputeNumberExponentialIntervals(double hhMaxRatio);

        bool ComputeVelocitiesSubIntervals(
            int s,
            int startGridLineIndex,
            int endGridLineIndex,
            int numHeights,
            int numOtherSideHeights,
            double firstHeight,
            const std::vector<int>& gridLineIndex,
            const std::vector<int>& otherGridLineIndex,
            std::vector<std::vector<int>>& numPerpendicularFacesOnSubintervalAndEdge,
            std::vector<double>& edgeVelocities,
            double& hh0MaxRatio);

        /// compute the grid heights at grid edges on the center spline
        ///comp_gridheights
        bool ComputeGridHeights();


        bool FindNearestCrossSplines(int s,
            int j,
            const std::vector<int>& numHeightsLeft,
            const std::vector<double>& edgesCenterPoints,
            const std::vector<std::vector<double>>& crossSplineLeftHeights,
            std::vector<int>& localValidSplineIndexes,
            std::vector<double>& localSplineDerivatives,
            std::vector<double>& crossingSplinesDimensionalCoordinates,
            std::vector<std::vector<double>>& heights);

        // GetValidSplineIndexses
        bool GetValidSplineIndexses(int s, int numValues, const std::vector<int>& v, std::vector<int>& validIndexses, int& numValid);


        /// SECT3R
        /// compute the intersection of two splines
        bool GetSplinesIntersection(int first, int second,
            const Projections& projection,
            double& crossProductIntersection,
            Point& intersectionPoint,
            double& firstSplineRatio,
            double& secondSplineRatio);

        /// get_crosssplines
        /// compute the intersection of two splines, one must have only two nodes
        bool GetSplineIntersections(int index);


        bool MakeAllGridLines(bool isSpacingCurvatureAdapeted);

        /// make_gridline, generate a gridline on a spline with a prescribed maximum mesh width
        /// generate a gridline on a spline with a prescribed maximum mesh width
        bool MakeGridLine(int splineIndex,
            int startingIndex,
            std::vector<Point>& gridLine,
            std::vector<double>& adimensionalCoordinates,
            int& numM);

        ///get_splineprops
        bool ComputeSplineProperties(bool restoreOriginalProperties);


        /// get the grid heights from the cross spline information
        /// get_heights
        bool ComputeHeights();


        ///comp_subheights, compute the height of the subintervals of grid layers on a cross spline, w.r.t. a center spline
        bool ComputeSubHeights(int centerSplineIndex, int crossingSplineLocalIndex);


        ///splinelength
        /// TODO: remove special treatment assign delta and calculate number of points 
        double GetSplineLength(int index,
            double beginFactor,
            double endFactor,
            int numSamples = 100,
            bool accountForCurvature = false,
            double height = 1.0,
            double assignedDelta = -1);


        /// comp_curv
        /// compute curvature in a point on a spline
        bool ComputeCurvatureOnSplinePoint(
            int splineIndex,
            double adimensionalPointCoordinate,
            double& curvatureFactor,
            Point& normalVector,
            Point& tangentialVector);

        int m_numSplines = 0;
        int m_numAllocatedSplines = 0;
        std::vector<int> m_numSplineNodes;
        std::vector<int> m_numAllocatedSplineNodes;
        std::vector<std::vector<Point>> m_splineCornerPoints;
        std::vector<std::vector<Point>> m_splineDerivatives;
        Projections m_projection;
        Polygons m_polygon;

        // OrthogonalCurvilinearGridFromSplines parameters
        double m_aspectRatioFirstLayer = 0.10;
        int m_maxNumM = 20;                                                             // mfacmax
        int m_maxNumN = 40;                                                             // N - refinement factor for regular grid generation.
        int m_numM = 0;                                                                 // mc
        int m_maxNumCenterSplineHeights = 10;                                           // Nsubmax, naz number of different heights a cross spline can have (is determined by how many crossing spline the user can input).
        double m_averageMeshWidth = 500.0;
        double m_dtolcos = 0.95;                                                        // minimum allowed absolute value of crossing - angle cosine
        double m_aspectRatio = 0.1;                                                     // daspect aspect ratio
        double m_aspectRatioGrowFactor = 1.1;                                           // dgrow grow factor of aspect ratio
        double m_gridLayerHeight0 = 10.0;                                               // dheight0 grid layer height
        double m_maxaspect = 1.0;                                                       // maxaspect maximum cell aspect ratio *inoperative*
        int m_maxNUniformPart = 5;                                                      // maximum number of layers in the uniform part
        bool m_growGridOutside = true;                                                  // grow the grid outside the prescribed grid height
        double m_onTopOfEachOtherTolerance = 1e-4;                                      // On - top - of - each - other tolerance *IMPORTANT*
        bool m_checkFrontCollisions = false;                                            // check front collisions
        bool m_isSpacingCurvatureAdapted = true;                                        // is curvature adapted

        // spline types
        enum class SplineTypes
        {
            central,
            crossing,
            arficial,
            lateral
        };

        // Spline properties (first index is the spline number) 
        std::vector<SplineTypes> m_type;
        std::vector<int> m_centralSplineIndex;                                          // for each spline the index to its central
        std::vector<int> m_numCrossingSplines;                                          // ncs num of cross splines
        std::vector<std::vector<int>> m_crossingSplinesIndexses;                        // ics for each cross spline, the indexses of the center splines
        std::vector<double> m_splinesLength;                                            // splinesLength spline path length
        std::vector<double> m_maximumGridHeights;                                       // hmax maximum grid height
        std::vector<std::vector<bool>> m_isLeftOriented;                                // isLeftOriented cross spline is left to right(.true.) or not (.false.) w.r.t.center spline
        std::vector<std::vector<double>>  m_crossSplineCoordinates;                     // t center spline coordinates of cross splines
        std::vector<std::vector<double>> m_cosCrossingAngle;                            // cosPhi cosine of crossing angle
        std::vector<std::vector<std::vector<double>>> m_crossSplineLeftHeights;         // hL left - hand side grid heights at cross spline locations for each grid layer subinterval, hL(1, :) being the height of the first subinterval, etc.
        std::vector<std::vector<std::vector<double>>> m_crossSplineRightHeights;        // hR right - hand side grid heights at cross spline locations for each grid layer subinterval, hR(1, :) being the height of the first subinterval, etc.
        std::vector<std::vector<int>>   m_numCrossSplineLeftHeights;                    // NsubL number of subintervals of grid layers at cross spline locations at the left - hand side of the spline, each having their own exponential grow factor
        std::vector<std::vector<int>>  m_numCrossSplineRightHeights;                    // NsubR number of subintervals of grid layers at cross spline locations at the right - hand side of the spline, each having their own exponential grow factor
        std::vector<int> m_numMSpline;                                                  // mfac number of grid intervals on the spline
        std::vector<std::vector<int>> m_nfacL;                                          // nfacL number of grid layers in each subinterval at the left - hand side of the spline * not used yet*
        std::vector<std::vector<int>> m_nfacR;                                          // nfacR number of grid layers in each subinterval at the right - hand side of the spline * not used yet*
        std::vector<int> m_leftGridLineIndex;                                           // iL index in the whole gridline array of the first grid point on the left - hand side of the spline
        std::vector<int> m_rightGridLineIndex;                                          // iR index in the whole gridline array of the first grid point on the right - hand side of the spline
        std::vector<Point> m_gridLine;                                                  // xg1, yg1 coordinates of the first gridline
        std::vector<double> m_gridLineDimensionalCoordinates;                           // sg1 center spline coordinates of the first gridline
        std::vector<std::vector<double>> m_gridHeights;                                 // heights of all grid elements        


        //original spline chaches
        int m_numOriginalSplines = 0;
        int m_allocationSize = 5;                                                       // allocation cache size

        std::vector<int> m_leftGridLineIndexOriginal;
        std::vector<int> m_rightGridLineIndexOriginal;
        std::vector<int> m_mfacOriginal;
        std::vector<double> m_maximumGridHeightsOriginal;
        std::vector<SplineTypes> m_originalTypes;


        //cache variables during iterations
        std::vector<double> m_edgeVelocities;
        std::vector<int> m_validFrontNodes;
        std::vector<std::vector<Point>> m_gridPoints;
        double m_timeStep = 1.0;
        std::vector<int> m_subLayerGridPoints;
        std::vector<std::vector<int>> m_numPerpendicularFacesOnSubintervalAndEdge;
        std::vector<std::vector<double>> m_growFactorOnSubintervalAndEdge;

        struct FuncDimensionalToAdimensionalDistance
        {
            FuncDimensionalToAdimensionalDistance(Splines& spline, int splineIndex, bool isSpacingCurvatureAdapted, double h) :
                m_spline(spline),
                m_splineIndex(splineIndex),
                m_isSpacingCurvatureAdapted(isSpacingCurvatureAdapted),
                m_h(h)
            {
            };

            void SetDimensionalDistance(const double distance)
            {
                m_DimensionalDistance = distance;
            }

            // this is the function we want to find the root
            double operator()(double const& adimensionalDistancereferencePoint)
            {
                double distanceFromReferencePoint = m_spline.GetSplineLength(m_splineIndex, 0, adimensionalDistancereferencePoint, m_numSamples, m_isSpacingCurvatureAdapted, m_h, 0.1);
                distanceFromReferencePoint = std::abs(distanceFromReferencePoint - m_DimensionalDistance);
                return distanceFromReferencePoint;
            }

            Splines& m_spline;
            int m_splineIndex;
            bool m_isSpacingCurvatureAdapted;
            double m_h;
            int m_numSamples = 10;
            double m_DimensionalDistance;
        };

    };

}
