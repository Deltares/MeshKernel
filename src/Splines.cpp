/// TODO: CHECK FOR FRONT COLLISION

#include <vector>
#include <algorithm>
#include <cassert>
#include "Operations.cpp"
#include "Entities.hpp"
#include "CurvilinearParametersNative.hpp"
#include "SplinesToCurvilinearParametersNative.hpp"
#include "Splines.hpp"
#include "CurvilinearGrid.hpp"

GridGeom::Splines::Splines() : m_projection(Projections::cartesian)
{
};

GridGeom::Splines::Splines(Projections projection, Polygons& polygon) : m_projection(projection), m_polygon(polygon)
{
};

/// add a new spline, return the index
bool GridGeom::Splines::AddSpline(const std::vector<Point>& splines, int start, int size)
{
    ResizeVectorIfNeededWithMinimumSize(m_numSplines + 1, m_splineCornerPoints, m_allocationSize, std::vector<Point>(10, { doubleMissingValue, doubleMissingValue }));

    m_numAllocatedSplines = m_splineCornerPoints.size();
    m_numAllocatedSplineNodes.resize(m_numAllocatedSplines, 10);

    m_numSplineNodes.resize(m_numAllocatedSplines, 0);
    m_numSplineNodes[m_numSplines] = size;

    m_splineDerivatives.resize(m_numAllocatedSplines);
    int index = 0;
    for (int n = start; n < start + size; ++n)
    {
        m_splineCornerPoints[m_numSplines][index] = splines[n];
        index++;
    }

    m_splinesLength.resize(m_numAllocatedSplines);
    m_type.resize(m_numAllocatedSplines);

    // compute basic properties
    SecondOrderDerivative(m_splineCornerPoints[m_numSplines], m_numSplineNodes[m_numSplines], m_splineDerivatives[m_numSplines]);
    m_splinesLength[m_numSplines] = GetSplineLength(m_numSplines, 0, m_numSplineNodes[m_numSplines] - 1);
    m_numSplines++;

    return true;
}

//Run-time settings
bool GridGeom::Splines::SetParameters(const GridGeomApi::CurvilinearParametersNative& curvilinearParametersNative,
    const GridGeomApi::SplinesToCurvilinearParametersNative& splinesToCurvilinearParametersNative)
{
    m_aspectRatio = splinesToCurvilinearParametersNative.AspectRatio;
    m_aspectRatioGrowFactor = splinesToCurvilinearParametersNative.AspectRatioGrowFactor;
    m_averageMeshWidth = splinesToCurvilinearParametersNative.AverageWidth;
    m_onTopOfEachOtherSquaredTolerance = splinesToCurvilinearParametersNative.GridsOnTopOfEachOtherTolerance * 
                                         splinesToCurvilinearParametersNative.GridsOnTopOfEachOtherTolerance;
    m_dtolcos = splinesToCurvilinearParametersNative.MinimumCosineOfCrossingAngles;
    m_checkFrontCollisions = splinesToCurvilinearParametersNative.CheckFrontCollisions;
    m_isSpacingCurvatureAdapted = splinesToCurvilinearParametersNative.CurvatureAdapetedGridSpacing;
    m_removeSkinnyTriangles = splinesToCurvilinearParametersNative.RemoveSkinnyTriangles == 1? true: false;
    
    // curvature adapted grid spacing? to add
    m_maxNumM = curvilinearParametersNative.MRefinement;
    m_maxNumN = curvilinearParametersNative.NRefinement;
    return true;
}

// set the reference to polygon
bool GridGeom::Splines::SetPolygon(const Polygons& polygon)
{
    m_polygon = polygon;
    return true;
}

bool GridGeom::Splines::DeleteSpline(int splineIndex)
{
    m_splineCornerPoints.erase(m_splineCornerPoints.begin() + splineIndex);
    m_numSplineNodes.erase(m_numSplineNodes.begin() + splineIndex);
    m_splineDerivatives.erase(m_splineDerivatives.begin() + splineIndex);
    m_splinesLength.erase(m_splinesLength.begin() + splineIndex);
    m_numSplines--;
    return true;
}

/// to be called after all splines have been stored
bool GridGeom::Splines::AllocateSplinesProperties()
{
    m_type.resize(m_numSplines);

    m_centralSplineIndex.resize(m_numSplines);
    std::fill(m_centralSplineIndex.begin(), m_centralSplineIndex.end(), intMissingValue);

    m_numCrossingSplines.resize(m_numSplines, 0);
    std::fill(m_numCrossingSplines.begin(), m_numCrossingSplines.end(), 0);

    m_maximumGridHeights.resize(m_numSplines);
    std::fill(m_maximumGridHeights.begin(), m_maximumGridHeights.end(), doubleMissingValue);

    // multi-dimensional arrays
    m_crossingSplinesIndexses.resize(m_numSplines);
    m_isLeftOriented.resize(m_numSplines);
    m_crossSplineCoordinates.resize(m_numSplines);
    m_cosCrossingAngle.resize(m_numSplines);
    m_crossSplineLeftHeights.resize(m_numSplines);
    m_crossSplineRightHeights.resize(m_numSplines);
    m_numCrossSplineLeftHeights.resize(m_numSplines);
    m_numCrossSplineRightHeights.resize(m_numSplines);
    m_nfacL.resize(m_numSplines);
    m_nfacR.resize(m_numSplines);
    for (int s = 0; s < m_numSplines; ++s)
    {
        m_crossingSplinesIndexses[s].resize(m_numSplines);
        std::fill(m_crossingSplinesIndexses[s].begin(), m_crossingSplinesIndexses[s].end(), -1);

        m_isLeftOriented[s].resize(m_numSplines, true);
        std::fill(m_isLeftOriented[s].begin(), m_isLeftOriented[s].end(), true);

        m_crossSplineCoordinates[s].resize(m_numSplines);
        std::fill(m_crossSplineCoordinates[s].begin(), m_crossSplineCoordinates[s].end(), doubleMissingValue);

        m_cosCrossingAngle[s].resize(m_numSplines, doubleMissingValue);
        std::fill(m_cosCrossingAngle[s].begin(), m_cosCrossingAngle[s].end(), doubleMissingValue);

        m_crossSplineLeftHeights[s].resize(m_numSplines);
        std::fill(m_crossSplineLeftHeights[s].begin(), m_crossSplineLeftHeights[s].end(), std::vector<double>(m_maxNumCenterSplineHeights, doubleMissingValue));

        m_crossSplineRightHeights[s].resize(m_numSplines);
        std::fill(m_crossSplineRightHeights[s].begin(), m_crossSplineRightHeights[s].end(), std::vector<double>(m_maxNumCenterSplineHeights, doubleMissingValue));

        m_numCrossSplineLeftHeights[s].resize(m_numSplines);
        std::fill(m_numCrossSplineLeftHeights[s].begin(), m_numCrossSplineLeftHeights[s].end(), 0);

        m_numCrossSplineRightHeights[s].resize(m_numSplines, doubleMissingValue);
        std::fill(m_numCrossSplineRightHeights[s].begin(), m_numCrossSplineRightHeights[s].end(), 0);

        m_nfacL[s].resize(m_numSplines);
        std::fill(m_nfacL[s].begin(), m_nfacL[s].end(), 0);

        m_nfacR[s].resize(m_numSplines);
        std::fill(m_nfacR[s].begin(), m_nfacR[s].end(), 0);
    }

    m_numMSpline.resize(m_numSplines);
    std::fill(m_numMSpline.begin(), m_numMSpline.end(), 0);
    m_leftGridLineIndex.resize(m_numSplines);
    std::fill(m_leftGridLineIndex.begin(), m_leftGridLineIndex.end(), intMissingValue);
    m_rightGridLineIndex.resize(m_numSplines);
    std::fill(m_rightGridLineIndex.begin(), m_rightGridLineIndex.end(), intMissingValue);

    return true;
}

/// add a new spline point in an existing spline
bool GridGeom::Splines::AddPointInExistingSpline(int splineIndex, const Point& point)
{
    if (splineIndex > m_numSplines)
    {
        return false;
    }
    ResizeVectorIfNeededWithMinimumSize(m_numSplineNodes[splineIndex] + 1, m_splineCornerPoints[splineIndex], m_allocationSize, { doubleMissingValue, doubleMissingValue });
    m_numAllocatedSplineNodes[splineIndex] = m_splineCornerPoints[splineIndex].size();

    m_splineCornerPoints[splineIndex][m_numSplineNodes[splineIndex]] = point;
    m_numSplineNodes[splineIndex]++;
    return true;
}

/// spline2curvi
/// 1. Eliminate spline that are not in polygon
/// 2. Compute properties (crossings)
/// 3. Make all grid lines of the central spline
/// 4. Add artificial splines
/// 5. Compute properties with artificial spline added
/// 6. Compute edge velocities
/// 7. Grow layers
/// 8. Remove skinny triangles
bool GridGeom::Splines::OrthogonalCurvilinearGridFromSplines(CurvilinearGrid& curvilinearGrid)
{

    bool successful = OrthogonalCurvilinearGridFromSplinesInitialize();
    if (!successful)
    {
        return false;
    }

    // Grow grid, from the second layer
    for (int layer = 1; layer <= m_maxNumN; ++layer)
    {
        successful = OrthogonalCurvilinearGridFromSplinesIteration(layer);
        if (!successful)
        {
            break;
        }
    }

    if (m_removeSkinnyTriangles)
    {
        RemoveSkinnyTriangles();
    }

    successful = OrthogonalCurvilinearGridFromSplinesRefreshMesh(curvilinearGrid);
    return successful;
}


bool GridGeom::Splines::RemoveSkinnyTriangles()
{
    int numMaxIterations = 10;
    int numN = m_gridPoints.size() - 2;
    const double squaredDistanceTolerance = 1e-4;
    const double cosineTolerance = 1e-2;
    const double maxCosine = 0.93969;
    for (int j = numN - 1; j >= 1; --j)
    {
        for (int iter = 0; iter < numMaxIterations; ++iter)
        {
            int numChanged = 0;
            
            int firstLeftIndex;
            int firstRightIndex = 0;
            int i = 0;

            while (firstRightIndex != m_numM - 1 || i != m_numM - 1)
            {
                if (firstRightIndex > i)
                {
                    i = firstRightIndex;
                }
                else
                {
                    i++;
                    if (i >= m_numM - 1)
                    {
                        break;
                    }
                }

                if (!m_gridPoints[j][i].IsValid())
                {
                    continue;
                }

                GetNeighbours(m_gridPoints[j], i, firstLeftIndex, firstRightIndex);

                double squaredRightDistance = ComputeSquaredDistance(m_gridPoints[j][i], m_gridPoints[j][firstRightIndex], m_projection);

                if (squaredRightDistance < squaredDistanceTolerance)
                {
                    continue;
                }

                // Detect triangular cell
                if (!m_gridPoints[j + 1][i].IsValid())
                {
                    continue;
                }

                //GetNeighbours(m_gridPoints[j+1], i, secondLeftIndex, secondRightIndex);
                double squaredLeftDistance = ComputeSquaredDistance(m_gridPoints[j][firstLeftIndex], m_gridPoints[j][i], m_projection);
                if (squaredLeftDistance < squaredDistanceTolerance)
                {
                    firstLeftIndex = i;
                }

                if (m_gridPoints[j + 1][firstRightIndex].IsValid())
                {
                    double squaredCurrentDistance = ComputeSquaredDistance(m_gridPoints[j + 1][i], m_gridPoints[j + 1][firstRightIndex], m_projection);
                    double currentCosPhi = NormalizedInnerProductTwoSegments(
                        m_gridPoints[j + 1][i],
                        m_gridPoints[j][i],
                        m_gridPoints[j + 1][i],
                        m_gridPoints[j][firstRightIndex],
                        m_projection);
                    if (squaredCurrentDistance < squaredDistanceTolerance && currentCosPhi > maxCosine)
                    {

                        //determine persistent node
                        double leftCosPhi = NormalizedInnerProductTwoSegments(
                            m_gridPoints[j - 1][i],
                            m_gridPoints[j][i],
                            m_gridPoints[j][i],
                            m_gridPoints[j + 1][i],
                            m_projection);

                        double rightCosPhi = NormalizedInnerProductTwoSegments(
                            m_gridPoints[j - 1][firstRightIndex],
                            m_gridPoints[j][firstRightIndex],
                            m_gridPoints[j][firstRightIndex],
                            m_gridPoints[j + 1][firstRightIndex],
                            m_projection);


                        int secondLeftIndex;
                        int secondRightIndex;
                        GetNeighbours(m_gridPoints[j], firstRightIndex, secondLeftIndex, secondRightIndex);


                        if ((secondRightIndex == firstRightIndex || leftCosPhi - rightCosPhi < -cosineTolerance) && firstLeftIndex != i)
                        {
                            //move left node
                            for (int k = i; k <= firstRightIndex - 1; ++k)
                            {
                                m_gridPoints[j][k] = m_gridPoints[j][firstRightIndex];
                            }
                            numChanged++;
                        }
                        else if ((firstLeftIndex == i || rightCosPhi - leftCosPhi < -cosineTolerance) && secondRightIndex != firstRightIndex)
                        {
                            //move right node
                            for (int k = firstRightIndex; k <= secondRightIndex - 1; ++k)
                            {
                                m_gridPoints[j][k] = m_gridPoints[j][i];
                            }
                            numChanged++;
                        }
                        else
                        {
                            //move both nodes
                            Point middle = (m_gridPoints[j][i] + m_gridPoints[j][firstRightIndex]) * 0.5;
                            for (int k = i; k <= firstRightIndex - 1; ++k)
                            {
                                m_gridPoints[j][k] = middle;
                            }
                            for (int k = firstRightIndex; k <= secondRightIndex - 1; ++k)
                            {
                                m_gridPoints[j][k] = middle;
                            }
                            numChanged++;
                        }
                    }
                }
            }

            if (numChanged == 0)
            {
                break;
            }
        }
    }

    return true;
}

bool GridGeom::Splines::OrthogonalCurvilinearGridFromSplinesInitialize()
{
    // no splines
    if (m_numSplines < 2)
    {
        return false;
    }

    // Delete the splines that are not fully inside the polygons
    for (int s = 0; s < m_numSplines; s++)
    {
        for (int n = 0; n < m_numSplineNodes[s]; n++)
        {
            bool isInPolygons = m_polygon.IsPointInPolygons(m_splineCornerPoints[s][n]);
            if (!isInPolygons)
            {
                DeleteSpline(s);
                break;
            }
        }
    }

    // compute properties
    bool success = ComputeSplineProperties(false);
    if (!success)
    {
        return false;
    }

    // get the properties of the center splines
    m_numM = 0;
    success = MakeAllGridLines(true);
    if (!success)
    {
        return false;
    }

    // Store original number of splines
    std::vector<Point> newCrossSpline(2);
    m_numOriginalSplines = m_numSplines;
    for (int s = 0; s < m_numOriginalSplines; ++s)
    {
        // mirrow only center splines
        if (m_type[s] != SplineTypes::central)
        {
            continue;
        }

        // construct the cross splines through the edges, along m discretization
        for (int i = m_leftGridLineIndex[s]; i < m_leftGridLineIndex[s] + m_numMSpline[s]; ++i)
        {
            Point normal;
            NormalVectorOutside(m_gridLine[i], m_gridLine[i + 1], normal, m_projection);

            double xMiddle = (m_gridLine[i].x + m_gridLine[i + 1].x) * 0.5;
            double yMiddle = (m_gridLine[i].y + m_gridLine[i + 1].y) * 0.5;
            
            double xs1;
            double xs2;
            double ys1;
            double ys2;

            if (m_projection == Projections::cartesian)
            {
                xs1 = xMiddle + 2.0 * m_maximumGridHeights[s] * -normal.x;
                xs2 = xMiddle + 2.0 * m_maximumGridHeights[s] * normal.x;
                ys1 = yMiddle + 2.0 * m_maximumGridHeights[s] * -normal.y;
                ys2 = yMiddle + 2.0 * m_maximumGridHeights[s] * normal.y;
            }
            if (m_projection == Projections::spherical)
            {
                const double factor = 1.0 / (earth_radius * degrad_hp);
                xs1 = xMiddle + 2.0 * m_maximumGridHeights[s] * -normal.x * factor;
                xs2 = xMiddle + 2.0 * m_maximumGridHeights[s] * normal.x * factor;
                ys1 = yMiddle + 2.0 * m_maximumGridHeights[s] * -normal.y * factor;
                ys2 = yMiddle + 2.0 * m_maximumGridHeights[s] * normal.y * factor;
            }

            newCrossSpline[0] = { xs1, ys1 };
            newCrossSpline[1] = { xs2, ys2 };
            AddSpline(newCrossSpline, 0, newCrossSpline.size());
            // flag the cross spline as artificially added
            m_type[m_numSplines - 1] = SplineTypes::arficial;
        }
    }

    // Backup original spline properties
    m_leftGridLineIndexOriginal.resize(m_numOriginalSplines);
    m_rightGridLineIndexOriginal.resize(m_numOriginalSplines);
    m_mfacOriginal.resize(m_numOriginalSplines);
    m_maximumGridHeightsOriginal.resize(m_numOriginalSplines);
    m_originalTypes.resize(m_numOriginalSplines);
    for (int s = 0; s < m_numOriginalSplines; ++s)
    {
        m_leftGridLineIndexOriginal[s] = m_leftGridLineIndex[s];
        m_rightGridLineIndexOriginal[s] = m_rightGridLineIndex[s];
        m_mfacOriginal[s] = m_numMSpline[s];
        m_maximumGridHeightsOriginal[s] = m_maximumGridHeights[s];
        m_originalTypes[s] = m_type[s];
    }

    // compute spline properties with artificial splines added
    ComputeSplineProperties(true);

    // artificial cross spline: remove the last part of the sub-intervals (since it makes no sense, 
    // as the artificial cross spline has an arbitrary, but sufficiently large, length)
    for (int s = 0; s < m_numOriginalSplines; ++s)
    {
        // Remove the last part of the sub-intervals
        if (m_type[s] != SplineTypes::central)
        {
            continue;
        }

        // For number of intersecting splines
        for (int i = 0; i < m_numCrossingSplines[s]; ++i)
        {
            int crossingSplineIndex = m_crossingSplinesIndexses[s][i];
            if (m_type[crossingSplineIndex] == SplineTypes::arficial)
            {
                m_numCrossSplineLeftHeights[s][i] = m_numCrossSplineLeftHeights[s][i] - 1;
                m_numCrossSplineRightHeights[s][i] = m_numCrossSplineRightHeights[s][i] - 1;
            }
        }
    }

    // Compute edge velocities
    m_edgeVelocities.resize(m_numM - 1, doubleMissingValue);
    m_growFactorOnSubintervalAndEdge.resize(m_maxNumCenterSplineHeights, std::vector<double>(m_numM - 1, 1.0));
    m_numPerpendicularFacesOnSubintervalAndEdge.resize(m_maxNumCenterSplineHeights, std::vector<int>(m_numM - 1, 0));
    success = ComputeEdgeVelocities(m_edgeVelocities, m_growFactorOnSubintervalAndEdge, m_numPerpendicularFacesOnSubintervalAndEdge);
    if (!success)
    {
        return false;
    }

    // Increase curvilinear grid
    const int numGridLayers = m_maxNumN + 1;
    // The layer by coordinate to grow
    m_gridPoints.resize(numGridLayers + 1, std::vector<Point>(m_numM + 1, { doubleMissingValue, doubleMissingValue }));
    m_validFrontNodes.resize(m_numM, 1);

    // Copy the first n in m_gridPoints
    for (int n = 0; n < m_numM; ++n)
    {
        m_gridPoints[0][n] = m_gridLine[n];
        if (m_gridLine[n].x == doubleMissingValue)
        {
            m_validFrontNodes[n] = 0;
        }
        int sumLeft = 0;
        int sumRight = 0;
        int leftColumn = std::max(n - 1, 0);
        int rightColumn = std::min(n, m_numM - 2);
        for (int j = 0; j < m_numPerpendicularFacesOnSubintervalAndEdge.size(); ++j)
        {
            sumLeft += m_numPerpendicularFacesOnSubintervalAndEdge[j][leftColumn];
            sumRight += m_numPerpendicularFacesOnSubintervalAndEdge[j][rightColumn];
        }
        if (sumLeft == 0 && sumRight == 0)
        {
            m_validFrontNodes[n] = 0;
        }
    }

    //compute maximum mesh width and get dtolLR in the proper dimension
    double squaredMaximumGridWidth = 0.0;
    for (int i = 0; i < m_gridPoints[0].size() - 1; i++)
    {
        if (!m_gridPoints[0][i].IsValid() || !m_gridPoints[0][i + 1].IsValid())
        {
            continue;
        }
        squaredMaximumGridWidth = std::max(squaredMaximumGridWidth, ComputeSquaredDistance(m_gridPoints[0][i], m_gridPoints[0][i + 1], m_projection));
    }
    m_onTopOfEachOtherSquaredTolerance = m_onTopOfEachOtherSquaredTolerance * squaredMaximumGridWidth;


    m_subLayerGridPoints.resize(m_numPerpendicularFacesOnSubintervalAndEdge.size());

    return true;
}

bool GridGeom::Splines::OrthogonalCurvilinearGridFromSplinesIteration(int layer)
{
    bool successful = GrowLayer(layer);
    if (!successful)
    {
        return false;
    }

    for (int j = 0; j < m_subLayerGridPoints.size(); ++j)
    {
        m_subLayerGridPoints[j] = m_numPerpendicularFacesOnSubintervalAndEdge[j][0];
    }

    int gridLayer;
    int subLayerRightIndex;

    successful = GetSubIntervalAndGridLayer(layer, gridLayer, subLayerRightIndex);
    if (!successful)
    {
        return false;
    }

    for (int i = 0; i < m_numM; i++)
    {
        int subLayerLeftIndex = subLayerRightIndex;
        int minRight = std::min(i, int(m_numPerpendicularFacesOnSubintervalAndEdge[0].size() - 1));
        for (int j = 0; j < m_subLayerGridPoints.size(); ++j)
        {
            m_subLayerGridPoints[j] = m_numPerpendicularFacesOnSubintervalAndEdge[j][minRight];
        }

        successful = GetSubIntervalAndGridLayer(layer, gridLayer, subLayerRightIndex);
        if (!successful)
        {
            return false;
        }

        if (subLayerRightIndex >= 0 && i < m_numM - 1 && gridLayer >= 0)
        {
            m_edgeVelocities[i] = m_growFactorOnSubintervalAndEdge[subLayerRightIndex][i] * m_edgeVelocities[i];
        }

        if (subLayerLeftIndex < 0 && subLayerRightIndex < 0)
        {
            m_validFrontNodes[i] = -1;
        }
    }

    if (m_timeStep < 1e-8)
    {
        return false;
    }

    return true;
}

bool GridGeom::Splines::OrthogonalCurvilinearGridFromSplinesRefreshMesh(CurvilinearGrid& curvilinearGrid)
{
    bool successful = ConvertSplineMeshToCurvilinearMesh(m_gridPoints, curvilinearGrid);
    return successful;
}

bool GridGeom::Splines::ConvertSplineMeshToCurvilinearMesh(const std::vector<std::vector<Point>>& gridPoints, CurvilinearGrid& curvilinearGrid)
{

    std::vector<std::vector<int>> mIndexsesThisSide(1, std::vector<int>(2));
    std::vector<std::vector<int>> mIndexsesOtherSide(1, std::vector<int>(2));
    std::vector<std::vector<int>> nIndexsesThisSide(1, std::vector<int>(2));
    std::vector<std::vector<int>> nIndexsesOtherSide(1, std::vector<int>(2));
    std::vector<std::vector<Point>> gridPointsNDirection(gridPoints[0].size(), std::vector<Point>(gridPoints.size()));
    std::vector<std::vector<Point>> curvilinearMeshPoints;
    double squaredDistanceTolerance = 1e-12;

    // get the grid sizes in j-direction
    for (int i = 0; i < gridPoints[0].size(); i++)
    {
        for (int j = 0; j < gridPoints.size(); j++)
        {
            gridPointsNDirection[i][j] = gridPoints[j][i];

        }
    }

    int startIndex = 0;
    int startGridLine = 0;
    while (startIndex < gridPoints[0].size())
    {
        int pos = FindIndexes(gridPoints[0], startIndex, m_numM, doubleMissingValue, mIndexsesThisSide);

        mIndexsesOtherSide[0][0] = mIndexsesThisSide[0][1] + 2;
        mIndexsesOtherSide[0][1] = mIndexsesOtherSide[0][0] + (mIndexsesThisSide[0][1] - mIndexsesThisSide[0][0]);
        bool isConnected = true;

        int minN = m_maxNumN;
        int maxN = 0;
        int minNOther = m_maxNumN;
        int maxNOther = 0;
        //check if this part is connected to another part
        for (int i = mIndexsesThisSide[0][0]; i < mIndexsesThisSide[0][1] + 1; ++i)
        {
            pos = FindIndexes(gridPointsNDirection[i], 0, gridPointsNDirection[i].size(), doubleMissingValue, nIndexsesThisSide);
            minN = std::min(minN, nIndexsesThisSide[0][0]);
            maxN = std::max(maxN, nIndexsesThisSide[0][1]);

            int mOther = mIndexsesThisSide[0][1] + 2 + (mIndexsesThisSide[0][1] - i);

            if (mOther > m_numM - 1)
            {
                // no more grid available
                isConnected = false;
            }
            else
            {
                double squaredDistance = ComputeSquaredDistance(gridPoints[0][i], gridPoints[0][mOther], m_projection);
                if (squaredDistance > squaredDistanceTolerance)
                {
                    isConnected = false;
                }
                else
                {
                    pos = FindIndexes(gridPointsNDirection[mOther], 0, gridPointsNDirection[mOther].size(), doubleMissingValue, nIndexsesOtherSide);
                    minNOther = std::min(minNOther, nIndexsesOtherSide[0][0]);
                    maxNOther = std::max(maxNOther, nIndexsesOtherSide[0][1]);
                }
            }
        }

        const int endGridlineIndex = startGridLine + mIndexsesThisSide[0][1] - mIndexsesThisSide[0][0];
        if (isConnected)
        {
            startIndex = mIndexsesOtherSide[0][1] + 2;
        }
        else
        {
            maxNOther = 1;
            startIndex = mIndexsesThisSide[0][1] + 2;
        }

        // increment points
        int oldMSize = curvilinearMeshPoints.size();
        curvilinearMeshPoints.resize(endGridlineIndex + 1);
        int NSize = std::max(int(curvilinearMeshPoints[0].size()), maxN + maxNOther + 1);
        for (int i = 0; i < curvilinearMeshPoints.size(); i++)
        {
            curvilinearMeshPoints[i].resize(NSize);
        }

        // fill first part
        int columnIncrement = 0;
        for (int i = startGridLine; i < endGridlineIndex + 1; ++i)
        {
            for (int j = 0; j < maxN + 1; ++j)
            {
                curvilinearMeshPoints[i][j + maxNOther] = gridPoints[j][mIndexsesThisSide[0][0] + columnIncrement];
            }
            columnIncrement++;
        }

        columnIncrement = 0;
        for (int i = startGridLine; i < endGridlineIndex + 1; ++i)
        {
            for (int j = 0; j < maxNOther + 1; ++j)
            {
                curvilinearMeshPoints[i][maxNOther - j] = gridPoints[j][mIndexsesOtherSide[0][1] - columnIncrement];
            }
            columnIncrement++;
        }

        startGridLine = endGridlineIndex + 2;
    }

    curvilinearGrid.IncreaseGrid(curvilinearMeshPoints.size(), curvilinearMeshPoints[0].size());
    curvilinearGrid.Set(curvilinearMeshPoints);

    return true;
}


///get_isub
bool GridGeom::Splines::GetSubIntervalAndGridLayer(int layer, int& gridLayer, int& subLayerIndex)
{

    gridLayer = layer - 1;
    int sum = std::accumulate(m_subLayerGridPoints.begin(), m_subLayerGridPoints.end(), 0);
    if (layer >= sum)
    {
        subLayerIndex = -1;
        return true;
    }

    subLayerIndex = 0;
    sum = m_subLayerGridPoints[0] + 1;
    while (sum <= layer && subLayerIndex < m_maxNumCenterSplineHeights)
    {
        subLayerIndex = subLayerIndex + 1;
        sum += m_subLayerGridPoints[subLayerIndex];
    }
    gridLayer = layer - sum + m_subLayerGridPoints[subLayerIndex];

    return true;
}

/// growlayer
bool GridGeom::Splines::GrowLayer(int layerIndex)
{
    assert(layerIndex - 1 >= 0);
    std::vector<Point> velocityVectorAtGridPoints(m_numM);
    bool success = ComputeVelocitiesAtGridPoints(layerIndex - 1, velocityVectorAtGridPoints);
    if (!success)
    {
        return false;
    }

    std::vector<Point> activeLayerPoints(m_gridPoints[layerIndex - 1]);
    for (int m = 0; m < velocityVectorAtGridPoints.size(); ++m)
    {
        if (!velocityVectorAtGridPoints[m].IsValid())
        {
            m_gridPoints[layerIndex - 1][m] = { doubleMissingValue, doubleMissingValue };
            activeLayerPoints[m] = { doubleMissingValue, doubleMissingValue };
        }
    }

    int numGridPoints = m_gridPoints.size() * m_gridPoints[0].size();
    std::vector<std::vector<int>> gridPointsIndexses(numGridPoints, std::vector<int>(2, -1));
    std::vector<Point> frontGridPoints(numGridPoints);
    int numFrontPoints;
    success = FindFront(gridPointsIndexses, frontGridPoints, numFrontPoints);
    if (!success)
    {
        return false;
    }

    std::vector<Point> frontVelocities(numGridPoints);
    success = CopyVelocitiesToFront(layerIndex - 1, velocityVectorAtGridPoints, numFrontPoints,
        gridPointsIndexses, frontGridPoints, frontVelocities);

    if (!success)
    {
        return false;
    }

    double totalTimeStep = 0.0;
    std::vector<Point> gridLine(m_gridPoints[layerIndex - 1]);
    double localTimeStep = 0.0;
    double otherTimeStep = std::numeric_limits<double>::max();
    std::vector<int> newValidFrontNodes(numGridPoints);

    while (totalTimeStep < m_timeStep)
    {
        // Copy old front velocities
        newValidFrontNodes = m_validFrontNodes;

        for (int i = 0; i < m_validFrontNodes.size(); ++i)
        {
            if(m_validFrontNodes[i]<=0)
            {
                activeLayerPoints[i] = { doubleMissingValue, doubleMissingValue };
            }
        }

        std::vector<double> maximumGridLayerGrowTime(newValidFrontNodes.size(), std::numeric_limits<double>::max());
        success = ComputeMaximumGridLayerGrowTime(activeLayerPoints, velocityVectorAtGridPoints, maximumGridLayerGrowTime);
        if (!success)
        {
            return false;
        }
        localTimeStep = std::min(m_timeStep - totalTimeStep, *std::min_element(maximumGridLayerGrowTime.begin(), maximumGridLayerGrowTime.end()));

        if (m_checkFrontCollisions)
        {
            //TODO: implement front collisions
            otherTimeStep = 0;
        }
        localTimeStep = std::min(localTimeStep, otherTimeStep);

        // remove isolated points at the start end end of the masl
        if (newValidFrontNodes[0] == 1 && newValidFrontNodes[1] == 0)
        {
            newValidFrontNodes[0] = 0;
        }

        if (newValidFrontNodes[m_numM - 1] == 1 && newValidFrontNodes[m_numM - 2] == 0)
        {
            newValidFrontNodes[m_numM - 1] = 0;
        }

        for (int i = 0; i < newValidFrontNodes.size() - 2; ++i)
        {
            if (newValidFrontNodes[i + 1] == 1 && newValidFrontNodes[i] == 0 && newValidFrontNodes[i + 2] == 0)
            {
                newValidFrontNodes[i + 1] = 0;
            }
        }

        m_validFrontNodes = newValidFrontNodes;

        for (int i = 0; i < velocityVectorAtGridPoints.size(); ++i)
        {
            if (m_validFrontNodes[i] == 1 && velocityVectorAtGridPoints[i].IsValid())
            {
                if (velocityVectorAtGridPoints[i].x == 0.0 && velocityVectorAtGridPoints[i].y == 0.0)
                {
                    continue;
                }
                activeLayerPoints[i].x = activeLayerPoints[i].x + localTimeStep * velocityVectorAtGridPoints[i].x;
                activeLayerPoints[i].y = activeLayerPoints[i].y + localTimeStep * velocityVectorAtGridPoints[i].y;
            }
            else
            {
                activeLayerPoints[i].x = doubleMissingValue;
                activeLayerPoints[i].y = doubleMissingValue;
            }
        }

        // update the grid points
        m_gridPoints[layerIndex] = activeLayerPoints;

        // update the time step
        totalTimeStep += localTimeStep;

        if (totalTimeStep < m_timeStep)
        {
            success = ComputeVelocitiesAtGridPoints(layerIndex, velocityVectorAtGridPoints);
            if (!success)
            {
                return false;
            }

            for (int i = 0; i < m_numM; ++i)
            {
                // Disable points that have no valid normal vector
                // Remove stationary points
                if (!frontVelocities[i].IsValid() || m_validFrontNodes[i] == 0)
                {
                    activeLayerPoints[i] = { doubleMissingValue, doubleMissingValue };
                }
            }

            success = FindFront(gridPointsIndexses, frontGridPoints, numFrontPoints);
            if (!success)
            {
                return false;
            }

            success = CopyVelocitiesToFront(layerIndex - 1, velocityVectorAtGridPoints, numFrontPoints,
                gridPointsIndexses, frontGridPoints, frontVelocities);

            if (!success)
            {
                return false;
            }
        }
    }

    if (layerIndex >= 2)
    {
        for (int i = 1; i < m_numM - 1; ++i)
        {

            if (!activeLayerPoints[i].IsValid())
            {
                continue;
            }
            double cosphi = NormalizedInnerProductTwoSegments(m_gridPoints[layerIndex - 2][i],
                m_gridPoints[layerIndex - 1][i],
                m_gridPoints[layerIndex - 1][i],
                activeLayerPoints[i],
                m_projection);
            if (cosphi < -0.5)
            {
                int currentLeftIndex;
                int currentRightIndex;
                GetNeighbours(frontGridPoints, i, currentLeftIndex, currentRightIndex);
                for (int j = currentLeftIndex + 1; j < currentRightIndex; ++j)
                {
                    newValidFrontNodes[j] = 0;
                    m_gridPoints[layerIndex - 1][j] = { doubleMissingValue, doubleMissingValue };
                }
            }
        }
    }

    m_validFrontNodes = newValidFrontNodes;


    return true;
}

bool GridGeom::Splines::ComputeMaximumGridLayerGrowTimeOtherFront()
{

    return true;
}


///comp_tmax_self
bool GridGeom::Splines::ComputeMaximumGridLayerGrowTime
(const std::vector<Point>& coordinates,
    const std::vector<Point>& velocities,
    std::vector<double>& maximumGridLayerGrowTime)
{
    std::vector<double> edgeWidth(coordinates.size() - 1);
    std::vector<double> edgeIncrement(coordinates.size() - 1);
    double minEdgeWidth = 1e-8;
    double dt = 1.0;
    for (int i = 0; i < coordinates.size() - 1; ++i)
    {
        if (!coordinates[i].IsValid() || !coordinates[i + 1].IsValid())
        {
            continue;
        }

        edgeWidth[i] = Distance(coordinates[i], coordinates[i + 1], m_projection);

        if (edgeWidth[i] < minEdgeWidth)
        {
            continue;
        }

        Point firstPointIncremented(coordinates[i] + velocities[i] * dt);
        Point secondPointIncremented(coordinates[i + 1] + velocities[i + 1] * dt);;
        edgeIncrement[i] = InnerProductTwoSegments
        (
            coordinates[i],
            coordinates[i + 1],
            firstPointIncremented,
            secondPointIncremented, m_projection) / edgeWidth[i] - edgeWidth[i];

        edgeIncrement[i] = edgeIncrement[i] / dt;
    }

    for (int i = 0; i < coordinates.size() - 1; ++i)
    {
        if (edgeIncrement[i] < 0.0)
        {
            maximumGridLayerGrowTime[i] = -edgeWidth[i] / edgeIncrement[i];
        }
    }

    return true;
}


/// copy growth velocities to the front, and add points in the front at corners
/// copy_vel_to_front
bool GridGeom::Splines::CopyVelocitiesToFront(
    const int layerIndex,
    const std::vector<Point>& previousVelocities,
    int& numFrontPoints,
    std::vector<std::vector<int>>& gridPointsIndexses,
    std::vector<Point>& frontGridPoints,
    std::vector<Point>& velocities)
{
    int numCornerNodes = 0;
    int p = -1;
    while (p < numFrontPoints)
    {
        p = p + 1;
        if (gridPointsIndexses[p][1] == layerIndex && m_validFrontNodes[gridPointsIndexses[p][0]] == 1)
        {
            velocities[p] = previousVelocities[gridPointsIndexses[p][0]];
            if (!velocities[p].IsValid())
            {
                velocities[p] = { 0.0,0.0 };
            }

            // Check for cornernodes
            int previous = std::max(p - 1, 0);
            std::vector<int> previousIndexses = gridPointsIndexses[previous];
            int next = std::min(p + 1, numFrontPoints);
            std::vector<int> nextIndexses = gridPointsIndexses[next];

            // Check corner nodes
            bool ll = previousIndexses[0] == gridPointsIndexses[p][0] - 1 &&
                previousIndexses[1] == gridPointsIndexses[p][1] &&
                m_validFrontNodes[previousIndexses[0]] == -1;

            bool lr = nextIndexses[0] == gridPointsIndexses[p][0] + 1 &&
                nextIndexses[1] == gridPointsIndexses[p][1] &&
                m_validFrontNodes[nextIndexses[0]] == -1;

            ll = ll || previousIndexses[0] == gridPointsIndexses[p][0] && previousIndexses[1] < gridPointsIndexses[p][1];
            lr = lr || nextIndexses[0] == gridPointsIndexses[p][0] && nextIndexses[1] < gridPointsIndexses[p][1];
            if (ll || lr)
            {
                numCornerNodes++;
                if (numFrontPoints + 1 > frontGridPoints.size())
                {
                    continue;
                }
                for (int i = numFrontPoints; i >= p; --i)
                {
                    frontGridPoints[i + 1] = frontGridPoints[i];
                    velocities[i + 1] = velocities[i];
                    gridPointsIndexses[i + 1] = gridPointsIndexses[i];
                }
                numFrontPoints++;

                if (ll)
                {
                    velocities[p] = { 0.0,0.0 };
                }
                else
                {
                    velocities[p + 1] = { 0.0,0.0 };
                }
                p = p + 1;
            }
        }
    }
    return true;
}


///find the frontline of the old (static) grid
///findfront
bool GridGeom::Splines::FindFront(
    std::vector<std::vector<int>>& gridPointsIndexses,
    std::vector<Point>& frontGridPoints,
    int& numFrontPoints)
{

    std::vector<int> frontPosition(m_gridPoints[0].size() - 2, m_gridPoints.size());
    for (int m = 0; m < frontPosition.size(); ++m)
    {
        for (int n = 0; n < m_gridPoints.size(); ++n)
        {
            if (!m_gridPoints[n][m].IsValid() || !m_gridPoints[n][m + 1].IsValid())
            {
                frontPosition[m] = n - 1;
                break;
            }
        }
    }

    numFrontPoints = 0;
    // check for circular connectivity
    int currentLeftIndex;
    int currentRightIndex;
    int previousFrontPosition = 0;
    GetNeighbours(m_gridPoints[0], 0, currentLeftIndex, currentRightIndex);
    if (currentLeftIndex == 0)
    {
        frontGridPoints[0] = m_gridPoints[0][0];
        // store front index
        gridPointsIndexses[numFrontPoints][0] = 0;
        gridPointsIndexses[numFrontPoints][1] = 0;
        numFrontPoints++;
    }
    else
    {
        previousFrontPosition = frontPosition[currentLeftIndex];
        frontGridPoints[numFrontPoints] = m_gridPoints[0][frontPosition[0]];
        gridPointsIndexses[numFrontPoints][0] = frontPosition[0];
        gridPointsIndexses[numFrontPoints][1] = 0;
        numFrontPoints++;
    }

    for (int m = 0; m < m_gridPoints[0].size() - 2; ++m)
    {
        GetNeighbours(m_gridPoints[0], m, currentLeftIndex, currentRightIndex);
        int currentFrontPosition = frontPosition[m];
        if (currentFrontPosition >= 0)
        {
            if (previousFrontPosition == -1)
            {
                frontGridPoints[numFrontPoints] = m_gridPoints[0][m];
                gridPointsIndexses[numFrontPoints][0] = m;
                gridPointsIndexses[numFrontPoints][1] = 0;
                numFrontPoints++;
            }
            for (int i = previousFrontPosition + 1; i <= currentFrontPosition; ++i)
            {
                frontGridPoints[numFrontPoints] = m_gridPoints[i][m];
                gridPointsIndexses[numFrontPoints][0] = m;
                gridPointsIndexses[numFrontPoints][1] = i;
                numFrontPoints++;
            }
            for (int i = previousFrontPosition; i > currentFrontPosition; --i)
            {
                frontGridPoints[numFrontPoints] = m_gridPoints[i][m];
                gridPointsIndexses[numFrontPoints][0] = m;
                gridPointsIndexses[numFrontPoints][1] = i;
                numFrontPoints++;
            }

            frontGridPoints[numFrontPoints] = m_gridPoints[currentFrontPosition][m + 1];
            gridPointsIndexses[numFrontPoints][0] = m + 1;
            gridPointsIndexses[numFrontPoints][1] = currentFrontPosition;
            numFrontPoints++;
        }
        else if (previousFrontPosition >= 0)
        {
            for (int i = previousFrontPosition - 1; i >= 0; --i)
            {
                frontGridPoints[numFrontPoints] = m_gridPoints[i][m];
                gridPointsIndexses[numFrontPoints][0] = m;
                gridPointsIndexses[numFrontPoints][1] = i;
                numFrontPoints++;
            }

            frontGridPoints[numFrontPoints] = { doubleMissingValue,doubleMissingValue };
            gridPointsIndexses[numFrontPoints][0] = m;
            gridPointsIndexses[numFrontPoints][1] = -1;
            numFrontPoints++;
        }

        previousFrontPosition = currentFrontPosition;
    }

    // add last j-edge, check for circular connectivity
    int lastPoint = m_gridPoints[0].size() - 2;
    GetNeighbours(m_gridPoints[0], lastPoint, currentLeftIndex, currentRightIndex);
    if (currentRightIndex == m_gridPoints[0].size() - 2)
    {
        for (int i = previousFrontPosition; i >= 0; --i)
        {
            frontGridPoints[numFrontPoints] = m_gridPoints[i][lastPoint];
            gridPointsIndexses[numFrontPoints][0] = lastPoint;
            gridPointsIndexses[numFrontPoints][1] = i;
            numFrontPoints++;
        }

    }

    return true;
}


//comp_vel
bool GridGeom::Splines::ComputeVelocitiesAtGridPoints
(
    int layerIndex,
    std::vector<Point>& velocityVector
)
{
    std::fill(velocityVector.begin(), velocityVector.end(), Point{ doubleMissingValue,doubleMissingValue });
    Point normalVectorLeft;
    Point normalVectorRight;
    const double cosTolerance = 1e-8;
    double eps = 1e-10;
    for (int m = 0; m < velocityVector.size(); ++m)
    {
        if (!m_gridPoints[layerIndex][m].IsValid())
        {
            continue;
        }

        int currentLeftIndex;
        int currentRightIndex;
        GetNeighbours(m_gridPoints[layerIndex], m, currentLeftIndex, currentRightIndex);

        double squaredLeftRightDistance = ComputeSquaredDistance(m_gridPoints[layerIndex][currentLeftIndex], m_gridPoints[layerIndex][currentRightIndex], m_projection);
        double squaredLeftDistance = ComputeSquaredDistance(m_gridPoints[layerIndex][currentLeftIndex], m_gridPoints[layerIndex][m], m_projection);
        double squaredRightDistance = ComputeSquaredDistance(m_gridPoints[layerIndex][currentRightIndex], m_gridPoints[layerIndex][m], m_projection);

        if (squaredLeftRightDistance <= m_onTopOfEachOtherSquaredTolerance)
        {
            continue;
        }

        if (squaredLeftDistance <= m_onTopOfEachOtherSquaredTolerance || squaredRightDistance <= m_onTopOfEachOtherSquaredTolerance)
        {
            NormalVectorOutside(m_gridPoints[layerIndex][currentRightIndex], m_gridPoints[layerIndex][currentLeftIndex], normalVectorLeft, m_projection);
            if (m_projection == Projections::spherical)
            {
                normalVectorLeft.x = normalVectorLeft.x * std::cos(degrad_hp * 0.5 * (m_gridPoints[layerIndex][currentLeftIndex].y + m_gridPoints[layerIndex][currentRightIndex].y));
            }
            normalVectorRight = normalVectorLeft;
        }
        else
        {
            NormalVectorOutside(m_gridPoints[layerIndex][m], m_gridPoints[layerIndex][currentLeftIndex], normalVectorLeft, m_projection);
            NormalVectorOutside(m_gridPoints[layerIndex][currentRightIndex], m_gridPoints[layerIndex][m], normalVectorRight, m_projection);

            if (m_projection == Projections::spherical)
            {
                normalVectorLeft.x = normalVectorLeft.x * std::cos(degrad_hp * 0.5 * (m_gridPoints[layerIndex][currentLeftIndex].y + m_gridPoints[layerIndex][m].y));
                normalVectorRight.x = normalVectorRight.x * std::cos(degrad_hp * 0.5 * (m_gridPoints[layerIndex][currentRightIndex].y + m_gridPoints[layerIndex][m].y));
            }
        }

        if (currentLeftIndex == velocityVector.size() - 1)
        {
            continue;
        }

        double cosphi = DotProduct(normalVectorLeft.x, normalVectorRight.x, normalVectorLeft.y, normalVectorRight.y);
        Point leftVelocity = normalVectorLeft * m_edgeVelocities[currentLeftIndex];
        Point rightVelocity = normalVectorRight * m_edgeVelocities[currentRightIndex - 1];
        double rightLeftVelocityRatio = m_edgeVelocities[currentRightIndex - 1] / m_edgeVelocities[currentLeftIndex];

        if (cosphi < -1.0 + cosTolerance)
        {
            continue;
        }

        if (rightLeftVelocityRatio - cosphi > eps && 1.0 / rightLeftVelocityRatio - cosphi > eps || cosphi <= cosTolerance)
        {
            velocityVector[m] = (leftVelocity * (1.0 - rightLeftVelocityRatio * cosphi) +
                rightVelocity * (1.0 - 1.0 / rightLeftVelocityRatio * cosphi)) / (1.0 - cosphi * cosphi);
        }
        else if (cosphi - rightLeftVelocityRatio > eps)
        {
            velocityVector[m] = leftVelocity * rightLeftVelocityRatio / cosphi;
        }
        else
        {
            velocityVector[m] = rightVelocity * 1.0 / (rightLeftVelocityRatio * cosphi);
        }

        if (m_projection == Projections::spherical)
        {
            velocityVector[m].x = velocityVector[m].x * one_over_earth_radius * raddeg_hp / std::cos(degrad_hp * m_gridPoints[layerIndex][m].y);
            velocityVector[m].y = velocityVector[m].y * one_over_earth_radius * raddeg_hp;
        }
    }
    return true;
}

/// get_LR
bool GridGeom::Splines::GetNeighbours(
    const std::vector<Point>& gridPoints,
    int index,
    int& currentLeftIndex,
    int& currentRightIndex)
{
    bool circularConnection = false;
    currentLeftIndex = index;
    currentRightIndex = index;
    int start = 0;
    int end = gridPoints.size() - 1;

    // left
    while (ComputeSquaredDistance(gridPoints[currentLeftIndex], gridPoints[index], m_projection) < m_onTopOfEachOtherSquaredTolerance)
    {
        if (!circularConnection)
        {
            if (currentLeftIndex - 1 < 0)
            {
                break;
            }
        }
        else if (currentLeftIndex - 1 < 0)
        {
            currentLeftIndex = end + 1;
            circularConnection = false;
        }
        if (currentLeftIndex - 1 < 0 || !gridPoints[currentLeftIndex - 1].IsValid())
        {
            break;
        }
        currentLeftIndex--;
    }

    // right
    while (ComputeSquaredDistance(gridPoints[currentRightIndex], gridPoints[index], m_projection) < m_onTopOfEachOtherSquaredTolerance)
    {
        if (!circularConnection)
        {
            if (currentRightIndex + 1 > gridPoints.size())
            {
                break;
            }
        }
        else if (currentRightIndex + 1 > gridPoints.size())
        {
            currentRightIndex = start - 1;
            circularConnection = false;
        }
        if (currentRightIndex + 1 >= gridPoints.size() || !gridPoints[currentRightIndex + 1].IsValid())
        {
            break;
        }
        currentRightIndex++;
    }
    return true;
}

///comp_edgevel: TODO: can this be splitted in compute heights and computeGrowFactors
bool GridGeom::Splines::ComputeEdgeVelocities(
    std::vector<double>& edgeVelocities, // edgevel
    std::vector<std::vector<double>>& growFactorOnSubintervalAndEdge, //dgrow1
    std::vector<std::vector<int>>& numPerpendicularFacesOnSubintervalAndEdge //nfac1
)
{

    bool success = ComputeGridHeights();
    if (!success)
    {
        return false;
    }

    for (int i = 0; i < numPerpendicularFacesOnSubintervalAndEdge[0].size(); ++i)
    {
        numPerpendicularFacesOnSubintervalAndEdge[0][i] = 1;
    }

    for (int s = 0; s < m_numSplines; s++)
    {
        double maxHeight = std::numeric_limits<double>::min();
        for (int i = 0; i < m_gridHeights[0].size(); ++i)
        {
            if (m_gridHeights[0][i] != doubleMissingValue && m_gridHeights[0][i] > maxHeight)
            {
                maxHeight = m_gridHeights[0][i];
            }
        }

        double firstHeight = std::min(maxHeight, m_aspectRatio * m_averageMeshWidth);

        // Get true crossing splines heights
        int numLeftHeights = m_maxNumCenterSplineHeights;
        int numRightHeights = m_maxNumCenterSplineHeights;
        int numTrueCrossings = 0;
        for (int i = 0; i < m_numCrossingSplines[s]; ++i)
        {
            if (m_type[m_crossingSplinesIndexses[s][i]] != SplineTypes::crossing)
            {
                // true crossing splines only
                continue;
            }
            numTrueCrossings++;
            numLeftHeights = std::min(numLeftHeights, m_numCrossSplineLeftHeights[s][i]);
            numRightHeights = std::min(numRightHeights, m_numCrossSplineRightHeights[s][i]);
        }

        // no true cross splines: exponentially growing grid only
        if (numTrueCrossings == 0)
        {
            numLeftHeights = 0;
            numRightHeights = 0;
        }


        int startGridLineLeft = m_leftGridLineIndex[s];
        int endGridLineLeft = startGridLineLeft + m_numMSpline[s];
        int startGridLineRight = m_rightGridLineIndex[s];
        int endGridLineRight = startGridLineRight + m_numMSpline[s];

        double hh0LeftMaxRatio;
        double hh0RightMaxRatio;
        const int numIterations = 2;
        for (int iter = 0; iter < numIterations; ++iter)
        {
            ComputeVelocitiesSubIntervals(s, startGridLineLeft, endGridLineLeft, numLeftHeights, numRightHeights, firstHeight,
                m_leftGridLineIndex, m_rightGridLineIndex, numPerpendicularFacesOnSubintervalAndEdge, edgeVelocities, hh0LeftMaxRatio);

            ComputeVelocitiesSubIntervals(s, startGridLineRight, endGridLineRight, numLeftHeights, numRightHeights, firstHeight,
                m_rightGridLineIndex, m_leftGridLineIndex, numPerpendicularFacesOnSubintervalAndEdge, edgeVelocities, hh0RightMaxRatio);
        }

        // re-evaluate if growing grid outside is needed

        if (numLeftHeights == 0 && numRightHeights <= 1 ||
            numRightHeights == 0 && numLeftHeights <= 1 ||
            numLeftHeights == numRightHeights == 1)
        {
            m_growGridOutside = true;
        }

        // left part
        int numNLeftExponential = 0;
        if (m_growGridOutside)
        {
            numNLeftExponential = std::min(ComputeNumberExponentialIntervals(hh0LeftMaxRatio), m_maxNumN);
        }
        for (int i = startGridLineLeft; i < endGridLineLeft; ++i)
        {
            numPerpendicularFacesOnSubintervalAndEdge[1][i] = numNLeftExponential;
        }

        // right part
        int numNRightExponential = 0;
        if (m_growGridOutside)
        {
            numNRightExponential = std::min(ComputeNumberExponentialIntervals(hh0RightMaxRatio), m_maxNumN);
        }
        for (int i = startGridLineRight; i < endGridLineRight; ++i)
        {
            numPerpendicularFacesOnSubintervalAndEdge[1][i] = numNRightExponential;
        }
    }

    // compute local grow factors
    for (int s = 0; s < m_numSplines; s++)
    {
        if (m_numMSpline[s] < 1)
        {
            continue;
        }

        for (int i = m_leftGridLineIndex[s]; i < m_rightGridLineIndex[s] + m_numMSpline[s]; ++i)
        {
            if (m_gridLine[i].x == doubleMissingValue || m_gridLine[i + 1].x == doubleMissingValue || numPerpendicularFacesOnSubintervalAndEdge[1][i] < 1)
            {
                continue;
            }
            bool successfull = ComputeGrowFactor(m_gridHeights[1][i],
                edgeVelocities[i],
                numPerpendicularFacesOnSubintervalAndEdge[1][i],
                growFactorOnSubintervalAndEdge[1][i]);
            if (!successfull)
            {
                growFactorOnSubintervalAndEdge[1][i] = 1.0;
            }
        }
    }

    return true;
}

///comp_dgrow: this is another root finding algorithm, could go in the general part
bool GridGeom::Splines::ComputeGrowFactor(
    double totalGridHeight,
    double firstGridLayerHeight,
    double numberOfGridLayers,
    double& result)
{
    // eheight m_gridHeights
    double aspectRatioGrowFactor = 1.0;
    double heightDifference = ComputeTotalExponentialHeight(aspectRatioGrowFactor, firstGridLayerHeight, numberOfGridLayers) - totalGridHeight;

    double deps = 0.01;
    double aspectRatioGrowFactorIncremented = 1.0 + deps;
    double heightDifferenceIncremented = ComputeTotalExponentialHeight(aspectRatioGrowFactorIncremented, firstGridLayerHeight, numberOfGridLayers) - totalGridHeight;

    const double tolerance = 1e-8;
    const int numIterations = 1000;
    const double relaxationFactor = 0.5;
    double oldAspectRatio;
    double oldHeightDifference = heightDifference;

    if (std::abs(heightDifferenceIncremented) > tolerance&& std::abs(heightDifferenceIncremented - heightDifference) > tolerance)
    {
        for (int i = 0; i < numIterations; ++i)
        {
            oldAspectRatio = aspectRatioGrowFactor;
            oldHeightDifference = heightDifference;

            aspectRatioGrowFactor = aspectRatioGrowFactorIncremented;
            heightDifference = heightDifferenceIncremented;

            aspectRatioGrowFactorIncremented = aspectRatioGrowFactor - relaxationFactor * heightDifference / (heightDifference - oldHeightDifference + 1e-16) * (aspectRatioGrowFactor - oldAspectRatio);
            heightDifferenceIncremented = ComputeTotalExponentialHeight(aspectRatioGrowFactorIncremented, firstGridLayerHeight, numberOfGridLayers) - totalGridHeight;

            if (std::abs(oldHeightDifference) < tolerance)
            {
                break;
            }
        }
    }

    if (oldHeightDifference > tolerance)
    {
        result = doubleMissingValue;
        return false;
    }

    result = aspectRatioGrowFactorIncremented;
    return true;
}

double GridGeom::Splines::ComputeTotalExponentialHeight(double aspectRatioGrowFactor, double firstGridLayerHeights, int numberOfGridLayers)
{
    double height;
    if (std::abs(aspectRatioGrowFactor - 1.0) > 1e-8)
    {
        height = (std::pow(aspectRatioGrowFactor, numberOfGridLayers) - 1.0) / (aspectRatioGrowFactor - 1.0) * firstGridLayerHeights;
    }
    else
    {
        height = firstGridLayerHeights * numberOfGridLayers;
    }
    return height;
}

///comp_nfac
///compute the number of grid layers for a given grow factor, first grid layer height and total grid height
int GridGeom::Splines::ComputeNumberExponentialIntervals(const double hhMaxRatio)
{
    int numIntervals = 0;
    if (m_aspectRatioGrowFactor - 1.0 > 1e-8)
    {
        numIntervals = std::floor(std::log((m_aspectRatioGrowFactor - 1.0) * hhMaxRatio + 1.0) / log(m_aspectRatioGrowFactor));
    }
    else
    {
        numIntervals = std::floor(0.999 + hhMaxRatio);
    }
    return numIntervals;
}



bool GridGeom::Splines::ComputeVelocitiesSubIntervals(
    const int s,
    const int startGridLineIndex,
    const int endGridLineIndex,
    const int numHeights,
    const int numOtherSideHeights,
    const double firstHeight,
    const std::vector<int>& gridLineIndex,
    const std::vector<int>& otherGridLineIndex,
    std::vector<std::vector<int>>& numPerpendicularFacesOnSubintervalAndEdge,
    std::vector<double>& edgeVelocities,
    double& hh0MaxRatio)
{
    hh0MaxRatio = 0.0;
    if (numHeights > 1 && numHeights == numOtherSideHeights || numHeights > numOtherSideHeights)
    {
        double maxHeight = *std::max_element(m_gridHeights[0].begin() + startGridLineIndex, m_gridHeights[0].begin() + endGridLineIndex);

        int numNUniformPart = std::floor(maxHeight / firstHeight + 0.99999);
        numNUniformPart = std::min(numNUniformPart, m_maxNUniformPart);

        for (int i = startGridLineIndex; i < endGridLineIndex; ++i)
        {
            numPerpendicularFacesOnSubintervalAndEdge[0][i] = numNUniformPart;
            edgeVelocities[i] = m_gridHeights[0][i] / numNUniformPart;
            hh0MaxRatio = std::max(hh0MaxRatio, m_gridHeights[1][i] / edgeVelocities[i]);
        }
    }
    else
    {
        // only one subinterval: no uniform part
        int numNUniformPart = 0;
        for (int i = startGridLineIndex; i < endGridLineIndex; ++i)
        {
            numPerpendicularFacesOnSubintervalAndEdge[0][i] = numNUniformPart;
            edgeVelocities[i] = firstHeight;

            //compare with other side of spline
            int otherSideIndex = otherGridLineIndex[s] + m_numMSpline[s] - (i - gridLineIndex[s] + 1);

            if (edgeVelocities[otherSideIndex] != doubleMissingValue)
            {
                if (numPerpendicularFacesOnSubintervalAndEdge[0][otherSideIndex] == 0)
                {
                    edgeVelocities[i] = std::max(edgeVelocities[i], edgeVelocities[otherSideIndex]);
                }
                else
                {
                    edgeVelocities[i] = edgeVelocities[otherSideIndex];
                }
            }

            for (int j = 1; j < m_maxNumCenterSplineHeights; ++j)
            {
                m_gridHeights[j][i] = m_gridHeights[j - 1][i];
            }

            for (int j = startGridLineIndex; j < endGridLineIndex; ++j)
            {

                hh0MaxRatio = std::max(hh0MaxRatio, m_gridHeights[1][j] / edgeVelocities[j]);
            }
        }
    }

    return true;
}

/// compute the grid heights at grid edges on the center spline
///comp_gridheights
bool GridGeom::Splines::ComputeGridHeights()
{
    m_gridHeights.resize(m_maxNumCenterSplineHeights, std::vector<double>(m_numM - 1));
    std::fill(m_gridHeights.begin(), m_gridHeights.end(), std::vector<double>(m_numM - 1, doubleMissingValue));

    std::vector<std::vector<double>> heightsLeft(m_maxNumCenterSplineHeights, std::vector<double>(m_maxNumM, 0.0));
    std::vector<std::vector<double>> heightsRight(m_maxNumCenterSplineHeights, std::vector<double>(m_maxNumM, 0.0));
    std::vector<double> edgesCenterPoints(m_numM, 0.0);
    std::vector<double> crossingSplinesDimensionalCoordinates(m_numSplines, 0.0);
    std::vector<int> numHeightsLeft(m_numSplines, 0.0);
    std::vector<int> numHeightsRight(m_numSplines, 0.0);
    std::vector<double> localSplineDerivatives(m_numSplines, 0.0);
    std::vector<int> localValidSplineIndexes(m_numSplines, 0.0);

    for (int s = 0; s < m_numSplines; s++)
    {
        if (m_type[s] != SplineTypes::central)
        {
            continue;
        }

        int numM = m_numMSpline[s];

        // Get the minimum number of sub-intervals in the cross splines for this center spline
        int minNumLeftIntervals = *std::min_element(m_numCrossSplineLeftHeights[s].begin(), m_numCrossSplineLeftHeights[s].end());
        int minNumRightIntervals = *std::min_element(m_numCrossSplineRightHeights[s].begin(), m_numCrossSplineRightHeights[s].end());

        std::fill(heightsLeft[0].begin(), heightsLeft[0].begin() + numM, m_maximumGridHeights[s]);
        std::fill(heightsRight[0].begin(), heightsRight[0].begin() + numM, m_maximumGridHeights[s]);

        if (m_numCrossingSplines[s] == 1)
        {
            // only one crossing spline present: 
            for (int i = 0; i < minNumLeftIntervals; ++i)
            {
                std::fill(heightsLeft[i].begin(), heightsLeft[i].begin() + numM, m_crossSplineRightHeights[s][i][0]);

            }
            for (int i = 0; i < minNumRightIntervals; ++i)
            {
                std::fill(heightsRight[i].begin(), heightsRight[i].begin() + numM, m_crossSplineLeftHeights[s][i][0]);
            }
        }
        else
        {
            int leftGridLineIndex = m_leftGridLineIndex[s];
            edgesCenterPoints[0] = GetSplineLength(s, 0, m_gridLineDimensionalCoordinates[leftGridLineIndex]);
            for (int i = 0; i < numM; ++i)
            {
                edgesCenterPoints[i + 1] = edgesCenterPoints[i] + GetSplineLength(s, m_gridLineDimensionalCoordinates[leftGridLineIndex + i], m_gridLineDimensionalCoordinates[leftGridLineIndex + i + 1]);
            }

            // compute at edge center points
            for (int i = 0; i < numM; ++i)
            {
                edgesCenterPoints[i] = 0.5 * (edgesCenterPoints[i] + edgesCenterPoints[i + 1]);
            }
            edgesCenterPoints[numM] = doubleMissingValue;

            //compute center spline path length of cross splines
            crossingSplinesDimensionalCoordinates[0] = GetSplineLength(s, 0.0, m_crossSplineCoordinates[s][0]);
            for (int i = 0; i < m_numCrossingSplines[s] - 1; ++i)
            {
                crossingSplinesDimensionalCoordinates[i + 1] = crossingSplinesDimensionalCoordinates[i] + GetSplineLength(s, m_crossSplineCoordinates[s][i], m_crossSplineCoordinates[s][i + 1]);
            }

            for (int i = 0; i < m_numCrossingSplines[s]; ++i)
            {
                numHeightsLeft[i] = m_numCrossSplineLeftHeights[s][i];
                numHeightsRight[i] = m_numCrossSplineRightHeights[s][i];
            }


            for (int j = 0; j < m_maxNumCenterSplineHeights; ++j)
            {
                AddValueToVector(numHeightsLeft, -1);
                AddValueToVector(numHeightsRight, -1);

                bool success = FindNearestCrossSplines(s, j,
                    numHeightsLeft,
                    edgesCenterPoints,
                    m_crossSplineLeftHeights[s],
                    localValidSplineIndexes,
                    localSplineDerivatives,
                    crossingSplinesDimensionalCoordinates,
                    heightsLeft);

                if (!success)
                {
                    return false;
                }

                success = FindNearestCrossSplines(s, j,
                    numHeightsRight,
                    edgesCenterPoints,
                    m_crossSplineRightHeights[s],
                    localValidSplineIndexes,
                    localSplineDerivatives,
                    crossingSplinesDimensionalCoordinates,
                    heightsRight);

                if (!success)
                {
                    return false;
                }
            }
        }


        // store grid height
        for (int j = 0; j < m_maxNumCenterSplineHeights; ++j)
        {
            for (int i = 0; i < m_numMSpline[s]; ++i)
            {
                m_gridHeights[j][m_leftGridLineIndex[s] + i] = heightsLeft[j][i];
                m_gridHeights[j][m_rightGridLineIndex[s] + m_numMSpline[s] - i - 1] = heightsRight[j][i];
            }
        }
    }

    return true;
}

bool GridGeom::Splines::FindNearestCrossSplines(const int s,
    const int j,
    const std::vector<int>& numHeightsLeft,
    const std::vector<double>& edgesCenterPoints,
    const std::vector<std::vector<double>>& crossSplineLeftHeights,
    std::vector<int>& localValidSplineIndexes,
    std::vector<double>& localSplineDerivatives,
    std::vector<double>& crossingSplinesDimensionalCoordinates,
    std::vector<std::vector<double>>& heights)
{
    int numValid;
    GetValidSplineIndexses(s, m_numCrossingSplines[s], numHeightsLeft, localValidSplineIndexes, numValid);

    // no sub-heights to compute 
    if (numValid == 0)
    {
        return true;
    }

    int numM = m_numMSpline[s];
    std::vector<double> localCornerPoints(numValid);

    // TODO: strided memory access
    for (int i = 0; i < numValid; ++i)
    {
        int index = localValidSplineIndexes[i];
        localCornerPoints[i] = crossSplineLeftHeights[index][j];
    }

    SecondOrderDerivative(localCornerPoints, numValid, localSplineDerivatives);

    crossingSplinesDimensionalCoordinates[0] = GetSplineLength(s, 0.0, m_crossSplineCoordinates[s][0]);
    for (int i = 0; i < numM; ++i)
    {
        int leftIndex = 0;
        double leftCoordinate = crossingSplinesDimensionalCoordinates[localValidSplineIndexes[leftIndex]];
        int rightIndex = std::min(1, numValid - 1);
        double rightCoordinate = crossingSplinesDimensionalCoordinates[localValidSplineIndexes[rightIndex]];
        // Find two nearest cross splines
        while (rightCoordinate < edgesCenterPoints[i] && rightIndex < numValid)
        {
            leftIndex = rightIndex;
            leftCoordinate = rightCoordinate;
            rightIndex++;
            rightCoordinate = crossingSplinesDimensionalCoordinates[localValidSplineIndexes[rightIndex]];
            if (rightIndex == numValid - 1)
            {
                break;
            }
        }

        double factor = 0.0;
        if (std::abs(rightCoordinate - leftCoordinate) > 1e-8)
        {
            factor = (edgesCenterPoints[i] - leftCoordinate) / (rightCoordinate - leftCoordinate);
        }
        else
        {
            rightIndex = leftIndex;
        }

        factor = std::max(std::min(double(leftIndex + 1) + factor - 1.0, double(numValid - 1)), 0.0);

        bool success = Interpolate(localCornerPoints, localSplineDerivatives, factor, heights[j][i]);
        if (!success)
        {
            return false;
        }
    }

    return true;
}


// GetValidSplineIndexses
bool GridGeom::Splines::GetValidSplineIndexses(const int s, const int numValues, const std::vector<int>& v, std::vector<int>& validIndexses, int& numValid)
{
    numValid = 0;
    for (int i = 0; i < numValues; ++i)
    {
        if (v[i] >= 0)
        {
            validIndexses[numValid] = i;
            numValid++;
        }
    }
    return true;
};

/// SECT3R
/// compute the intersection of two splines
bool GridGeom::Splines::GetSplinesIntersection(const int first, const int second,
    const Projections& projection,
    double& crossProductIntersection,
    Point& intersectionPoint,
    double& firstSplineRatio,
    double& secondSplineRatio)
{
    double minimumCrossingDistance = std::numeric_limits<double>::max();
    double crossingDistance;
    int numCrossing = 0;
    double firstCrossingRatio;
    double secondCrossingRatio;
    int firstCrossingIndex = 0;
    int secondCrossingIndex = 0;
    Point closestIntersection;

    // First find a valid crossing, the closest to spline central point
    for (int n = 0; n < m_numSplineNodes[first] - 1; n++)
    {
        for (int nn = 0; nn < m_numSplineNodes[second] - 1; nn++)
        {
            Point intersection;
            double crossProduct;
            double firstRatio;
            double secondRatio;
            bool areCrossing = AreLinesCrossing(m_splineCornerPoints[first][n],
                m_splineCornerPoints[first][n + 1],
                m_splineCornerPoints[second][nn],
                m_splineCornerPoints[second][nn + 1],
                false,
                intersection,
                crossProduct,
                firstRatio,
                secondRatio,
                projection);


            if (areCrossing)
            {
                if (m_numSplineNodes[first] == 2)
                {
                    crossingDistance = std::min(minimumCrossingDistance, std::abs(firstRatio - 0.5));
                }
                else if (m_numSplineNodes[second] == 2)
                {
                    crossingDistance = std::abs(secondRatio - 0.5);
                }
                else
                {
                    crossingDistance = minimumCrossingDistance;
                }

                if (crossingDistance < minimumCrossingDistance || numCrossing == 0)
                {
                    minimumCrossingDistance = crossingDistance;
                    numCrossing = 1;
                    firstCrossingIndex = n;             //TI0
                    secondCrossingIndex = nn;           //TJ0
                    firstCrossingRatio = firstRatio;    //SL
                    secondCrossingRatio = secondRatio;  //SM
                }
            }
            closestIntersection = intersection;
        }
    }

    // if no crossing found, return
    if (numCrossing == 0)
    {
        return false;
    }

    double firstCrossing = firstCrossingRatio == -1 ? 0 : firstCrossingIndex + firstCrossingRatio;
    double secondCrossing = secondCrossingRatio == -1 ? 0 : secondCrossingIndex + secondCrossingRatio;

    // use bisection to find the intersection 
    double squaredDistanceBetweenCrossings = std::numeric_limits<double>::max();
    double maxSquaredDistanceBetweenCrossings = 1e-12;
    double maxDistanceBetweenVertices = 0.0001;
    double firstRatioIterations = 1.0;
    double secondRatioIterations = 1.0;
    double previousFirstCrossing;
    double previousSecondCrossing;
    int numIterations = 0;
    while (squaredDistanceBetweenCrossings > maxSquaredDistanceBetweenCrossings&& numIterations < 20)
    {
        // increment counter
        numIterations++;

        if (firstCrossingRatio > 0 && firstCrossingRatio < 1.0)
        {
            firstRatioIterations = 0.5 * firstRatioIterations;
        }
        if (secondCrossingRatio > 0 && secondCrossingRatio < 1.0)
        {
            secondRatioIterations = 0.5 * secondRatioIterations;
        }

        firstCrossing = std::max(0.0, std::min(firstCrossing, double(m_numSplineNodes[first])));
        secondCrossing = std::max(0.0, std::min(secondCrossing, double(m_numSplineNodes[second])));

        double firstLeft = std::max(0.0, std::min(double(m_numSplineNodes[first] - 1), firstCrossing - firstRatioIterations / 2.0));
        double firstRight = std::max(0.0, std::min(double(m_numSplineNodes[first] - 1), firstCrossing + firstRatioIterations / 2.0));

        double secondLeft = std::max(0.0, std::min(double(m_numSplineNodes[second] - 1), secondCrossing - secondRatioIterations / 2.0));
        double secondRight = std::max(0.0, std::min(double(m_numSplineNodes[second] - 1), secondCrossing + secondRatioIterations / 2.0));

        firstRatioIterations = firstRight - firstLeft;
        secondRatioIterations = secondRight - secondLeft;

        Point firstLeftSplinePoint;
        Interpolate(m_splineCornerPoints[first], m_splineDerivatives[first], firstLeft, firstLeftSplinePoint);
        Point firstRightSplinePoint;
        Interpolate(m_splineCornerPoints[first], m_splineDerivatives[first], firstRight, firstRightSplinePoint);

        Point secondLeftSplinePoint;
        Interpolate(m_splineCornerPoints[second], m_splineDerivatives[second], secondLeft, secondLeftSplinePoint);
        Point secondRightSplinePoint;
        Interpolate(m_splineCornerPoints[second], m_splineDerivatives[second], secondRight, secondRightSplinePoint);

        Point oldIntersection = closestIntersection;

        double crossProduct;
        double firstRatio;
        double secondRatio;
        bool areCrossing = AreLinesCrossing(firstLeftSplinePoint, firstRightSplinePoint,
            secondLeftSplinePoint, secondRightSplinePoint,
            true,
            closestIntersection,
            crossProduct,
            firstRatio,
            secondRatio,
            projection);

        // search close by
        if (-2.0 < firstRatio < 3.0 && -2.0 < secondRatio < 3.0)
        {
            previousFirstCrossing = firstCrossing;
            previousSecondCrossing = secondCrossing;

            firstCrossing = firstLeft + firstRatio * (firstRight - firstLeft);
            secondCrossing = secondLeft + secondRatio * (secondRight - secondLeft);

            firstCrossing = std::max(0.0, std::min(m_numSplineNodes[first] - 1.0, firstCrossing));
            secondCrossing = std::max(0.0, std::min(m_numSplineNodes[second] - 1.0, secondCrossing));

            if (areCrossing)
            {
                numCrossing = 1;
                crossProductIntersection = crossProduct;
            }

            if (std::abs(firstCrossing - previousFirstCrossing) > maxDistanceBetweenVertices ||
                std::abs(secondCrossing - previousSecondCrossing) > maxDistanceBetweenVertices)
            {
                squaredDistanceBetweenCrossings = ComputeSquaredDistance(oldIntersection, closestIntersection, projection);
            }
            else
            {
                break;
            }
        }
    }

    if (numCrossing == 1)
    {
        intersectionPoint = closestIntersection;
        firstSplineRatio = firstCrossing;
        secondSplineRatio = secondCrossing;
        return true;
    }

    //not crossing
    return false;
}

/// get_crosssplines
/// compute the intersection of two splines, one must have only two nodes
bool GridGeom::Splines::GetSplineIntersections(const int index)
{
    m_numCrossingSplines[index] = 0;
    std::fill(m_crossingSplinesIndexses[index].begin(), m_crossingSplinesIndexses[index].end(), -1);
    std::fill(m_isLeftOriented[index].begin(), m_isLeftOriented[index].end(), true);
    std::fill(m_crossSplineCoordinates[index].begin(), m_crossSplineCoordinates[index].end(), std::numeric_limits<double>::max());
    std::fill(m_cosCrossingAngle[index].begin(), m_cosCrossingAngle[index].end(), doubleMissingValue);

    for (int s = 0; s < m_numSplines; ++s)
    {
        // a crossing is a spline with 2 nodes and another with more than 2 nodes
        if (m_numSplineNodes[s] == 2 && m_numSplineNodes[index] == 2 ||
            m_numSplineNodes[s] > 2 && m_numSplineNodes[index] > 2)
        {
            continue;
        }

        double crossProductIntersection;
        Point intersectionPoint;
        double firstSplineRatio;
        double secondSplineRatio;
        bool crossing = GetSplinesIntersection(index, s, m_projection, crossProductIntersection, intersectionPoint, firstSplineRatio, secondSplineRatio);

        if (std::abs(crossProductIntersection) < m_dtolcos)
        {
            crossing = false;
        }


        if (crossing)
        {
            m_numCrossingSplines[index]++;
            m_crossingSplinesIndexses[index][s] = s;
            if (crossProductIntersection > 0.0)
            {
                m_isLeftOriented[index][s] = false;
            }
            m_crossSplineCoordinates[index][s] = firstSplineRatio;
            m_cosCrossingAngle[index][s] = crossProductIntersection;
        }
    }

    const auto sortedIndexses = SortedIndexes(m_crossSplineCoordinates[index]);
    ReorderVector(m_crossSplineCoordinates[index], sortedIndexses);
    ReorderVector(m_crossingSplinesIndexses[index], sortedIndexses);
    ReorderVector(m_isLeftOriented[index], sortedIndexses);

    return true;
}


///make_wholegridline
bool GridGeom::Splines::MakeAllGridLines(bool isSpacingCurvatureAdapeted)
{

    int numCenterSplines = 0;
    for (int s = 0; s < m_numSplines; ++s)
    {
        //center splines only
        if (m_type[s] != SplineTypes::central)
        {
            continue;
        }
        numCenterSplines += 1;
    }

    if (numCenterSplines == 0)
    {
        return false;
    }


    int gridLineIndex = 0;
    for (int s = 0; s < m_numSplines; ++s)
    {
        //center splines only
        if (m_type[s] != SplineTypes::central)
        {
            continue;
        }

        // upper bound of m_gridLine, with two sides of spline and two missing values added
        int sizeGridLine = gridLineIndex + 1 + 2 * (m_maxNumM + 1) + 2;
        // increase size
        ResizeVectorIfNeeded(sizeGridLine, m_gridLine, { doubleMissingValue,doubleMissingValue });
        m_gridLineDimensionalCoordinates.resize(sizeGridLine);

        if (gridLineIndex > 0)
        {
            gridLineIndex++;
            m_gridLine[gridLineIndex] = { doubleMissingValue,doubleMissingValue };
            m_gridLineDimensionalCoordinates[gridLineIndex] = doubleMissingValue;
        }

        m_leftGridLineIndex[s] = gridLineIndex;

        int numM = 0;
        bool success = MakeGridLine(s, gridLineIndex, m_gridLine, m_gridLineDimensionalCoordinates, numM);
        if (!success)
        {
            return false;
        }

        gridLineIndex = gridLineIndex + numM + 1;
        m_gridLine[gridLineIndex] = Point{ doubleMissingValue,  doubleMissingValue };
        m_gridLineDimensionalCoordinates[gridLineIndex] = doubleMissingValue;
        gridLineIndex++;

        //add other side of gridline
        m_rightGridLineIndex[s] = gridLineIndex;
        for (int i = m_rightGridLineIndex[s] - 1, j = m_rightGridLineIndex[s] - 1; j >= m_leftGridLineIndex[s]; ++i, --j)
        {
            m_gridLine[i] = m_gridLine[j];
            m_gridLineDimensionalCoordinates[i] = m_gridLineDimensionalCoordinates[j];
        }

        //compute new (actual) grid size
        //new size   old size   both sides of spline   DMISS between both sides
        gridLineIndex = gridLineIndex + numM + 1;

        m_gridLine[gridLineIndex] = Point{ doubleMissingValue,  doubleMissingValue };
        m_gridLineDimensionalCoordinates[gridLineIndex] = doubleMissingValue;
        m_numMSpline[s] = numM;
        m_numM = gridLineIndex;
    }

    return true;
}

/// make_gridline, generate a gridline on a spline with a prescribed maximum mesh width
/// generate a gridline on a spline with a prescribed maximum mesh width
bool GridGeom::Splines::MakeGridLine(int splineIndex,
    int startingIndex,
    std::vector<Point>& gridLine,
    std::vector<double>& adimensionalCoordinates,
    int& numM)
{
    // first estimation of nodes along m
    numM = 1 + std::floor(m_splinesLength[splineIndex] / m_averageMeshWidth);
    numM = std::min(numM, m_maxNumM);

    double endSplineAdimensionalCoordinate = m_numSplineNodes[splineIndex] - 1;
    double splineLength = GetSplineLength(splineIndex, 0.0, endSplineAdimensionalCoordinate, 10, m_isSpacingCurvatureAdapted, m_maximumGridHeights[splineIndex]);

    gridLine[startingIndex] = m_splineCornerPoints[splineIndex][0];
    FuncDimensionalToAdimensionalDistance func(*this, splineIndex, m_isSpacingCurvatureAdapted, m_maximumGridHeights[splineIndex]);

    double currentMaxWidth = std::numeric_limits<double>::max();
    while (currentMaxWidth > m_averageMeshWidth)
    {
        currentMaxWidth = 0.0;
        for (int n = 1; n <= numM; ++n)
        {
            int index = startingIndex + n;
            double dimensionalDistance = splineLength * double(n) / double(numM);
            func.SetDimensionalDistance(dimensionalDistance);
            adimensionalCoordinates[index] = FindFunctionRootWithGoldenSectionSearch(func, 0, endSplineAdimensionalCoordinate);
            Interpolate(m_splineCornerPoints[splineIndex], m_splineDerivatives[splineIndex], adimensionalCoordinates[index], gridLine[index]);
            currentMaxWidth = std::max(currentMaxWidth, Distance(gridLine[index - 1], gridLine[index], m_projection));
        }

        // a gridline is computed
        if (currentMaxWidth < m_averageMeshWidth || numM == m_maxNumM)
        {
            break;
        }

        // room for sub-division
        if (currentMaxWidth > m_averageMeshWidth)
        {
            numM = std::min(std::max(int(m_maxNumM / m_maximumGridHeights[splineIndex] * numM), numM + 1), m_maxNumM);
        }
    }
    return true;
}

///get_splineprops
bool GridGeom::Splines::ComputeSplineProperties(const bool restoreOriginalProperties)
{
    // cross spline is a spline with only two points, others are non-cross splines
    bool successful = AllocateSplinesProperties();
    if (!successful)
    {
        return false;
    }

    for (int s = 0; s < m_numSplines; ++s)
    {
        successful = GetSplineIntersections(s);

        if (!successful)
        {
            return false;
        }
    }
    // select all non-cross splines only
    for (int s = 0; s < m_numSplines; ++s)
    {
        m_type[s] = SplineTypes::crossing;
        // select all non-cross splines for growing the grid
        if (m_numSplineNodes[s] > 2)
        {
            m_type[s] = SplineTypes::central;
        }
    }
    // check the cross splines. The center spline is the middle spline that crosses the cross spline
    for (int s = 0; s < m_numSplines; ++s)
    {
        // only crossing splines with one or more center spline
        if (m_numSplineNodes[s] != 2 || m_numCrossingSplines[s] < 1)
        {
            continue;
        }

        int middleCrossingSpline = std::min(m_numCrossingSplines[s] / 2, m_numCrossingSplines[s]);
        int crossingSplineIndex = m_crossingSplinesIndexses[s][middleCrossingSpline];

        // if m_numIntersectingSplines[s] is even, check if the middle spline has already been assigned as a bounding spline
        if (m_type[crossingSplineIndex] != SplineTypes::central && 2 * crossingSplineIndex == m_numCrossingSplines[s])
        {
            middleCrossingSpline = std::min(middleCrossingSpline + 1, m_numCrossingSplines[s] - 1);
            crossingSplineIndex = m_crossingSplinesIndexses[s][middleCrossingSpline];
        }

        if (m_type[crossingSplineIndex] == SplineTypes::central)
        {
            // associate bounding splines with the middle spline
            for (int i = 0; i < middleCrossingSpline; ++i)
            {
                int index = m_crossingSplinesIndexses[s][i];
                m_type[index] = SplineTypes::lateral; // lateral spline
                m_centralSplineIndex[index] = -crossingSplineIndex;

            }
            for (int i = middleCrossingSpline + 1; i < m_numCrossingSplines[s]; ++i)
            {
                int index = m_crossingSplinesIndexses[s][i];
                m_type[index] = SplineTypes::lateral; // lateral spline
                m_centralSplineIndex[index] = -crossingSplineIndex;
            }
        }
    }

    if (restoreOriginalProperties)
    {
        // restore original spline properties
        for (int s = 0; s < m_numOriginalSplines; ++s)
        {
            m_leftGridLineIndex[s] = m_leftGridLineIndexOriginal[s];
            m_rightGridLineIndex[s] = m_rightGridLineIndexOriginal[s];
            m_numMSpline[s] = m_mfacOriginal[s];
            m_maximumGridHeights[s] = m_maximumGridHeightsOriginal[s];
            m_type[s] = m_originalTypes[s];
        }

        //mark new splines as artificial cross splines
        for (int s = m_numOriginalSplines; s < m_numSplines; ++s)
        {
            m_type[s] = SplineTypes::arficial;
        }
    }

    bool successfull = ComputeHeights();
    return successfull;
}

/// get the grid heights from the cross spline information
/// get_heights
bool GridGeom::Splines::ComputeHeights()
{
    for (int i = 0; i < m_numSplines; ++i)
    {
        // Heights should be computed only for center splines
        if (m_numSplineNodes[i] <= 2)
        {
            continue;
        }
        for (int j = 0; j < m_numCrossingSplines[i]; ++j)
        {
            int intersectingSplineIndex = m_crossingSplinesIndexses[i][j];
            bool success = ComputeSubHeights(i, j);
            if (!success)
            {
                return false;
            }
        }
    }

    // compute m_maximumGridHeight
    for (int s = 0; s < m_numSplines; ++s)
    {
        if (m_numCrossingSplines[s] == 0)
        {
            m_maximumGridHeights[s] = m_aspectRatioFirstLayer * m_splinesLength[s];
            continue;
        }
        double maximumHeight = 0.0;
        for (int c = 0; c < m_numCrossingSplines[s]; ++c)
        {
            double sumLeftHeights = 0.0;
            for (int ss = 0; ss < m_numCrossSplineLeftHeights[s][c]; ++ss)
            {
                sumLeftHeights += m_crossSplineLeftHeights[s][c][ss];
            }
            double sumRightHeights = 0.0;
            for (int ss = 0; ss < m_numCrossSplineRightHeights[s][c]; ++ss)
            {
                sumRightHeights += m_crossSplineRightHeights[s][c][ss];
            }
            maximumHeight = std::max(maximumHeight, std::max(sumLeftHeights, sumRightHeights));
        }

        m_maximumGridHeights[s] = maximumHeight;
    }
    return true;
}

///comp_subheights, compute the height of the subintervals of grid layers on a cross spline, w.r.t. a center spline
bool GridGeom::Splines::ComputeSubHeights( int centerSplineIndex,  int crossingSplineLocalIndex)
{
    // find center spline index
    int centerSplineLocalIndex = 0;
    int crossingSplineIndex = m_crossingSplinesIndexses[centerSplineIndex][crossingSplineLocalIndex]; //js
    //m_intersectingSplinesIndexses[intersectingSplinesIndex] // ics
    for (int s = 0; s < m_numCrossingSplines[crossingSplineIndex]; ++s)
    {
        if (m_crossingSplinesIndexses[crossingSplineIndex][s] == centerSplineIndex)
        {
            centerSplineLocalIndex = s;
            break;
        }
    }

    // right part
    int numSubIntervalsRight = 0;
    int rightCenterSplineIndex = centerSplineLocalIndex;
    int leftCenterSplineIndex;
    m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex].resize(m_maxNumCenterSplineHeights, 0);
    for (int s = centerSplineLocalIndex; s < m_numCrossingSplines[crossingSplineIndex] - 1; ++s)
    {
        if (numSubIntervalsRight >= m_maxNumCenterSplineHeights)
        {
            break;
        }
        if (m_centralSplineIndex[m_crossingSplinesIndexses[crossingSplineIndex][s + 1]] != -centerSplineIndex)
        {
            continue;
        }
        leftCenterSplineIndex = rightCenterSplineIndex;
        rightCenterSplineIndex = s + 1;
        m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex][numSubIntervalsRight] = GetSplineLength(crossingSplineIndex,
            m_crossSplineCoordinates[crossingSplineIndex][leftCenterSplineIndex], m_crossSplineCoordinates[crossingSplineIndex][rightCenterSplineIndex]);
        numSubIntervalsRight++;
    }

    m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex][numSubIntervalsRight] = GetSplineLength(crossingSplineIndex,
        m_crossSplineCoordinates[crossingSplineIndex][rightCenterSplineIndex], m_numSplineNodes[crossingSplineIndex] - 1);

    numSubIntervalsRight++;
    std::fill(m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex].begin() + numSubIntervalsRight, m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex].end(), 0.0);
    m_numCrossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex] = numSubIntervalsRight;

    // left part
    int numSubIntervalsLeft = 0;
    leftCenterSplineIndex = centerSplineLocalIndex;
    m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex].resize(m_maxNumCenterSplineHeights, 0);
    for (int s = centerSplineLocalIndex; s >= 1; --s)
    {
        if (numSubIntervalsLeft >= m_maxNumCenterSplineHeights)
        {
            break;
        }
        if (m_centralSplineIndex[m_crossingSplinesIndexses[crossingSplineIndex][s - 1]] != -centerSplineIndex)
        {
            continue;
        }
        rightCenterSplineIndex = leftCenterSplineIndex;
        leftCenterSplineIndex = s - 1;
        m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex][numSubIntervalsLeft] = GetSplineLength(crossingSplineIndex,
            m_crossSplineCoordinates[crossingSplineIndex][leftCenterSplineIndex], m_crossSplineCoordinates[crossingSplineIndex][rightCenterSplineIndex]);
        numSubIntervalsLeft++;
    }

    m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex][numSubIntervalsLeft] = GetSplineLength(crossingSplineIndex,
        0.0, m_crossSplineCoordinates[crossingSplineIndex][leftCenterSplineIndex]);

    numSubIntervalsLeft++;
    std::fill(m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex].begin() + numSubIntervalsLeft, m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex].end(), 0.0);
    m_numCrossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex] = numSubIntervalsLeft;

    // if not left oriented, swap
    if (!m_isLeftOriented[centerSplineIndex][centerSplineLocalIndex])
    {
        m_numCrossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex] = numSubIntervalsRight;
        m_numCrossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex] = numSubIntervalsLeft;

        std::vector<double> leftSubIntervalsTemp(m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex]);
        m_crossSplineLeftHeights[centerSplineIndex][crossingSplineLocalIndex] = m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex];
        m_crossSplineRightHeights[centerSplineIndex][crossingSplineLocalIndex] = leftSubIntervalsTemp;
    }

    return true;
}

///splinelength
/// TODO: remove special treatment assign delta and calculate number of points 
double GridGeom::Splines::GetSplineLength(int index,
    double beginFactor,
    double endFactor,
    int numSamples,
    bool accountForCurvature,
    double height,
    double assignedDelta)
{
    double delta = assignedDelta;
    int numPoints = endFactor / delta + 1;
    if (delta < 0.0)
    {
        delta = 1.0 / numSamples;
        numPoints = std::max(std::floor(0.9999 + (endFactor - beginFactor) / delta), 10.0);
        delta = (endFactor - beginFactor) / numPoints;
    }

    // first point
    Point leftPoint;
    Interpolate(m_splineCornerPoints[index], m_splineDerivatives[index], beginFactor, leftPoint);

    double splineLength = 0.0;

    double rightPointCoordinateOnSpline = beginFactor;
    double leftPointCoordinateOnSpline;
    for (int p = 0; p < numPoints; ++p)
    {
        leftPointCoordinateOnSpline = rightPointCoordinateOnSpline;
        rightPointCoordinateOnSpline += delta;
        if (rightPointCoordinateOnSpline > endFactor)
        {
            rightPointCoordinateOnSpline = endFactor;
        }

        Point rightPoint;
        Interpolate(m_splineCornerPoints[index], m_splineDerivatives[index], rightPointCoordinateOnSpline, rightPoint);
        double curvatureFactor = 0.0;
        if (accountForCurvature)
        {
            Point normalVector;
            Point tangentialVector;
            ComputeCurvatureOnSplinePoint(index, 0.5 * (rightPointCoordinateOnSpline + leftPointCoordinateOnSpline), curvatureFactor, normalVector, tangentialVector);
        }
        splineLength = splineLength + Distance(leftPoint, rightPoint, m_projection) * (1.0 + curvatureFactor * height);
        leftPoint = rightPoint;
    }

    return splineLength;
}

/// comp_curv
/// compute curvature in a point on a spline
bool GridGeom::Splines::ComputeCurvatureOnSplinePoint(
    int splineIndex,
    double adimensionalPointCoordinate,
    double& curvatureFactor,
    Point& normalVector,
    Point& tangentialVector)
{
    auto const leftCornerPoint = int(std::max(std::min(double(std::floor(adimensionalPointCoordinate)), double(m_numSplineNodes[splineIndex] - 1)), 0.0));
    auto const rightCornerPoint = int(std::max(double(leftCornerPoint + 1.0), 0.0));

    double leftSegment = rightCornerPoint - adimensionalPointCoordinate;
    double rightSegment = adimensionalPointCoordinate - leftCornerPoint;

    Point pointCoordinate;
    Interpolate(m_splineCornerPoints[splineIndex], m_splineDerivatives[splineIndex], adimensionalPointCoordinate, pointCoordinate);

    Point p = m_splineCornerPoints[splineIndex][rightCornerPoint] - m_splineCornerPoints[splineIndex][leftCornerPoint]
        + (m_splineDerivatives[splineIndex][leftCornerPoint] * (-3.0 * leftSegment * leftSegment + 1.0)
            + m_splineDerivatives[splineIndex][rightCornerPoint] * (3.0 * rightSegment * rightSegment - 1.0)) / 6.0;

    Point pp = m_splineDerivatives[splineIndex][leftCornerPoint] * leftSegment + m_splineDerivatives[splineIndex][rightCornerPoint] * rightSegment;

    if (m_projection == Projections::spherical)
    {
        p.TransformSphericalToCartesian(pointCoordinate.y);
        pp.TransformSphericalToCartesian(pointCoordinate.y);
    }

    curvatureFactor = std::abs(pp.x * p.y - pp.y * p.x) / std::pow((p.x * p.x + p.y * p.y + 1e-8), 1.5);

    Point incremenetedPointCoordinate = pointCoordinate + p * 1e-4;
    NormalVectorOutside(pointCoordinate, incremenetedPointCoordinate, normalVector, m_projection);

    double distance = Distance(pointCoordinate, incremenetedPointCoordinate, m_projection);
    double dx = GetDx(pointCoordinate, incremenetedPointCoordinate, m_projection);
    double dy = GetDy(pointCoordinate, incremenetedPointCoordinate, m_projection);

    tangentialVector.x = dx / distance;
    tangentialVector.y = dy / distance;

    return true;
}


/// second order derivative of spline coordinates
bool GridGeom::Splines::SecondOrderDerivative(const std::vector<Point>& coordinates, int numNodes, std::vector<Point>& coordinatesDerivatives)
{
    std::vector<Point> u(numNodes);
    u[0] = { 0.0, 0.0 };
    coordinatesDerivatives.resize(coordinates.size(), { 0.0, 0.0 });
    coordinatesDerivatives[0] = { 0.0, 0.0 };

    for (int i = 1; i < numNodes - 1; i++)
    {
        const Point p = coordinatesDerivatives[i - 1] * 0.5 + 2.0;
        coordinatesDerivatives[i].x = -0.5 / p.x;
        coordinatesDerivatives[i].y = -0.5 / p.y;

        const Point delta = coordinates[i + 1] - coordinates[i] - (coordinates[i] - coordinates[i - 1]);
        u[i] = (delta * 6.0 / 2.0 - u[i - 1] * 0.5) / p;
    }

    coordinatesDerivatives[numNodes - 1] = { 0.0, 0.0 };
    for (int i = numNodes - 2; i >= 0; i--)
    {
        coordinatesDerivatives[i] = coordinatesDerivatives[i] * coordinatesDerivatives[i + 1] + u[i];
    }

    return true;
}

bool GridGeom::Splines::SecondOrderDerivative(const std::vector<double>& coordinates, int numNodes, std::vector<double>& coordinatesDerivatives)
{
    std::vector<double> u(numNodes);
    u[0] = 0.0;
    //coordinatesDerivatives.resize(coordinates.size());
    coordinatesDerivatives[0] = 0.0;

    for (int i = 1; i < numNodes - 1; i++)
    {
        const double p = coordinatesDerivatives[i - 1] * 0.5 + 2.0;
        coordinatesDerivatives[i] = -0.5 / p;

        const double delta = coordinates[i + 1] - coordinates[i] - (coordinates[i] - coordinates[i - 1]);
        u[i] = (delta * 6.0 / 2.0 - u[i - 1] * 0.5) / p;
    }

    coordinatesDerivatives[numNodes - 1] = 0.0;
    for (int i = numNodes - 2; i >= 0; i--)
    {
        coordinatesDerivatives[i] = coordinatesDerivatives[i] * coordinatesDerivatives[i + 1] + u[i];
    }

    return true;
}

/// splint
template<typename T>
bool GridGeom::Splines::Interpolate(const std::vector<T>& coordinates, const std::vector<T>& coordinatesDerivatives, double pointAdimensionalCoordinate, T& pointCoordinate)
{
    if (pointAdimensionalCoordinate < 0)
    {
        return false;
    }

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

    pointCoordinate = coordinates[low] * a + coordinates[high] * b +
        (coordinatesDerivatives[low] * (pow(a, 3) - a) + coordinatesDerivatives[high] * (pow(b, 3) - b)) / 6.0 * splFac;

    return true;
}