//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2021.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// contact: delft3d.support@deltares.nl
// Stichting Deltares
// P.O. Box 177
// 2600 MH Delft, The Netherlands
//
// All indications and logos of, and references to, "Delft3D" and "Deltares"
// are registered trademarks of Stichting Deltares, and remain the property of
// Stichting Deltares. All rights reserved.
//
//------------------------------------------------------------------------------

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridFromSplines.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Parameters.hpp>
#include <MeshKernel/SplineAlgorithms.hpp>
#include <MeshKernel/Splines.hpp>

#include <span>

// How to get the constants?

// extern const int MaxDegree;
// extern const int MaxDegreeP1;

static constexpr int MaxDegree = 4;
static const int MaxDegreeP1 = MaxDegree + 1;

extern "C"
{
    /// @brief Function of the rpoly library
    // void rpoly_ak1(const std::array<double, MaxDegreeP1>& op,
    //                int& degree,
    //                std::array<double, MaxDegree>& zeror,
    //                std::array<double, MaxDegree>& zeroi);
    void rpoly(double op[MaxDegreeP1],
               int* degree,
               double zeror[MaxDegree],
               double zeroi[MaxDegree]);
}

namespace meshkernel
{
    void RootFinder(const std::array<double, MaxDegreeP1>& coefficients,
                    int& degree,
                    std::array<double, MaxDegree>& realRoots,
                    std::array<double, MaxDegree>& imaginaryRoots)
    {
        double coeffs[MaxDegreeP1];
        double rootsR[MaxDegree];
        double rootsI[MaxDegree];

        for (size_t i = 0; i < MaxDegreeP1; ++i)
        {
            coeffs[i] = coefficients[i];
        }

        for (size_t i = 0; i < MaxDegree; ++i)
        {
            rootsR[i] = realRoots[i];
            rootsI[i] = imaginaryRoots[i];
        }

        rpoly(coeffs, &degree, rootsR, rootsI);

        for (size_t i = 0; i < MaxDegree; ++i)
        {
            realRoots[i] = rootsR[i];
            imaginaryRoots[i] = rootsI[i];
        }
    }

    CurvilinearGridFromSplines::CurvilinearGridFromSplines(std::shared_ptr<Splines> splines,
                                                           const CurvilinearParameters& curvilinearParameters,
                                                           const SplinesToCurvilinearParameters& splinesToCurvilinearParameters)
        : m_splines(splines),
          m_curvilinearParameters(curvilinearParameters),
          m_splinesToCurvilinearParameters(splinesToCurvilinearParameters)
    {
        CheckCurvilinearParameters(curvilinearParameters);
        CheckSplinesToCurvilinearParameters(splinesToCurvilinearParameters);

        for (size_t i = 0; i < splines->m_splineNodes.size(); ++i)
        {
            std::cout << "separation distance of node for spline: " << i << " -- ";

            for (size_t j = 1; j < splines->m_splineNodes[i].size(); ++j)
            {
                double distance = ComputeDistance(splines->m_splineNodes[i][j - 1], splines->m_splineNodes[i][j], splines->m_projection);
                std::cout << distance << "  ";
            }

            std::cout << std::endl;
        }

        m_onTopOfEachOtherSquaredTolerance = m_splinesToCurvilinearParameters.nodes_on_top_of_each_other_tolerance *
                                             m_splinesToCurvilinearParameters.nodes_on_top_of_each_other_tolerance;
    }

    /// to be called after all splines have been stored
    void CurvilinearGridFromSplines::AllocateSplinesProperties()
    {
        const auto numSplines = m_splines->GetNumSplines();
        m_type.resize(numSplines);

        m_centralSplineIndex.resize(numSplines);
        std::fill(m_centralSplineIndex.begin(), m_centralSplineIndex.end(), constants::missing::intValue);

        m_numCrossingSplines.resize(numSplines, 0);
        std::fill(m_numCrossingSplines.begin(), m_numCrossingSplines.end(), 0);

        m_maximumGridHeights.resize(numSplines, -999.0);
        std::fill(m_maximumGridHeights.begin(), m_maximumGridHeights.end(), constants::missing::doubleValue);

        lin_alg::ResizeAndFillMatrix(m_crossingSplinesIndices, numSplines, numSplines, false, 0U);
        lin_alg::ResizeAndFillMatrix(m_isLeftOriented, numSplines, numSplines, false, true);
        lin_alg::ResizeAndFillMatrix(m_crossSplineCoordinates, numSplines, numSplines, false, constants::missing::doubleValue);
        lin_alg::ResizeAndFillMatrix(m_cosCrossingAngle, numSplines, numSplines, false, constants::missing::doubleValue);
        lin_alg::ResizeAndFillMatrix(m_crossSplineLeftHeights, numSplines, numSplines, false,
                                     std::vector<double>(m_maxNumCenterSplineHeights, constants::missing::doubleValue));
        lin_alg::ResizeAndFillMatrix(m_crossSplineRightHeights, numSplines, numSplines, false,
                                     std::vector<double>(m_maxNumCenterSplineHeights, constants::missing::doubleValue));
        lin_alg::ResizeAndFillMatrix(m_numCrossSplineLeftHeights, numSplines, numSplines, false, 0U);
        lin_alg::ResizeAndFillMatrix(m_numCrossSplineRightHeights, numSplines, numSplines, false, 0U);

        m_numMSplines.resize(numSplines);
        std::fill(m_numMSplines.begin(), m_numMSplines.end(), 0);
        m_leftGridLineIndex.resize(numSplines);
        std::fill(m_leftGridLineIndex.begin(), m_leftGridLineIndex.end(), constants::missing::uintValue);
        m_rightGridLineIndex.resize(numSplines);
        std::fill(m_rightGridLineIndex.begin(), m_rightGridLineIndex.end(), constants::missing::uintValue);
    }

    std::unique_ptr<CurvilinearGrid> CurvilinearGridFromSplines::Compute()
    {
        Initialize();

        // Grow grid, from the second layer
        // for (auto layer = 1; layer <= 1; ++layer)
        for (auto layer = 1; layer <= m_curvilinearParameters.n_refinement; ++layer)
        {
            Iterate(layer);
        }

        const auto deleteSkinnyTriangles = m_splinesToCurvilinearParameters.remove_skinny_triangles == 1;

        if (deleteSkinnyTriangles)
        {
            DeleteSkinnyTriangles();
        }

        return ComputeCurvilinearGridFromGridPoints();
    }

    void CurvilinearGridFromSplines::DeleteSkinnyTriangles()
    {
        constexpr UInt numMaxIterations = 10;
        const UInt numN = static_cast<UInt>(m_gridPoints.rows()) - 2;
        constexpr double squaredDistanceTolerance = 1e-4;
        constexpr double cosineTolerance = 1e-2;
        constexpr double maxCosine = 0.93969;
        for (UInt j = numN - 1; j >= 1; --j)
        {
            for (UInt iter = 0; iter < numMaxIterations; ++iter)
            {
                UInt numChanged = 0;
                UInt firstRightIndex = 0;
                UInt i = 0;

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

                    if (!m_gridPoints(j, i).IsValid())
                    {
                        continue;
                    }

                    auto [firstLeftIndex, firstRightIndex] = GetNeighbours(m_gridPoints.row(j), i);

                    const auto squaredRightDistance = ComputeSquaredDistance(m_gridPoints(j, i),
                                                                             m_gridPoints(j, firstRightIndex),
                                                                             m_splines->m_projection);

                    if (squaredRightDistance < squaredDistanceTolerance)
                    {
                        continue;
                    }

                    // Detect triangular cell
                    if (!m_gridPoints(j + 1, i).IsValid())
                    {
                        continue;
                    }

                    const auto squaredLeftDistance = ComputeSquaredDistance(m_gridPoints(j, firstLeftIndex),
                                                                            m_gridPoints(j, i),
                                                                            m_splines->m_projection);
                    if (squaredLeftDistance < squaredDistanceTolerance)
                    {
                        firstLeftIndex = i;
                    }

                    if (m_gridPoints(j + 1, firstRightIndex).IsValid())
                    {
                        const auto squaredCurrentDistance = ComputeSquaredDistance(m_gridPoints(j + 1, i),
                                                                                   m_gridPoints(j + 1, firstRightIndex),
                                                                                   m_splines->m_projection);
                        const auto currentCosPhi = NormalizedInnerProductTwoSegments(
                            m_gridPoints(j + 1, i),
                            m_gridPoints(j, i),
                            m_gridPoints(j + 1, i),
                            m_gridPoints(j, firstRightIndex),
                            m_splines->m_projection);
                        if (squaredCurrentDistance < squaredDistanceTolerance && currentCosPhi > maxCosine)
                        {

                            // determine persistent node
                            const auto leftCosPhi = NormalizedInnerProductTwoSegments(
                                m_gridPoints(j - 1, i),
                                m_gridPoints(j, i),
                                m_gridPoints(j, i),
                                m_gridPoints(j + 1, i),
                                m_splines->m_projection);

                            const auto rightCosPhi = NormalizedInnerProductTwoSegments(
                                m_gridPoints(j - 1, firstRightIndex),
                                m_gridPoints(j, firstRightIndex),
                                m_gridPoints(j, firstRightIndex),
                                m_gridPoints(j + 1, firstRightIndex),
                                m_splines->m_projection);

                            const auto [secondLeftIndex, secondRightIndex] = GetNeighbours(m_gridPoints.row(j), firstRightIndex);

                            if ((secondRightIndex == firstRightIndex || leftCosPhi - rightCosPhi < -cosineTolerance) && firstLeftIndex != i)
                            {
                                // move left node
                                for (auto k = i; k <= firstRightIndex - 1; ++k)
                                {
                                    m_gridPoints(j, k) = m_gridPoints(j, firstRightIndex);
                                }
                                numChanged++;
                            }
                            else if ((firstLeftIndex == i || rightCosPhi - leftCosPhi < -cosineTolerance) && secondRightIndex != firstRightIndex)
                            {
                                // move right node
                                for (auto k = firstRightIndex; k <= secondRightIndex - 1; ++k)
                                {
                                    m_gridPoints(j, k) = m_gridPoints(j, i);
                                }
                                numChanged++;
                            }
                            else
                            {
                                // move both nodes
                                const Point middle = (m_gridPoints(j, i) + m_gridPoints(j, firstRightIndex)) * 0.5;
                                for (auto k = i; k <= firstRightIndex - 1; ++k)
                                {
                                    m_gridPoints(j, k) = middle;
                                }
                                for (auto k = firstRightIndex; k <= secondRightIndex - 1; ++k)
                                {
                                    m_gridPoints(j, k) = middle;
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
    }

    void CurvilinearGridFromSplines::Initialize()
    {
        // with only 1 spline, can't continue
        if (m_splines->GetNumSplines() <= 1)
        {
            throw std::invalid_argument("CurvilinearGridFromSplines::Initialize: Not enough splines to create a curvilinear grid.");
        }

        // compute properties
        ComputeSplineProperties(false);

        // get the properties of the center splines
        MakeAllGridLines();

        // Store original number of splines
        std::vector<Point> newCrossSpline(2);
        m_numOriginalSplines = m_splines->GetNumSplines();
        int count = 0;

        std::cout << "m_numOriginalSplines = " << m_numOriginalSplines << "  " << std::endl;

        for (UInt s = 0; s < m_numOriginalSplines; ++s)
        {
            // mirror only center splines
            if (m_type[s] != SplineTypes::central)
            {
                continue;
            }

            std::cout << "m_leftGridLineIndex[s] = " << m_leftGridLineIndex[s] << "  " << m_leftGridLineIndex[s] << "  " << m_numMSplines[s] << "  " << std::endl;

            // construct the cross splines through the edges, along m discretization
            for (auto i = m_leftGridLineIndex[s]; i < m_leftGridLineIndex[s] + m_numMSplines[s]; ++i)
            {
                const auto normal = NormalVectorOutside(m_gridLine[i], m_gridLine[i + 1], m_splines->m_projection);

                const double xMiddle = (m_gridLine[i].x + m_gridLine[i + 1].x) * 0.5;
                const double yMiddle = (m_gridLine[i].y + m_gridLine[i + 1].y) * 0.5;

                double xs1 = constants::missing::doubleValue;
                double xs2 = constants::missing::doubleValue;
                double ys1 = constants::missing::doubleValue;
                double ys2 = constants::missing::doubleValue;

                if (m_splines->m_projection == Projection::cartesian)
                {
                    xs1 = xMiddle + 2.0 * m_maximumGridHeights[s] * -normal.x;
                    xs2 = xMiddle + 2.0 * m_maximumGridHeights[s] * normal.x;
                    ys1 = yMiddle + 2.0 * m_maximumGridHeights[s] * -normal.y;
                    ys2 = yMiddle + 2.0 * m_maximumGridHeights[s] * normal.y;
                }

                if (m_splines->m_projection == Projection::spherical)
                {
                    const double factor = 1.0 / (constants::geometric::earth_radius * constants::conversion::degToRad);
                    xs1 = xMiddle + 2.0 * m_maximumGridHeights[s] * -normal.x * factor;
                    xs2 = xMiddle + 2.0 * m_maximumGridHeights[s] * normal.x * factor;
                    ys1 = yMiddle + 2.0 * m_maximumGridHeights[s] * -normal.y * factor;
                    ys2 = yMiddle + 2.0 * m_maximumGridHeights[s] * normal.y * factor;
                }

                ++count;

                if (count == 16)
                {
                    break;
                }

                std::cout << "adding spline: " << count << "  "
                          << "{ " << xMiddle << ", " << yMiddle << " }"
                          << " -- "
                          << "{ " << xs1 << ", " << ys1 << " }"
                          << " -- "
                          << "{ " << xs2 << ", " << ys2 << " }" << std::endl;

                newCrossSpline[0] = {xs1, ys1};
                newCrossSpline[1] = {xs2, ys2};
                m_splines->AddSpline(newCrossSpline);
                // flag the cross spline as artificially added
                m_type.emplace_back(SplineTypes::artificial);
            }
        }

        // Backup original spline properties
        m_leftGridLineIndexOriginal.resize(m_numOriginalSplines);
        m_rightGridLineIndexOriginal.resize(m_numOriginalSplines);
        m_mfacOriginal.resize(m_numOriginalSplines);
        m_maximumGridHeightsOriginal.resize(m_numOriginalSplines);
        m_originalTypes.resize(m_numOriginalSplines);

        for (UInt s = 0; s < m_numOriginalSplines; ++s)
        {
            m_leftGridLineIndexOriginal[s] = m_leftGridLineIndex[s];
            m_rightGridLineIndexOriginal[s] = m_rightGridLineIndex[s];
            m_mfacOriginal[s] = m_numMSplines[s];
            m_maximumGridHeightsOriginal[s] = m_maximumGridHeights[s];
            m_originalTypes[s] = m_type[s];
        }

        // compute spline properties with artificial splines added
        ComputeSplineProperties(true);

        // artificial cross spline: remove the last part of the sub-intervals (since it makes no sense,
        // as the artificial cross spline has an arbitrary, but sufficiently large, length)
        for (UInt s = 0; s < m_numOriginalSplines; ++s)
        {
            // Remove the last part of the sub-intervals
            if (m_type[s] != SplineTypes::central)
            {
                continue;
            }

            // For number of intersecting splines
            for (UInt i = 0; i < m_numCrossingSplines[s]; ++i)
            {
                const auto crossingSplineIndex = m_crossingSplinesIndices(s, i);
                if (m_type[crossingSplineIndex] == SplineTypes::artificial)
                {
                    m_numCrossSplineLeftHeights(s, i)--;
                    m_numCrossSplineRightHeights(s, i)--;
                }
            }
        }

        // Compute edge velocities
        ComputeEdgeVelocities();

        // Increase curvilinear grid
        const auto numGridLayers = m_curvilinearParameters.n_refinement + 1;
        // The layer by coordinate to grow
        // ResizeAndFill2DVector(m_gridPoints, numGridLayers + 1, m_numM + 1);
        lin_alg::ResizeAndFillMatrix(m_gridPoints, numGridLayers + 1, m_numM + 1, true);
        m_validFrontNodes.resize(m_numM, 1);

        // Copy the first m point in m_gridPoints
        for (UInt n = 0; n < m_numM; ++n)
        {
            m_gridPoints(0, n) = m_gridLine[n];

            std::cout << "gridlinex (" << n + 1 << " ) = " << m_gridLine[n].x << ";" << std::endl;
            std::cout << "gridliney (" << n + 1 << " ) = " << m_gridLine[n].y << ";" << std::endl;
            // std::cout << m_gridLine[n].x << ", " << m_gridLine[n].y << " -- ";

            if (!m_gridLine[n].IsValid())
            {
                m_validFrontNodes[n] = 0;
            }
            UInt sumLeft = 0;
            UInt sumRight = 0;
            const auto leftColumn = n == 0 ? 0 : n - 1;
            const auto rightColumn = n <= m_numM - 2 ? n : m_numM - 2;

            for (const auto& numFaces : m_numPerpendicularFacesOnSubintervalAndEdge.rowwise()) // NOT SO SURE ABOUT THIS
            {
                sumLeft += numFaces[leftColumn];
                sumRight += numFaces[rightColumn];
            }

            if (sumLeft == 0 && sumRight == 0)
            {
                m_validFrontNodes[n] = 0;
            }
        }

        std::cout << std::endl;

        // compute maximum mesh width and get dtolLR in the proper dimension
        double squaredMaximumGridWidth = 0.0;
        for (UInt i = 0; i < m_gridPoints.cols() - 1; i++)
        {
            if (!m_gridPoints(0, i).IsValid() || !m_gridPoints(0, i + 1).IsValid())
            {
                continue;
            }
            squaredMaximumGridWidth = std::max(squaredMaximumGridWidth,
                                               ComputeSquaredDistance(m_gridPoints(0, i),
                                                                      m_gridPoints(0, i + 1),
                                                                      m_splines->m_projection));
        }

        std::cout << " m_onTopOfEachOtherSquaredTolerance " << m_onTopOfEachOtherSquaredTolerance << "  " << squaredMaximumGridWidth << std::endl;

        m_onTopOfEachOtherSquaredTolerance = m_onTopOfEachOtherSquaredTolerance * squaredMaximumGridWidth;

        m_subLayerGridPoints.resize(m_numPerpendicularFacesOnSubintervalAndEdge.rows());
    }

    void CurvilinearGridFromSplines::Iterate(UInt layer)
    {
        GrowLayer(layer);

        for (UInt j = 0; j < m_subLayerGridPoints.size(); ++j)
        {
            m_subLayerGridPoints[j] = m_numPerpendicularFacesOnSubintervalAndEdge(j, 0);
        }

        auto results = ComputeGridLayerAndSubLayer(layer);
        auto subLayerRightIndex = std::get<1>(results);

        for (UInt i = 0; i < m_numM; i++)
        {
            const auto subLayerLeftIndex = subLayerRightIndex;
            const auto minRight = std::min(i, static_cast<UInt>(m_numPerpendicularFacesOnSubintervalAndEdge.cols()) - 1);

            for (UInt j = 0; j < m_subLayerGridPoints.size(); ++j)
            {
                m_subLayerGridPoints[j] = m_numPerpendicularFacesOnSubintervalAndEdge(j, minRight);
            }

            results = ComputeGridLayerAndSubLayer(layer);
            const auto gridLayer = std::get<0>(results);
            subLayerRightIndex = std::get<1>(results);

            if (subLayerRightIndex != constants::missing::uintValue && i < m_numM - 1 && gridLayer != constants::missing::uintValue)
            {
                m_edgeVelocities[i] = m_growFactorOnSubintervalAndEdge(subLayerRightIndex, i) * m_edgeVelocities[i];
            }

            if (subLayerLeftIndex == constants::missing::uintValue && subLayerRightIndex == constants::missing::uintValue)
            {
                m_validFrontNodes[i] = constants::missing::uintValue;
            }
        }

        if (m_timeStep <= 1e-8)
        {
            throw AlgorithmError("time step is smaller than 1e-8 !");
        }
    }

    std::unique_ptr<CurvilinearGrid> CurvilinearGridFromSplines::ComputeCurvilinearGridFromGridPoints()
    {
        std::vector<std::vector<Point>> gridPointsNDirection(m_gridPoints.cols(),
                                                             std::vector<Point>(m_gridPoints.rows(),
                                                                                {constants::missing::doubleValue, constants::missing::doubleValue}));
        std::vector<std::vector<Point>> curvilinearMeshPoints;
        const double squaredDistanceTolerance = 1e-12;

        // get the grid sizes in j-direction
        for (UInt i = 0; i < m_gridPoints.cols(); i++)
        {
            for (UInt j = 0; j < m_gridPoints.rows(); j++)
            {

                if (IsEqual(m_gridPoints(j, i).x, 1607.12, 0.001))
                {
                    std::cout << " m_gridPoints(j, i) " << i << "  " << j << "  " << m_gridPoints(j, i).x << ", " << m_gridPoints(j, i).y << std::endl;
                }

                gridPointsNDirection[i][j] = m_gridPoints(j, i);
            }
        }

        UInt startIndex = 0;
        UInt startGridLine = 0;
        while (startIndex < m_gridPoints.cols())
        {
            auto mIndicesThisSide = FindIndices(lin_alg::MatrixRowToSTLVector(m_gridPoints, 0),
                                                startIndex,
                                                m_numM,
                                                constants::missing::doubleValue);

            const auto& [mStartIndexThisSide, mEndIndexThisSide] = mIndicesThisSide[0];

            const auto mStartIndexOtherSide = mEndIndexThisSide + 2;
            const auto mEndIndexOtherSide = mStartIndexOtherSide + (mEndIndexThisSide - mStartIndexThisSide);

            bool isConnected = true;

            UInt minN = m_curvilinearParameters.n_refinement;
            UInt maxN = 0;
            UInt minNOther = m_curvilinearParameters.n_refinement;
            UInt maxNOther = 0;
            // check if this part is connected to another part
            for (auto i = mStartIndexThisSide; i < mEndIndexThisSide + 1; ++i)
            {
                const auto nIndicesThisSide = FindIndices(gridPointsNDirection[i],
                                                          0,
                                                          static_cast<UInt>(gridPointsNDirection[i].size()),
                                                          constants::missing::doubleValue);
                const auto& [nStartIndexThisSide, nEndIndexThisSide] = nIndicesThisSide[0];
                minN = std::min(minN, nStartIndexThisSide);
                maxN = std::max(maxN, nEndIndexThisSide);

                const UInt mOther = mEndIndexThisSide + 2 + (mEndIndexThisSide - i);

                if (mOther > m_numM - 1)
                {
                    // no more grid available
                    isConnected = false;
                }
                else
                {
                    const double squaredDistance = ComputeSquaredDistance(m_gridPoints(0, i),
                                                                          m_gridPoints(0, mOther),
                                                                          m_splines->m_projection);
                    if (squaredDistance > squaredDistanceTolerance)
                    {
                        isConnected = false;
                    }
                    else
                    {
                        const auto nIndicesOtherSide = FindIndices(gridPointsNDirection[mOther],
                                                                   0,
                                                                   static_cast<UInt>(gridPointsNDirection[mOther].size()),
                                                                   constants::missing::doubleValue);
                        const auto& [nStartIndexOtherSide, nEndIndexOtherSide] = nIndicesOtherSide[0];

                        minNOther = std::min(minNOther, nStartIndexOtherSide);
                        maxNOther = std::max(maxNOther, nEndIndexOtherSide);
                    }
                }
            }

            const auto endGridlineIndex = startGridLine + mEndIndexThisSide - mStartIndexThisSide;
            if (isConnected)
            {
                startIndex = mEndIndexOtherSide + 2;
            }
            else
            {
                maxNOther = 1;
                startIndex = mEndIndexThisSide + 2;
            }

            // increment points
            curvilinearMeshPoints.resize(endGridlineIndex + 1);
            const auto NSize = std::max(static_cast<UInt>(curvilinearMeshPoints[0].size()), maxN + maxNOther + 1);
            for (auto& element : curvilinearMeshPoints)
            {
                element.resize(NSize);
            }

            // fill first part
            UInt columnIncrement = 0;
            for (auto i = startGridLine; i < endGridlineIndex + 1; ++i)
            {
                for (UInt j = 0; j < maxN + 1; ++j)
                {
                    curvilinearMeshPoints[i][j + maxNOther] = m_gridPoints(j, mStartIndexThisSide + columnIncrement);
                }
                columnIncrement++;
            }

            columnIncrement = 0;
            for (auto i = startGridLine; i < endGridlineIndex + 1; ++i)
            {
                for (UInt j = 0; j < maxNOther + 1; ++j)
                {
                    curvilinearMeshPoints[i][maxNOther - j] = m_gridPoints(j, mEndIndexOtherSide - columnIncrement);
                }
                columnIncrement++;
            }

            startGridLine = endGridlineIndex + 2;
        }

        return std::make_unique<CurvilinearGrid>(lin_alg::STLVectorOfVectorsToMatrix(curvilinearMeshPoints), m_splines->m_projection);
    }

    std::tuple<UInt, UInt>
    CurvilinearGridFromSplines::ComputeGridLayerAndSubLayer(UInt layerIndex)
    {

        if (layerIndex == 0)
        {
            return {constants::missing::uintValue, constants::missing::uintValue};
        }

        UInt gridLayer = layerIndex - 1;
        auto sum = std::accumulate(m_subLayerGridPoints.begin(), m_subLayerGridPoints.end(), UInt(0));

        UInt subLayerIndex;
        if (layerIndex >= sum)
        {
            subLayerIndex = constants::missing::uintValue;
        }
        else
        {
            subLayerIndex = 0;
            sum = m_subLayerGridPoints[0] + 1;
            while (sum <= layerIndex && subLayerIndex < m_maxNumCenterSplineHeights)
            {
                subLayerIndex = subLayerIndex + 1;
                sum += m_subLayerGridPoints[subLayerIndex];
            }
            gridLayer = layerIndex - sum + m_subLayerGridPoints[subLayerIndex];
        }

        return {gridLayer, subLayerIndex};
    }

    void CurvilinearGridFromSplines::GetNeighbouringLayerPoints(const std::vector<Point>& activeLayerPoints,
                                                                const UInt layerPoint,
                                                                UInt& layerPointLeft,
                                                                UInt& layerPointRight) const
    {
        const double dtolLR = 1.0e-4;

        // Find left node
        // double dL1 = ComputeDistance(activeLayerPoints[layerPoint], activeLayerPoints[layerPoint + 1], m_splines->m_projection);

        layerPointLeft = layerPoint;

        while (ComputeDistance(activeLayerPoints[layerPointLeft], activeLayerPoints[layerPoint], m_splines->m_projection) <= dtolLR)
        {
            if (layerPointLeft < 1)
            {
                break;
            }

            if (!activeLayerPoints[layerPointLeft - 1].IsValid())
            {
                break;
            }

            --layerPointLeft;
        }

        layerPointRight = layerPoint;

        while (ComputeDistance(activeLayerPoints[layerPointRight], activeLayerPoints[layerPoint], m_splines->m_projection) <= dtolLR)
        {
            if (layerPointRight > 1)
            {
                break;
            }

            if (!activeLayerPoints[layerPointRight + 1].IsValid())
            {
                break;
            }

            ++layerPointRight;
        }
    }

    double CurvilinearGridFromSplines::CompCrossTime1(const Point& x1, const Point& x3, const Point& x4,
                                                      const Point& v1, const Point& v3, const Point& v4,
                                                      const double clearance) const
    {
        Point x13 = x3 - x1;
        Point x34 = x4 - x3;
        Point v13 = v3 - v1;
        Point v34 = v4 - v3;

        double a = cross(v13, v34);
        double b = cross(x13, v34) - cross(x34, v13);
        double c = cross(x13, x34);
        double e = 0.0;
        double f = 0.0;
        double g = 0.0;

        std::array<double, 5> coeffs{0.0, 0.0, a, b, c};
        std::array<double, 4> roots{0.0, 0.0, 0.0, 0.0};
        std::array<double, 4> beta;
        beta.fill(constants::missing::doubleValue);
        roots.fill(constants::missing::doubleValue);

        if (clearance > 0.0)
        {
            coeffs = {a * a, 2.0 * a * b, 2.0 * a * c + b * b, 2.0 * b * c, c * c};
            e = dot(v34, v34);
            f = 2.0 * dot(v34, v34);
            g = dot(x34, x34);

            coeffs[2] -= clearance * clearance * e;
            coeffs[3] -= clearance * clearance * f;
            coeffs[4] -= clearance * clearance * g;
        }

        SolveQuartic(coeffs, roots);

        for (UInt i = 0; i < 4; ++i)
        {

            if (roots[i] == constants::missing::doubleValue)
            {
                continue;
            }

            if (std::abs(roots[i]) < tolerance)
            {
                continue;
            }
            Point xs = x4 - x3 * (v4 - v3) * roots[i];
            double det = length(xs);

            if (std::abs(det) > tolerance)
            {
                beta[i] = dot(x3 - x1 + (v3 - v1) * roots[i], xs) / det;
            }
        }

        double time = 1.0e99;
        double DdDt;

        for (UInt i = 0; i < 4; ++i)
        {
            if (beta[i] >= 0.0 && beta[i] <= 1.0 && roots[i] >= 0.0 && roots[i] != constants::missing::doubleValue)
            {
                DdDt = (2.0 * (a * roots[i] * roots[i] + b * roots[i] + c) * (2.0 * a * roots[i] + b) - clearance * clearance * (2.0 * e * roots[i] + f)) / (2.0 * clearance * (e * roots[i] * roots[i] + f * roots[i] + g));
            }
            else
            {
                DdDt = -1.0e99;
            }

            if (DdDt < 0.0)
            {
                time = std::min(time, roots[i]);
            }
        }

        return time;
    }

    bool CurvilinearGridFromSplines::SegmentCrossesCentreSpline(const Point& x1, const Point& x2) const
    {
        const std::vector<Point>& splinePoints = m_splines->m_splineNodes[GetCentralSplineIndex()];

        for (UInt i = 0; i < splinePoints.size() - 1; ++i)
        {
            Point x3 = splinePoints[i];
            Point x4 = splinePoints[i + 1];

            if (!x3.IsValid() || !x4.IsValid())
            {
                continue;
            }

            const auto [segmentsDoCross,
                        intersection,
                        crossProduct,
                        firstRatio,
                        secondRatio] = AreSegmentsCrossing(x1, x2, x3, x4, false, m_splines->m_projection);

            if (segmentsDoCross)
            {
                return true;
            }
        }

        return false;
    }

    void CurvilinearGridFromSplines::SolveQuartic(const std::array<double, 5>& coefficients,
                                                  std::array<double, 4>& roots) const
    {
        int degree = 4;
        std::array<double, 5> coeffs(coefficients);
        std::array<double, 4> realRoots{0.0, 0.0, 0.0, 0.0};
        std::array<double, 4> imagRoots{0.0, 0.0, 0.0, 0.0};
        roots.fill(constants::missing::doubleValue);

        for (UInt i = 4; i >= 1; --i)
        {
            degree = i;

            if (coeffs[4 - degree] < tolerance)
            {
                continue;
            }

            // shuffle coefficients up

            if (degree == 2)
            {
                coeffs[0] = coeffs[2];
                coeffs[1] = coeffs[3];
                coeffs[2] = coeffs[4];
                coeffs[3] = 0.0;
                coeffs[4] = 0.0;
            }
            else if (degree == 3)
            {
                coeffs[0] = coeffs[1];
                coeffs[1] = coeffs[2];
                coeffs[2] = coeffs[3];
                coeffs[3] = coeffs[4];
                coeffs[4] = 0.0;
            }

            // for (UInt j = 0; j < degree; ++j)
            // {
            //     coeffs[j] = coeffs [];
            // }

            std::cout << " degree before " << degree << std::endl;

            RootFinder(coeffs, degree, realRoots, imagRoots);
            std::cout << " degree after " << degree << std::endl;

            std::cout << "solve quartic real: " << realRoots[0] << "  " << realRoots[1] << "  " << realRoots[2] << "  " << realRoots[3] << "  " << std::endl;
            std::cout << "solve quartic imag: " << imagRoots[0] << "  " << imagRoots[1] << "  " << imagRoots[2] << "  " << imagRoots[3] << "  " << std::endl;

            break;
        }

        if (degree <= 0)
        {
            return;
        }

        for (UInt i = 0; i < static_cast<UInt>(degree); ++i)
        {
            if (std::fabs(imagRoots[i]) < 1.0e-4)
            {
                roots[i] = realRoots[i];
            }
        }

        std::cout << "solve quartic all roots: " << roots[0] << "  " << roots[1] << "  " << roots[2] << "  " << roots[3] << "  " << std::endl;
    }

    double CurvilinearGridFromSplines::CompCrossTime2(const Point& x1, const Point& x3, const Point& x4,
                                                      const Point& v1, const Point& v3, const Point& v4,
                                                      const double clearance) const
    {
        auto [dnow, intersectionPoint, ratio] = DistanceFromLine(x1, x3, x4, m_splines->m_projection);

        double result = 1.0e99;
        double t2 = 1.0e99;

        // if ( -(x1(1)-x3(1))*(x4(2)-x3(2)) + (x1(2)-x3(2))*(x4(1)-x3(1)).lt.0d0 )

        if (-dot(x1 - x3, x4 - x3) < 0.0)
        {
            return result;
        }

        if (dnow <= clearance && clearance > 0.0)
        {
            t2 = CompCrossTime1(x1, x3, x4, v1, v3, v4, 0.0);

            if (t2 < 1.0e99)
            {
                // check if distance is increasing
                double dteps = 1e-2;
                auto [deps, intersectionPoint, ratio] = DistanceFromLine(x1 + v1 * dteps, x3 + v3 * dteps, x4 + v4 * dteps, m_splines->m_projection);

                // dlinedis(x1.x + v1.x * dteps, x1.y + v1.y * dteps, x3.x + v3.x * dteps, x3.y + v3.y * dteps, x4.x + v4.x * dteps, x4.y + v4.y * dteps, ja, dnow, xc, yc, jsferic, jasfer3D, dmiss);
                double DdDt = (deps - dnow) / dteps;

                if (DdDt < -1e-4)
                {
                    // t2 = comp_cross_time_1(x1,x3,x4,v1,v3,v4,0d0);
                    t2 = 0.0;
                }
                else
                {
                    t2 = CompCrossTime1(x1, x3, x4, v1, v3, v4, 0.0);
                }
            }

            std::cout << " time t2 " << t2 << std::endl;

            return t2;
        }

        double t1 = CompCrossTime1(x1, x3, x4, v1, v3, v4, clearance);

        std::cout << " time t1 before " << t1 << std::endl;

        if (t1 == constants::missing::doubleValue || t1 <= 0.0)
        {
            t1 = 1.0e99;
        }

        std::cout << " time t1 after " << t1 << std::endl;

        double a = dot(v1 - v3, v1 - v3);
        double b = 2.0 * dot(v1 - v3, x1 - x3);
        double c = dot(x1 - x3, x1 - x3);

        std::array<double, 5> coeffs = {0.0, 0.0, a, b, c - clearance * clearance};
        std::array<double, 4> roots{-2000.0, -2000.0, -2000.0, -2000.0};

        SolveQuartic(coeffs, roots);

        std::cout << "coeffs: " << coeffs[0] << "  " << coeffs[1] << "  " << coeffs[2] << "  " << coeffs[3] << "  " << coeffs[4] << "  " << std::endl;
        std::cout << "roots:  " << roots[0] << "  " << roots[1] << "  " << roots[2] << "  " << roots[3] << "  " << std::endl;

        for (int i = 0; i < 4; ++i)
        {
            if (roots[i] == constants::missing::doubleValue || roots[i] <= 0.0 || roots[i] > t1)
            {
                continue;
            }

            if (dot(x1 - x3 + (v1 - v3) * roots[i], x4 - x3 + (v4 - v3) * roots[i]) > 0.0)
            {
                continue;
            }

            // check if distance is decreasing
            double DdDt = 1.0e99;

            if (clearance > 0.0 && roots[i] > 0.0)
            {
                // check if the new connecting line does not cross the center spline gridline
                Point xdum1 = x1 + v1 * roots[i]; //{x1.x + v1.x * roots[i], x1.y + v1.y * roots[i]};
                Point xdum2 = x3 + v3 * roots[i]; //{x3.x + v3.x * roots[i], x3.y + v3.y * roots[i]};

                if (!SegmentCrossesCentreSpline(xdum1, xdum2))
                {
                    DdDt = (2.0 * a * roots[i] + b) / (2.0 * clearance);
                }
            }

            if (roots[i] != constants::missing::doubleValue && roots[i] > 0.0 && DdDt < 0.0)
            {
                t1 = std::min(t1, roots[i]);
            }
        }

        //--------------------------------

        a = dot(v1 - v4, v1 - v4);
        b = 2.0 * dot(v1 - v4, x1 - x4);
        c = dot(x1 - x4, x1 - x4);

        coeffs = std::array<double, 5>{0.0, 0.0, a, b, c - clearance * clearance};

        SolveQuartic(coeffs, roots);

        for (int i = 0; i < 4; ++i)
        {
            if (roots[i] == constants::missing::doubleValue || roots[i] <= 0.0 || roots[i] > t1)
            {
                continue;
            }

            if (dot(x1 - x4 + (v1 - v4) * roots[i], x4 - x3 + (v4 - v3) * roots[i]) > 0.0)
            {
                continue;
            }

            // check if distance is decreasing
            double DdDt = 1.0e99;

            if (clearance > 0.0 && roots[i] > 0.0)
            {
                // check if the new connecting line does not cross the center spline gridline
                Point xdum1 = x1 + v1 * roots[i]; //{x1.x + v1.x * roots[i], x1.y + v1.y * roots[i]};
                Point xdum2 = x3 + v3 * roots[i]; //{x3.x + v3.x * roots[i], x3.y + v3.y * roots[i]};

                if (!SegmentCrossesCentreSpline(xdum1, xdum2))
                {
                    DdDt = (2.0 * a * roots[i] + b) / (2.0 * clearance);
                }
            }

            if (roots[i] != constants::missing::doubleValue && roots[i] > 0.0 && DdDt < 0.0)
            {
                t1 = std::min(t1, roots[i]);
            }
        }

        //--------------------------------

        return std::min(t1, t2);
    }

    double CurvilinearGridFromSplines::ComputeMaximumTimeStep(const UInt layerIndex,
                                                              const std::vector<Point>& activeLayerPoints,
                                                              const std::vector<Point>& velocityVectorAtGridPoints,
                                                              const std::vector<Point>& frontVelocities,
                                                              const double timeStep) const
    {
        double maximumTimeStep = std::numeric_limits<double>::max();
        // Should be in header
        const double dtolLR = 1.0e-4;

        UInt nummax = 2 * static_cast<UInt>(activeLayerPoints.size());

        std::vector<double> tmax(activeLayerPoints.size());
        std::ranges::fill(tmax, timeStep);

        Point x1;
        Point x2;
        Point x3;
        Point x4;

        Point v1;
        Point v2;
        Point v3;
        Point v4;

        UInt iMax = 0;
        UInt iMin = constants::missing::uintValue;

        // Only need frontGridPoints
        auto [gridPointsIndices, frontGridPointsEigen, numFrontPoints] = FindFront();

        std::span frontGridPoints(frontGridPointsEigen.data(), frontGridPointsEigen.size());

        for (UInt i = 0; i < activeLayerPoints.size() - 1; ++i)
        {
            if (!activeLayerPoints[i].IsValid() || !activeLayerPoints[i + 1].IsValid())
            {
                continue;
            }

            x1 = activeLayerPoints[i];
            x2 = activeLayerPoints[i + 1];

            v1 = velocityVectorAtGridPoints[i];
            v2 = velocityVectorAtGridPoints[i];

            double dL1 = ComputeDistance(x1, activeLayerPoints[i + 1], m_splines->m_projection);
            UInt iL;
            UInt iLL;
            UInt iR;
            UInt iRR;
            UInt dummy;

            // TODO combine calculation of the iL and iR, removing the need for the dummy and two calls
            GetNeighbouringLayerPoints(activeLayerPoints, i, iL, dummy);
            GetNeighbouringLayerPoints(activeLayerPoints, i + 1, dummy, iR);

            GetNeighbouringLayerPoints(activeLayerPoints, iL, iLL, dummy);
            GetNeighbouringLayerPoints(activeLayerPoints, iR, dummy, iRR);

            dummy = iL;

            for (UInt j = 0; j < nummax; ++j)
            {
                UInt i1;
                GetNeighbouringLayerPoints(activeLayerPoints, dummy, iMin, i1);

                if (iMin == dummy)
                {
                    break;
                }

                dummy = iMin;
            }

            dummy = iR;

            for (UInt j = 0; j < nummax; ++j)
            {
                UInt i1;
                GetNeighbouringLayerPoints(activeLayerPoints, dummy, i1, iMax);

                if (iMax == dummy)
                {
                    break;
                }

                dummy = iMax;
            }

            for (UInt j = 0; j < frontGridPoints.size() - 1; ++j)
            {
                if (!frontGridPoints[j].IsValid() || !frontGridPoints[j + 1].IsValid())
                {
                    continue;
                }

                x3 = frontGridPoints[j];
                x4 = frontGridPoints[j + 1];

                v3 = frontVelocities[j];
                v4 = frontVelocities[j + 1];

                double dL2 = ComputeDistance(x3, frontGridPoints[j + 1], m_splines->m_projection);

                if (ComputeDistance(x1, x3, m_splines->m_projection) < dtolLR ||
                    ComputeDistance(activeLayerPoints[i + 1], frontGridPoints[j + 1], m_splines->m_projection) < dtolLR)
                {
                    continue;
                }

                if (ComputeDistance(activeLayerPoints[i + 1], x3, m_splines->m_projection) < dtolLR ||
                    ComputeDistance(x1, frontGridPoints[j + 1], m_splines->m_projection) < dtolLR)
                {
                    continue;
                }

                double d1 = ComputeDistance(x1, x3, m_splines->m_projection);
                double d2 = ComputeDistance(x2, x3, m_splines->m_projection);
                double d3 = ComputeDistance(x1, x4, m_splines->m_projection);
                double d4 = ComputeDistance(x2, x4, m_splines->m_projection);

                if (d1 < tolerance || d2 < tolerance || d3 < tolerance || d4 < tolerance)
                {
                    continue;
                }

                double clearance = 0.0;

                UInt i1 = gridPointsIndices(j, 0);
                UInt i2 = gridPointsIndices(j + 1, 0);
                UInt j1 = gridPointsIndices(j, 1);
                UInt j2 = gridPointsIndices(j + 1, 1);

                if ((iMin <= i1 && i1 <= iMax) || (iMin <= i2 && i2 <= iMax))
                {
                    clearance = 0.0;
                }

                if (iRR > iLL)
                {
                    if (((i1 > iLL && i1 < iRR) || (i2 > iLL && i2 < iRR)) && j1 >= (layerIndex - 1) && j2 >= (layerIndex - 1))
                    {
                        continue;
                    }
                }
                else
                {
                    if ((!(i1 >= iRR && i1 <= iLL) || !(i2 >= iRR && i2 <= iLL)) && j1 >= (layerIndex - 1) && j2 >= (layerIndex - 1))
                    {
                        continue;
                    }
                }

                double dmin = std::min({d1, d2, d3, d4});

                // get a lower bound for the cross time
                double hlow2 = 0.25 * std::max(dmin * dmin - 0.5 * std::pow(std::max(dL1, dL2), 2), 0.0);

                // check if the lower bounds is larger than the minimum found so far
                double vv1 = std::sqrt(length(v2 - v1));
                double vv2 = std::sqrt(length(v3 - v2));
                double vv3 = std::sqrt(length(v4 - v1));
                double vv4 = std::sqrt(length(v4 - v2));

                double maxvv = std::max(std::max(vv1, vv2), std::max(vv3, vv4));

                if (std::sqrt(hlow2) - clearance > maxvv * std::min(tmax[i], tmax[i + 1]))
                {
                    continue;
                }

                double t1 = CompCrossTime2(x1, x3, x4, v1, v3, v4, clearance);
                double t2 = CompCrossTime2(x2, x3, x4, v2, v3, v4, clearance);
                double t3 = CompCrossTime2(x3, x1, x2, v3, v1, v2, clearance);
                double t4 = CompCrossTime2(x4, x1, x2, v4, v1, v1, clearance);

                std::cout << "crooss time : " << t1 << "  " << t2 << "  " << t3 << "  " << t4 << "  " << std::endl;

                // Get the name right
                double tmax1234 = std::min(std::min(t1, t2), std::min(t3, t4));

                if (t1 == tmax1234)
                {
                    tmax[i] = std::min(tmax[i], tmax1234);
                }
                else if (t2 == tmax1234)
                {
                    tmax[i + 1] = std::min(tmax[i + 1], tmax1234);
                }
                else if (t3 == tmax1234 || t4 == tmax1234)
                {
                    tmax[i] = std::min(tmax[i], tmax1234);
                    tmax[i + 1] = std::min(tmax[i + 1], tmax1234);
                }

                if (tmax1234 == 0.0)
                {
                    break;
                }
            }
        }

        return maximumTimeStep;
    }

    UInt CurvilinearGridFromSplines::GetCentralSplineIndex() const
    {

        for (UInt i = 0; i < m_type.size(); ++i)
        {

            if (m_type[i] == SplineTypes::central)
            {
                return i;
            }
        }

        return constants::missing::uintValue;
    }

    void CurvilinearGridFromSplines::GrowLayer(UInt layerIndex)
    {
        auto velocityVectorAtGridPoints = ComputeVelocitiesAtGridPoints(layerIndex - 1);

        std::vector<Point> activeLayerPoints(lin_alg::MatrixRowToSTLVector(m_gridPoints, layerIndex - 1));
        for (UInt m = 0; m < velocityVectorAtGridPoints.size(); ++m)
        {
            if (!velocityVectorAtGridPoints[m].IsValid())
            {
                m_gridPoints(layerIndex - 1, m) = {constants::missing::doubleValue, constants::missing::doubleValue};
                activeLayerPoints[m] = {constants::missing::doubleValue, constants::missing::doubleValue};
            }
        }

        auto frontVelocities = CopyVelocitiesToFront(layerIndex - 1, velocityVectorAtGridPoints);

        double totalTimeStep = 0.0;
        double localTimeStep = 0.0;
        double otherTimeStep = std::numeric_limits<double>::max();
        const auto numGridPoints = static_cast<UInt>(m_gridPoints.size());
        std::vector<UInt> newValidFrontNodes(numGridPoints);

        while (totalTimeStep < m_timeStep)
        {
            std::cout << "time step: " << localTimeStep << "  " << m_timeStep << "  " << totalTimeStep << std::endl;

            // Copy old front velocities
            newValidFrontNodes = m_validFrontNodes;

            for (UInt i = 0; i < m_validFrontNodes.size(); ++i)
            {
                if (m_validFrontNodes[i] == constants::missing::uintValue)
                {
                    activeLayerPoints[i] = {constants::missing::doubleValue, constants::missing::doubleValue};
                }
            }

            const auto maximumGridLayerGrowTime = ComputeMaximumEdgeGrowTime(activeLayerPoints, velocityVectorAtGridPoints);
            localTimeStep = std::min(m_timeStep - totalTimeStep, *std::min_element(maximumGridLayerGrowTime.begin(), maximumGridLayerGrowTime.end()));

            if (m_splinesToCurvilinearParameters.check_front_collisions)
            {
                // TODO: implement front collisions
                // otherTimeStep = 0.0;
                otherTimeStep = ComputeMaximumTimeStep(layerIndex, activeLayerPoints, velocityVectorAtGridPoints, frontVelocities, localTimeStep + 1.0);
                std::cout << " otherTimeStep " << otherTimeStep << std::endl;
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

            for (UInt i = 0; i < newValidFrontNodes.size() - 2; ++i)
            {
                if (newValidFrontNodes[i + 1] == 1 && newValidFrontNodes[i] == 0 && newValidFrontNodes[i + 2] == 0)
                {
                    newValidFrontNodes[i + 1] = 0;
                }
            }

            m_validFrontNodes = newValidFrontNodes;

            for (UInt i = 0; i < velocityVectorAtGridPoints.size(); ++i)
            {

                if (m_validFrontNodes[i] == 1 && velocityVectorAtGridPoints[i].IsValid())
                {
                    activeLayerPoints[i].x = activeLayerPoints[i].x + localTimeStep * velocityVectorAtGridPoints[i].x;
                    activeLayerPoints[i].y = activeLayerPoints[i].y + localTimeStep * velocityVectorAtGridPoints[i].y;
                }
                else
                {
                    activeLayerPoints[i].x = constants::missing::doubleValue;
                    activeLayerPoints[i].y = constants::missing::doubleValue;
                }
            }

            // update the grid points
            m_gridPoints.row(layerIndex) = lin_alg::RowVector<Point>::Map(activeLayerPoints.data(),
                                                                          1,
                                                                          activeLayerPoints.size());

            // update the time step
            totalTimeStep += localTimeStep;

            if (totalTimeStep < m_timeStep)
            {
                velocityVectorAtGridPoints = ComputeVelocitiesAtGridPoints(layerIndex);

                for (UInt i = 0; i < m_numM; ++i)
                {
                    // Disable points that have no valid normal vector
                    // Remove stationary points
                    if (!frontVelocities[i].IsValid() || m_validFrontNodes[i] == 0)
                    {
                        activeLayerPoints[i] = {constants::missing::doubleValue, constants::missing::doubleValue};
                    }
                }
                frontVelocities = CopyVelocitiesToFront(layerIndex - 1, velocityVectorAtGridPoints);
            }
        }

        auto [gridPointsIndices, frontGridPoints, numFrontPoints] = FindFront();
        if (layerIndex >= 2)
        {
            for (UInt i = 1; i < m_numM - 1; ++i)
            {

                if (!activeLayerPoints[i].IsValid())
                {
                    continue;
                }

                const auto cosphi = NormalizedInnerProductTwoSegments(m_gridPoints(layerIndex - 2, i),
                                                                      m_gridPoints(layerIndex - 1, i),
                                                                      m_gridPoints(layerIndex - 1, i),
                                                                      activeLayerPoints[i],
                                                                      m_splines->m_projection);

                if (cosphi < -0.5)
                {
                    const auto [currentLeftIndex, currentRightIndex] = GetNeighbours(frontGridPoints, i);

                    for (auto j = currentLeftIndex + 1; j < currentRightIndex; ++j)
                    {
                        newValidFrontNodes[j] = 0;
                        m_gridPoints(layerIndex - 1, j) = {constants::missing::doubleValue, constants::missing::doubleValue};
                    }
                }
            }
        }

        m_validFrontNodes = newValidFrontNodes;
    }

    std::vector<double> CurvilinearGridFromSplines::ComputeMaximumEdgeGrowTime(const std::vector<Point>& coordinates,
                                                                               const std::vector<Point>& velocities) const
    {
        std::vector<double> maximumGridLayerGrowTime(m_validFrontNodes.size(), std::numeric_limits<double>::max());
        std::vector<double> edgeWidth(coordinates.size() - 1);
        std::vector<double> edgeIncrement(coordinates.size() - 1);
        const double minEdgeWidth = 1e-8;
        const double dt = 1.0;
        for (UInt i = 0; i < coordinates.size() - 1; ++i)
        {
            if (!coordinates[i].IsValid() || !coordinates[i + 1].IsValid())
            {
                continue;
            }

            edgeWidth[i] = ComputeDistance(coordinates[i], coordinates[i + 1], m_splines->m_projection);

            if (edgeWidth[i] < minEdgeWidth)
            {
                continue;
            }

            Point firstPointIncremented(coordinates[i] + velocities[i] * dt);
            Point secondPointIncremented(coordinates[i + 1] + velocities[i + 1] * dt);

            edgeIncrement[i] = InnerProductTwoSegments(coordinates[i], coordinates[i + 1], firstPointIncremented, secondPointIncremented, m_splines->m_projection) / edgeWidth[i] - edgeWidth[i];
            edgeIncrement[i] = edgeIncrement[i] / dt;
        }

        for (UInt i = 0; i < coordinates.size() - 1; ++i)
        {
            if (edgeIncrement[i] < 0.0)
            {
                maximumGridLayerGrowTime[i] = -edgeWidth[i] / edgeIncrement[i];
            }
        }

        return maximumGridLayerGrowTime;
    }

    std::vector<Point> CurvilinearGridFromSplines::CopyVelocitiesToFront(UInt layerIndex,
                                                                         const std::vector<Point>& previousFrontVelocities)
    {
        const auto numGridPoints = m_gridPoints.size();
        std::vector<Point> velocities(numGridPoints, {0.0, 0.0});

        UInt p = 0;
        auto [gridPointsIndices, frontGridPoints, numFrontPoints] = FindFront();
        while (p < numFrontPoints)
        {
            if (gridPointsIndices(p, 1) == layerIndex && m_validFrontNodes[gridPointsIndices(p, 0)] == 1)
            {
                velocities[p] = previousFrontVelocities[gridPointsIndices(p, 0)];
                if (!velocities[p].IsValid())
                {
                    velocities[p] = {0.0, 0.0};
                }

                // Check for corner nodes
                const auto previous = p == 0 ? 0 : p - 1;
                const auto previousIndices = gridPointsIndices.row(previous);
                const auto next = std::min(p + 1, numFrontPoints);
                const auto nextIndices = gridPointsIndices.row(next);

                // Check corner nodes
                bool ll = previousIndices[0] == gridPointsIndices(p, 0) - 1 &&
                          previousIndices[1] == gridPointsIndices(p, 1) &&
                          m_validFrontNodes[previousIndices[0]] == constants::missing::uintValue;

                bool lr = nextIndices[0] == gridPointsIndices(p, 0) + 1 &&
                          nextIndices[1] == gridPointsIndices(p, 1) &&
                          m_validFrontNodes[nextIndices[0]] == constants::missing::uintValue;

                ll = ll || (previousIndices[0] == gridPointsIndices(p, 0) && previousIndices[1] < gridPointsIndices(p, 1));
                lr = lr || (nextIndices[0] == gridPointsIndices(p, 0) && nextIndices[1] < gridPointsIndices(p, 1));
                if (ll || lr)
                {
                    if (numFrontPoints + 1 > frontGridPoints.size())
                    {
                        continue;
                    }
                    for (auto i = numFrontPoints; i >= p; --i)
                    {
                        frontGridPoints[i + 1] = frontGridPoints[i];
                        velocities[i + 1] = velocities[i];
                        gridPointsIndices.row(i + 1) = gridPointsIndices.row(i);
                    }
                    numFrontPoints++;

                    if (ll)
                    {
                        velocities[p] = {0.0, 0.0};
                    }
                    else
                    {
                        velocities[p + 1] = {0.0, 0.0};
                    }
                    p = p + 1;
                }
            }
            p = p + 1;
        }

        return velocities;
    }

    std::tuple<lin_alg::Matrix<UInt>,
               lin_alg::RowVector<Point>,
               UInt>
    CurvilinearGridFromSplines::FindFront() const
    {
        UInt const numGridPoints = static_cast<UInt>(m_gridPoints.size());
        lin_alg::Matrix<UInt> gridPointsIndices(numGridPoints, 2);
        gridPointsIndices.fill(constants::missing::uintValue);
        lin_alg::RowVector<Point> frontGridPoints(numGridPoints);
        frontGridPoints.fill(Point{0.0, 0.0});
        UInt numFrontPoints;

        std::vector<int> frontPosition(m_gridPoints.cols() - 2,
                                       static_cast<int>(m_gridPoints.rows()));
        for (UInt m = 0; m < frontPosition.size(); ++m)
        {
            for (UInt n = 0; n < m_gridPoints.rows(); ++n)
            {
                if (!m_gridPoints(n, m).IsValid() || !m_gridPoints(n, m + 1).IsValid())
                {
                    frontPosition[m] = n - 1;
                    break;
                }
            }
        }

        numFrontPoints = 0;
        // check for circular connectivity
        int previousFrontPosition = 0;
        const auto [leftNode, rightNode] = GetNeighbours(m_gridPoints.row(0), 0);
        if (leftNode == 0)
        {
            frontGridPoints[0] = m_gridPoints(0, 0);
            // store front index
            gridPointsIndices(numFrontPoints, 0) = 0;
            gridPointsIndices(numFrontPoints, 1) = 0;
            numFrontPoints++;
        }
        else
        {
            previousFrontPosition = frontPosition[leftNode];
            frontGridPoints[numFrontPoints] = m_gridPoints(0, frontPosition[0]);
            gridPointsIndices(numFrontPoints, 0) = frontPosition[0];
            gridPointsIndices(numFrontPoints, 1) = 0;
            numFrontPoints++;
        }

        for (UInt m = 0; m < m_gridPoints.cols() - 2; ++m)
        {
            const auto currentFrontPosition = frontPosition[m];
            if (currentFrontPosition >= 0)
            {
                if (previousFrontPosition == -1)
                {
                    frontGridPoints[numFrontPoints] = m_gridPoints(0, m);
                    gridPointsIndices(numFrontPoints, 0) = m;
                    gridPointsIndices(numFrontPoints, 1) = 0;
                    numFrontPoints++;
                }
                for (auto i = previousFrontPosition + 1; i <= currentFrontPosition; ++i)
                {
                    frontGridPoints[numFrontPoints] = m_gridPoints(i, m);
                    gridPointsIndices(numFrontPoints, 0) = m;
                    gridPointsIndices(numFrontPoints, 1) = i;
                    numFrontPoints++;
                }
                for (auto i = previousFrontPosition; i > currentFrontPosition; --i)
                {
                    frontGridPoints[numFrontPoints] = m_gridPoints(i, m);
                    gridPointsIndices(numFrontPoints, 0) = m;
                    gridPointsIndices(numFrontPoints, 1) = i;
                    numFrontPoints++;
                }

                frontGridPoints[numFrontPoints] = m_gridPoints(currentFrontPosition, m + 1);
                gridPointsIndices(numFrontPoints, 0) = m + 1;
                gridPointsIndices(numFrontPoints, 1) = currentFrontPosition;
                numFrontPoints++;
            }
            else if (previousFrontPosition >= 0)
            {
                for (auto i = previousFrontPosition - 1; i >= 0; --i)
                {
                    frontGridPoints[numFrontPoints] = m_gridPoints(i, m);
                    gridPointsIndices(numFrontPoints, 0) = m;
                    gridPointsIndices(numFrontPoints, 1) = i;
                    numFrontPoints++;
                }

                frontGridPoints[numFrontPoints] = {constants::missing::doubleValue, constants::missing::doubleValue};
                gridPointsIndices(numFrontPoints, 0) = m;
                gridPointsIndices(numFrontPoints, 1) = constants::missing::uintValue;
                numFrontPoints++;
            }

            previousFrontPosition = currentFrontPosition;
        }

        // add last j-edge, check for circular connectivity
        const auto lastPoint = static_cast<UInt>(m_gridPoints.cols()) - 2;
        const auto [currentLeftIndex, currentRightIndex] = GetNeighbours(m_gridPoints.row(0), lastPoint);
        if (currentRightIndex == m_gridPoints.cols() - 2)
        {
            for (auto i = previousFrontPosition; i >= 0; --i)
            {
                frontGridPoints[numFrontPoints] = m_gridPoints(i, lastPoint);
                gridPointsIndices(numFrontPoints, 0) = lastPoint;
                gridPointsIndices(numFrontPoints, 1) = i;
                numFrontPoints++;
            }
        }

        return {gridPointsIndices, frontGridPoints, numFrontPoints};
    }

    std::vector<Point>
    CurvilinearGridFromSplines::ComputeVelocitiesAtGridPoints(UInt layerIndex)
    {
        std::vector<Point> velocityVector(m_numM);
        std::fill(velocityVector.begin(), velocityVector.end(), Point());
        Point normalVectorLeft;
        Point normalVectorRight;
        const double cosTolerance = 1e-8;
        const double eps = 1e-10;
        for (UInt m = 0; m < velocityVector.size(); ++m)
        {
            if (!m_gridPoints(layerIndex, m).IsValid())
            {
                continue;
            }

            const auto [currentLeftIndex, currentRightIndex] = GetNeighbours(m_gridPoints.row(layerIndex), m);
            const auto squaredLeftRightDistance = ComputeSquaredDistance(m_gridPoints(layerIndex, currentLeftIndex),
                                                                         m_gridPoints(layerIndex, currentRightIndex),
                                                                         m_splines->m_projection);

            if (squaredLeftRightDistance <= m_onTopOfEachOtherSquaredTolerance)
            {
                continue;
            }

            const auto squaredLeftDistance = ComputeSquaredDistance(m_gridPoints(layerIndex, currentLeftIndex),
                                                                    m_gridPoints(layerIndex, m),
                                                                    m_splines->m_projection);
            const auto squaredRightDistance = ComputeSquaredDistance(m_gridPoints(layerIndex, currentRightIndex),
                                                                     m_gridPoints(layerIndex, m),
                                                                     m_splines->m_projection);

            if (squaredLeftDistance <= m_onTopOfEachOtherSquaredTolerance || squaredRightDistance <= m_onTopOfEachOtherSquaredTolerance)
            {
                std::cout << "here 1" << std::endl;
                normalVectorLeft = NormalVectorOutside(m_gridPoints(layerIndex, currentRightIndex),
                                                       m_gridPoints(layerIndex, currentLeftIndex),
                                                       m_splines->m_projection);
                if (m_splines->m_projection == Projection::spherical)
                {
                    normalVectorLeft.x = normalVectorLeft.x * std::cos(constants::conversion::degToRad * 0.5 *
                                                                       (m_gridPoints(layerIndex, currentLeftIndex).y +
                                                                        m_gridPoints(layerIndex, currentRightIndex).y));
                }
                normalVectorRight = normalVectorLeft;
            }
            else
            {
                std::cout << "here 2" << std::endl;
                normalVectorLeft = NormalVectorOutside(m_gridPoints(layerIndex, m),
                                                       m_gridPoints(layerIndex, currentLeftIndex),
                                                       m_splines->m_projection);
                normalVectorRight = NormalVectorOutside(m_gridPoints(layerIndex, currentRightIndex),
                                                        m_gridPoints(layerIndex, m),
                                                        m_splines->m_projection);

                if (m_splines->m_projection == Projection::spherical)
                {
                    normalVectorLeft.x = normalVectorLeft.x * std::cos(constants::conversion::degToRad * 0.5 *
                                                                       (m_gridPoints(layerIndex, currentLeftIndex).y +
                                                                        m_gridPoints(layerIndex, m).y));
                    normalVectorRight.x = normalVectorRight.x * std::cos(constants::conversion::degToRad * 0.5 *
                                                                         (m_gridPoints(layerIndex, currentRightIndex).y +
                                                                          m_gridPoints(layerIndex, m).y));
                }
            }

            if (currentLeftIndex == velocityVector.size() - 1)
            {
                continue;
            }

            const auto cosphi = DotProduct(normalVectorLeft.x, normalVectorRight.x, normalVectorLeft.y, normalVectorRight.y);
            if (cosphi < -1.0 + cosTolerance)
            {
                continue;
            }
            const auto leftVelocity = normalVectorLeft * m_edgeVelocities[currentLeftIndex];
            const auto rightVelocity = normalVectorRight * m_edgeVelocities[currentRightIndex - 1];
            double rightLeftVelocityRatio = m_edgeVelocities[currentRightIndex - 1] / m_edgeVelocities[currentLeftIndex];

            int value = 0;

            // TODO Probably should be : std::abs(cosphi) <= cosTolerance
            if ((rightLeftVelocityRatio - cosphi > eps && 1.0 / rightLeftVelocityRatio - cosphi > eps)) // || std::abs(cosphi) <= cosTolerance)
            {
                value = 1;
                velocityVector[m] = (leftVelocity * (1.0 - rightLeftVelocityRatio * cosphi) +
                                     rightVelocity * (1.0 - 1.0 / rightLeftVelocityRatio * cosphi)) /
                                    (1.0 - cosphi * cosphi);
            }
            else if (cosphi - rightLeftVelocityRatio > eps)
            {
                value = 2;
                velocityVector[m] = leftVelocity * rightLeftVelocityRatio / cosphi;
            }
            else
            {
                value = 3;
                velocityVector[m] = rightVelocity * 1.0 / (rightLeftVelocityRatio * cosphi);
            }

            // if (m == 28 || m == 29 || m == 30)
            {
                std::cout << " m = " << m << "  " << value
                          << " velocity: "
                          << velocityVector[m].x << "   " << velocityVector[m].y << " ---  "
                          << m_edgeVelocities[currentRightIndex - 1] << "  " << m_edgeVelocities[currentLeftIndex] << " --- "
                          << cosphi << "  " << (1.0 - cosphi * cosphi) << "  " << rightLeftVelocityRatio << " --- "
                          << leftVelocity.x << ", " << leftVelocity.y << " * " << (1.0 - rightLeftVelocityRatio * cosphi) << " --- "
                          << rightVelocity.x << ", " << rightVelocity.y << " * " << (1.0 - 1.0 / rightLeftVelocityRatio * cosphi) << " --- "
                          << rightVelocity.x * (1.0 - 1.0 / rightLeftVelocityRatio * cosphi) / (1.0 - cosphi * cosphi) << " -- "
                          << m_gridPoints(layerIndex, m).x << ", " << m_gridPoints(layerIndex, m).y << " -- "
                          << std::endl;
            }

            if (m_splines->m_projection == Projection::spherical)
            {
                velocityVector[m].x = velocityVector[m].x * constants::geometric::inverse_earth_radius * constants::conversion::radToDeg /
                                      std::cos(constants::conversion::degToRad * m_gridPoints(layerIndex, m).y);
                velocityVector[m].y = velocityVector[m].y * constants::geometric::inverse_earth_radius * constants::conversion::radToDeg;
            }
        }

        return velocityVector;
    }

    std::pair<UInt, UInt> CurvilinearGridFromSplines::GetNeighbours(lin_alg::RowVector<Point> const& gridPoints,
                                                                    UInt index) const
    {

        if (gridPoints.size() == 0)
        {
            return {constants::missing::uintValue, constants::missing::uintValue};
        }

        bool circularConnection = false;
        auto localLeftIndex = static_cast<int>(index);
        auto localRightIndex = static_cast<int>(index);

        // left
        while (ComputeSquaredDistance(gridPoints[localLeftIndex], gridPoints[index], m_splines->m_projection) < m_onTopOfEachOtherSquaredTolerance)
        {
            if (!circularConnection)
            {
                if (localLeftIndex == 0)
                {
                    break;
                }
            }
            else if (localLeftIndex == 0)
            {
                localLeftIndex = static_cast<int>(gridPoints.size());
                circularConnection = false;
            }

            if (localLeftIndex == 0 || !gridPoints[localLeftIndex - 1].IsValid())
            {
                break;
            }
            localLeftIndex--;
        }

        // right
        while (ComputeSquaredDistance(gridPoints[localRightIndex], gridPoints[index], m_splines->m_projection) < m_onTopOfEachOtherSquaredTolerance)
        {
            if (!circularConnection)
            {
                if (localRightIndex == static_cast<int>(gridPoints.size()))
                {
                    break;
                }
            }
            else if (localRightIndex == static_cast<int>(gridPoints.size()))
            {
                localRightIndex = -1;
                circularConnection = false;
            }

            if (localRightIndex == static_cast<int>(gridPoints.size()) || !gridPoints[localRightIndex + 1].IsValid())
            {
                break;
            }
            localRightIndex++;
        }

        return {localLeftIndex, localRightIndex};
    }

    void CurvilinearGridFromSplines::ComputeEdgeVelocities()
    {
        m_edgeVelocities.resize(m_numM - 1, constants::missing::doubleValue);

        lin_alg::ResizeAndFillMatrix(m_growFactorOnSubintervalAndEdge,
                                     m_maxNumCenterSplineHeights,
                                     m_numM - 1,
                                     false,
                                     1.0);

        lin_alg::ResizeAndFillMatrix(m_numPerpendicularFacesOnSubintervalAndEdge,
                                     m_maxNumCenterSplineHeights,
                                     m_numM - 1,
                                     false,
                                     static_cast<UInt>(0));

        ComputeGridHeights();

        std::fill(m_numPerpendicularFacesOnSubintervalAndEdge.row(0).begin(),
                  m_numPerpendicularFacesOnSubintervalAndEdge.row(0).end(),
                  1);

        for (UInt s = 0; s < m_splines->GetNumSplines(); s++)
        {

            if (m_type[s] != SplineTypes::central)
            {
                continue;
            }

            // Get true crossing splines heights
            auto numLeftHeights = m_maxNumCenterSplineHeights;
            auto numRightHeights = m_maxNumCenterSplineHeights;
            UInt numTrueCrossings = 0;
            for (UInt i = 0; i < m_numCrossingSplines[s]; ++i)
            {
                if (m_type[m_crossingSplinesIndices(s, i)] != SplineTypes::crossing)
                {
                    // true crossing splines only
                    continue;
                }
                numTrueCrossings++;
                numLeftHeights = std::min(numLeftHeights, m_numCrossSplineLeftHeights(s, i));
                numRightHeights = std::min(numRightHeights, m_numCrossSplineRightHeights(s, i));
            }

            // no true cross splines: exponentially growing grid only
            if (numTrueCrossings == 0)
            {
                numLeftHeights = 0;
                numRightHeights = 0;
            }

            const auto startGridLineLeft = m_leftGridLineIndex[s];
            const auto endGridLineLeft = startGridLineLeft + m_numMSplines[s];
            const auto startGridLineRight = m_rightGridLineIndex[s];
            const auto endGridLineRight = startGridLineRight + m_numMSplines[s];

            double hh0LeftMaxRatio;
            double hh0RightMaxRatio;
            const UInt numIterations = 2;

            double maxHeight = std::numeric_limits<double>::lowest();
            for (const auto& e : m_gridHeights.row(0))
            {
                if (!IsEqual(e, constants::missing::doubleValue) && e > maxHeight)
                {
                    maxHeight = e;
                }
            }
            const double firstHeight = std::min(maxHeight, m_splinesToCurvilinearParameters.aspect_ratio * m_splinesToCurvilinearParameters.average_width);

            for (UInt iter = 0; iter < numIterations; ++iter)
            {
                ComputeVelocitiesSubIntervals(s,
                                              startGridLineLeft,
                                              endGridLineLeft,
                                              numLeftHeights,
                                              numRightHeights,
                                              firstHeight,
                                              m_leftGridLineIndex,
                                              m_rightGridLineIndex,
                                              m_numPerpendicularFacesOnSubintervalAndEdge,
                                              m_edgeVelocities,
                                              hh0LeftMaxRatio);

                ComputeVelocitiesSubIntervals(s,
                                              startGridLineRight,
                                              endGridLineRight,
                                              numLeftHeights,
                                              numRightHeights,
                                              firstHeight,
                                              m_rightGridLineIndex,
                                              m_leftGridLineIndex,
                                              m_numPerpendicularFacesOnSubintervalAndEdge,
                                              m_edgeVelocities,
                                              hh0RightMaxRatio);
            }

            // re-evaluate if growing grid outside is needed

            if ((numLeftHeights == 0 && numRightHeights <= 1) ||
                (numRightHeights == 0 && numLeftHeights <= 1) ||
                (numLeftHeights == 1 && numRightHeights == 1))
            {
                m_splinesToCurvilinearParameters.grow_grid_outside = 1;
            }

            // left part
            UInt numNLeftExponential = 0;
            if (m_splinesToCurvilinearParameters.grow_grid_outside == 1)
            {
                numNLeftExponential = std::min(ComputeNumberExponentialLayers(hh0LeftMaxRatio), static_cast<UInt>(m_curvilinearParameters.n_refinement));
            }
            for (auto i = startGridLineLeft; i < endGridLineLeft; ++i)
            {
                m_numPerpendicularFacesOnSubintervalAndEdge(1, i) = numNLeftExponential;
            }

            // right part
            UInt numNRightExponential = 0;
            if (m_splinesToCurvilinearParameters.grow_grid_outside == 1)
            {
                numNRightExponential = std::min(ComputeNumberExponentialLayers(hh0RightMaxRatio), static_cast<UInt>(m_curvilinearParameters.n_refinement));
            }
            for (auto i = startGridLineRight; i < endGridLineRight; ++i)
            {
                m_numPerpendicularFacesOnSubintervalAndEdge(1, i) = numNRightExponential;
            }
        }

        // compute local grow factors
        for (UInt s = 0; s < m_splines->GetNumSplines(); s++)
        {
            if (m_numMSplines[s] < 1)
            {
                continue;
            }

            for (auto i = m_leftGridLineIndex[s]; i < m_rightGridLineIndex[s] + m_numMSplines[s]; ++i)
            {
                if (!m_gridLine[i].IsValid() || !m_gridLine[i + 1].IsValid() || m_numPerpendicularFacesOnSubintervalAndEdge(1, i) < 1)
                {
                    continue;
                }
                m_growFactorOnSubintervalAndEdge(1, i) = ComputeGrowFactor(i);
            }
        }
    }

    double CurvilinearGridFromSplines::ComputeGrowFactor(UInt splineIndex) const
    {

        // eheight m_gridHeights
        double aspectRatioGrowFactor = 1.0;
        auto heightDifference = ComputeTotalExponentialHeight(aspectRatioGrowFactor, m_edgeVelocities[splineIndex], m_numPerpendicularFacesOnSubintervalAndEdge(1, splineIndex)) - m_gridHeights(1, splineIndex);

        const double deps = 0.01;
        double aspectRatioGrowFactorIncremented = 1.0 + deps;
        auto heightDifferenceIncremented = ComputeTotalExponentialHeight(aspectRatioGrowFactorIncremented, m_edgeVelocities[splineIndex], m_numPerpendicularFacesOnSubintervalAndEdge(1, splineIndex)) - m_gridHeights(1, splineIndex);

        const UInt numIterations = 1000;
        const double relaxationFactor = 0.5;
        double oldAspectRatio;
        double oldHeightDifference = heightDifference;

        if (std::abs(heightDifferenceIncremented) > tolerance && std::abs(heightDifferenceIncremented - heightDifference) > tolerance)
        {
            for (UInt i = 0; i < numIterations; ++i)
            {
                oldAspectRatio = aspectRatioGrowFactor;
                oldHeightDifference = heightDifference;

                aspectRatioGrowFactor = aspectRatioGrowFactorIncremented;
                heightDifference = heightDifferenceIncremented;

                aspectRatioGrowFactorIncremented = aspectRatioGrowFactor - relaxationFactor * heightDifference /
                                                                               (heightDifference - oldHeightDifference + 1e-16) *
                                                                               (aspectRatioGrowFactor - oldAspectRatio);

                heightDifferenceIncremented = ComputeTotalExponentialHeight(aspectRatioGrowFactorIncremented,
                                                                            m_edgeVelocities[splineIndex],
                                                                            m_numPerpendicularFacesOnSubintervalAndEdge(1, splineIndex)) -
                                              m_gridHeights(1, splineIndex);

                if (std::abs(oldHeightDifference) < tolerance)
                {
                    break;
                }
            }
        }

        if (oldHeightDifference < tolerance)
        {
            return aspectRatioGrowFactorIncremented;
        }
        return 1.0;
    }

    double CurvilinearGridFromSplines::ComputeTotalExponentialHeight(double aspectRatio, double height, UInt numLayers) const
    {

        const auto absAspectRatio = aspectRatio - 1.0;
        if (absAspectRatio < -1e-8 || absAspectRatio > 1e-8)
        {
            return (std::pow(aspectRatio, static_cast<int>(numLayers)) - 1.0) / absAspectRatio * height;
        }
        return height * static_cast<double>(numLayers);
    }

    UInt CurvilinearGridFromSplines::ComputeNumberExponentialLayers(const double heightRatio) const
    {
        if (m_splinesToCurvilinearParameters.aspect_ratio_grow_factor - 1.0 > 1e-8)
        {
            return UInt(std::floor(std::log((m_splinesToCurvilinearParameters.aspect_ratio_grow_factor - 1.0) * heightRatio + 1.0) /
                                   log(m_splinesToCurvilinearParameters.aspect_ratio_grow_factor)));
        }
        return UInt(std::floor(0.999 + heightRatio));
    }

    void CurvilinearGridFromSplines::ComputeVelocitiesSubIntervals(UInt s,
                                                                   UInt startGridLineIndex,
                                                                   UInt endGridLineIndex,
                                                                   UInt numHeights,
                                                                   UInt numOtherSideHeights,
                                                                   const double firstHeight,
                                                                   const std::vector<UInt>& gridLineIndex,
                                                                   const std::vector<UInt>& otherGridLineIndex,
                                                                   lin_alg::Matrix<UInt>& numPerpendicularFacesOnSubintervalAndEdge,
                                                                   std::vector<double>& edgeVelocities,
                                                                   double& hh0MaxRatio)
    {

        hh0MaxRatio = 0.0;
        if ((numHeights > 1 && numHeights == numOtherSideHeights) || numHeights > numOtherSideHeights)
        {
            const auto maxHeight = *std::max_element(m_gridHeights.row(0).begin() + startGridLineIndex,
                                                     m_gridHeights.row(0).begin() + endGridLineIndex);

            auto numNUniformPart = static_cast<UInt>(std::floor(maxHeight / firstHeight + 0.99999));
            numNUniformPart = std::min(numNUniformPart, m_maxNUniformPart);

            for (auto i = startGridLineIndex; i < endGridLineIndex; ++i)
            {
                numPerpendicularFacesOnSubintervalAndEdge(0, i) = numNUniformPart;
                edgeVelocities[i] = m_gridHeights(0, i) / static_cast<double>(numNUniformPart);
                hh0MaxRatio = std::max(hh0MaxRatio, m_gridHeights(1, i) / edgeVelocities[i]);
            }
        }
        else
        {
            // only one subinterval: no uniform part
            const UInt numNUniformPart = 0;
            for (auto i = startGridLineIndex; i < endGridLineIndex; ++i)
            {
                numPerpendicularFacesOnSubintervalAndEdge(0, i) = numNUniformPart;
                edgeVelocities[i] = firstHeight;

                // compare with other side of spline
                const auto otherSideIndex = otherGridLineIndex[s] + m_numMSplines[s] - (i - gridLineIndex[s] + 1);

                if (edgeVelocities[otherSideIndex] != constants::missing::doubleValue)
                {
                    if (numPerpendicularFacesOnSubintervalAndEdge(0, otherSideIndex) == 0)
                    {
                        edgeVelocities[i] = std::max(edgeVelocities[i], edgeVelocities[otherSideIndex]);
                    }
                    else
                    {
                        edgeVelocities[i] = edgeVelocities[otherSideIndex];
                    }
                }
                for (UInt j = 1; j < m_maxNumCenterSplineHeights; ++j)
                {
                    m_gridHeights(j, i) = m_gridHeights(j - 1, i);
                }

                for (auto j = startGridLineIndex; j < endGridLineIndex; ++j)
                {

                    hh0MaxRatio = std::max(hh0MaxRatio, m_gridHeights(1, j) / edgeVelocities[j]);
                }
            }
        }
    }

    void CurvilinearGridFromSplines::ComputeGridHeights()
    {
        const auto numSplines = m_splines->GetNumSplines();
        lin_alg::ResizeAndFillMatrix(m_gridHeights,
                                     m_maxNumCenterSplineHeights,
                                     m_numM - 1,
                                     false,
                                     constants::missing::doubleValue);

        lin_alg::Matrix<double> heightsLeft(m_maxNumCenterSplineHeights, m_curvilinearParameters.m_refinement);
        heightsLeft.fill(0.0); // not necessary
        lin_alg::Matrix<double> heightsRight(m_maxNumCenterSplineHeights, m_curvilinearParameters.m_refinement);
        heightsRight.fill(0.0); // not necessary
        std::vector<double> edgesCenterPoints(m_numM, 0.0);
        std::vector<double> crossingSplinesDimensionalCoordinates(numSplines, 0.0);
        std::vector<double> localSplineDerivatives(numSplines, 0.0);
        std::vector<UInt> localValidSplineIndices(numSplines, 0);

        for (UInt s = 0; s < numSplines; s++)
        {
            if (m_type[s] != SplineTypes::central)
            {
                continue;
            }

            const auto numM = m_numMSplines[s];

            // Get the minimum number of sub-intervals in the cross splines for this center spline
            // can also use matrix.row(row_index).minCoeff()
            const auto minNumLeftIntervals = *std::min_element(m_numCrossSplineLeftHeights.row(s).begin(),
                                                               m_numCrossSplineLeftHeights.row(s).end());
            const auto minNumRightIntervals = *std::min_element(m_numCrossSplineRightHeights.row(s).begin(),
                                                                m_numCrossSplineRightHeights.row(s).end());

            std::fill(heightsLeft.row(0).begin(), heightsLeft.row(0).begin() + numM, m_maximumGridHeights[s]);
            std::fill(heightsRight.row(0).begin(), heightsRight.row(0).begin() + numM, m_maximumGridHeights[s]);

            if (m_numCrossingSplines[s] == 1)
            {
                // only one crossing spline present:
                for (UInt i = 0; i < minNumLeftIntervals; ++i)
                {
                    std::fill(heightsLeft.row(i).begin(), heightsLeft.row(i).begin() + numM, m_crossSplineRightHeights(s, i)[0]);
                }
                for (UInt i = 0; i < minNumRightIntervals; ++i)
                {
                    std::fill(heightsRight.row(i).begin(), heightsRight.row(i).begin() + numM, m_crossSplineLeftHeights(s, i)[0]);
                }
            }
            else
            {
                const auto leftGridLineIndex = m_leftGridLineIndex[s];
                edgesCenterPoints[0] = m_splines->ComputeSplineLength(s, 0, m_gridLineDimensionalCoordinates[leftGridLineIndex]);
                for (UInt i = 0; i < numM; ++i)
                {
                    edgesCenterPoints[i + 1] = edgesCenterPoints[i] +
                                               m_splines->ComputeSplineLength(s,
                                                                              m_gridLineDimensionalCoordinates[leftGridLineIndex + i],
                                                                              m_gridLineDimensionalCoordinates[leftGridLineIndex + i + 1]);
                }

                // compute at edge center points
                for (UInt i = 0; i < numM; ++i)
                {
                    edgesCenterPoints[i] = 0.5 * (edgesCenterPoints[i] + edgesCenterPoints[i + 1]);
                }
                edgesCenterPoints[numM] = constants::missing::doubleValue;

                // compute center spline path length of cross splines
                crossingSplinesDimensionalCoordinates[0] = m_splines->ComputeSplineLength(s, 0.0, m_crossSplineCoordinates(s, 0));
                for (UInt i = 0; i < m_numCrossingSplines[s] - 1; ++i)
                {
                    crossingSplinesDimensionalCoordinates[i + 1] = crossingSplinesDimensionalCoordinates[i] +
                                                                   m_splines->ComputeSplineLength(s, m_crossSplineCoordinates(s, i), m_crossSplineCoordinates(s, i + 1));
                }

                for (UInt j = 0; j < m_maxNumCenterSplineHeights; ++j)
                {

                    FindNearestCrossSplines(s,
                                            j,
                                            m_numCrossSplineLeftHeights,
                                            m_crossSplineLeftHeights,
                                            edgesCenterPoints,
                                            localValidSplineIndices,
                                            localSplineDerivatives,
                                            crossingSplinesDimensionalCoordinates,
                                            heightsLeft);

                    FindNearestCrossSplines(s,
                                            j,
                                            m_numCrossSplineRightHeights,
                                            m_crossSplineRightHeights,
                                            edgesCenterPoints,
                                            localValidSplineIndices,
                                            localSplineDerivatives,
                                            crossingSplinesDimensionalCoordinates,
                                            heightsRight);
                }
            }

            // store grid height
            for (UInt j = 0; j < m_maxNumCenterSplineHeights; ++j)
            {
                for (UInt i = 0; i < m_numMSplines[s]; ++i)
                {
                    m_gridHeights(j, m_leftGridLineIndex[s] + i) = heightsLeft(j, i);
                    m_gridHeights(j, m_rightGridLineIndex[s] + m_numMSplines[s] - i - 1) = heightsRight(j, i);
                }
            }
        }
    }

    void CurvilinearGridFromSplines::FindNearestCrossSplines(UInt s,
                                                             UInt j,
                                                             const lin_alg::Matrix<UInt>& numHeightsLeft,
                                                             const lin_alg::Matrix<std::vector<double>>& crossSplineLeftHeights,
                                                             const std::vector<double>& edgesCenterPoints,
                                                             std::vector<UInt>& localValidSplineIndices,
                                                             std::vector<double>& localSplineDerivatives,
                                                             std::vector<double>& crossingSplinesDimensionalCoordinates,
                                                             lin_alg::Matrix<double>& heights)
    {
        UInt numValid = 0;
        for (UInt i = 0; i < m_numCrossingSplines[s]; ++i)
        {
            if (numHeightsLeft(s, i) != 0)
            {
                localValidSplineIndices[numValid] = i;
                numValid++;
            }
        }

        // no sub-heights to compute
        if (numValid == 0)
        {
            return;
        }

        const auto numM = m_numMSplines[s];
        std::vector<double> localCornerPoints(numValid);

        // TODO: strided memory access
        for (UInt i = 0; i < numValid; ++i)
        {
            const auto index = localValidSplineIndices[i];
            localCornerPoints[i] = crossSplineLeftHeights(s, index)[j];
        }

        localSplineDerivatives = SplineAlgorithms::SecondOrderDerivative(localCornerPoints, 0, static_cast<UInt>(localCornerPoints.size()) - 1);

        crossingSplinesDimensionalCoordinates[0] = m_splines->ComputeSplineLength(s, 0.0, m_crossSplineCoordinates(s, 0));
        for (UInt i = 0; i < numM; ++i)
        {
            UInt leftIndex = 0;
            double leftCoordinate = crossingSplinesDimensionalCoordinates[localValidSplineIndices[leftIndex]];
            auto rightIndex = std::min(UInt(1), numValid - 1);
            double rightCoordinate = crossingSplinesDimensionalCoordinates[localValidSplineIndices[rightIndex]];
            // Find two nearest cross splines
            while (rightCoordinate < edgesCenterPoints[i] && rightIndex < numValid)
            {
                leftIndex = rightIndex;
                leftCoordinate = rightCoordinate;
                rightIndex++;
                rightCoordinate = crossingSplinesDimensionalCoordinates[localValidSplineIndices[rightIndex]];
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

            factor = std::max(std::min(double(leftIndex + 1) + factor - 1.0, double(numValid - 1)), 0.0);

            heights(j, i) = ComputePointOnSplineAtAdimensionalDistance(localCornerPoints, localSplineDerivatives, factor);
        }
    }

    void CurvilinearGridFromSplines::GetSplineIntersections(UInt splineIndex)
    {
        std::cout << " CurvilinearGridFromSplines::GetSplineIntersections " << splineIndex << std::endl;

        m_numCrossingSplines[splineIndex] = 0;
        const auto numSplines = m_splines->GetNumSplines();
        std::fill(m_crossingSplinesIndices.row(splineIndex).begin(), m_crossingSplinesIndices.row(splineIndex).end(), 0);
        std::fill(m_isLeftOriented.row(splineIndex).begin(), m_isLeftOriented.row(splineIndex).end(), true);
        std::fill(m_crossSplineCoordinates.row(splineIndex).begin(), m_crossSplineCoordinates.row(splineIndex).end(), std::numeric_limits<double>::max());
        std::fill(m_cosCrossingAngle.row(splineIndex).begin(), m_cosCrossingAngle.row(splineIndex).end(), constants::missing::doubleValue);

        for (UInt s = 0; s < numSplines; ++s)
        {
            // a crossing is a spline with 2 nodes and another with more than 2 nodes
            const auto numSplineNodesS = static_cast<UInt>(m_splines->m_splineNodes[s].size());
            const auto numSplineNodesI = static_cast<UInt>(m_splines->m_splineNodes[splineIndex].size());
            if ((numSplineNodesS == 2 && numSplineNodesI == 2) ||
                (numSplineNodesS > 2 && numSplineNodesI > 2))
            {
                continue;
            }

            double crossProductIntersection;
            Point intersectionPoint;
            double firstSplineRatio;
            double secondSplineRatio;
            bool crossing = m_splines->GetSplinesIntersection(splineIndex, s, crossProductIntersection, intersectionPoint, firstSplineRatio, secondSplineRatio);

            std::cout << " is crossing " << std::boolalpha << crossing << "  " << s << "  " << crossProductIntersection << "  "
                      << m_splinesToCurvilinearParameters.min_cosine_crossing_angles << "  "
                      << std::endl;

            if (std::abs(crossProductIntersection) < m_splinesToCurvilinearParameters.min_cosine_crossing_angles)
            {
                crossing = false;
            }

            std::cout << " is crossing (after) " << std::boolalpha << crossing << std::endl;

            if (crossing)
            {
                m_numCrossingSplines[splineIndex]++;
                m_crossingSplinesIndices(splineIndex, s) = s;
                if (crossProductIntersection > 0.0)
                {
                    m_isLeftOriented(splineIndex, s) = false;
                }
                m_crossSplineCoordinates(splineIndex, s) = firstSplineRatio;
                m_cosCrossingAngle(splineIndex, s) = crossProductIntersection;
            }
        }

        auto const sortedIndices = lin_alg::SortRow(m_crossSplineCoordinates.row(splineIndex));
        lin_alg::ReorderRow(m_crossSplineCoordinates.row(splineIndex), sortedIndices);
        lin_alg::ReorderRow(m_crossingSplinesIndices.row(splineIndex), sortedIndices);
        lin_alg::ReorderRow(m_isLeftOriented.row(splineIndex), sortedIndices);
    }

    void CurvilinearGridFromSplines::MakeAllGridLines()
    {
        std::cout << " CurvilinearGridFromSplines::MakeAllGridLines " << std::endl;
        m_numM = 0;
        UInt numCenterSplines = 0;
        for (UInt s = 0; s < m_splines->GetNumSplines(); ++s)
        {
            // center splines only
            if (m_type[s] != SplineTypes::central)
            {
                continue;
            }
            numCenterSplines += 1;
        }

        if (numCenterSplines == 0)
        {
            throw std::invalid_argument("CurvilinearGridFromSplines::MakeAllGridLines: There are no center splines.");
        }

        UInt gridLineIndex = 0;

        std::cout << " m_splines->GetNumSplines() " << m_splines->GetNumSplines() << std::endl;

        for (UInt s = 0; s < m_splines->GetNumSplines(); ++s)
        {
            // center splines only
            if (m_type[s] != SplineTypes::central)
            {
                continue;
            }

            // upper bound of m_gridLine size, with two sides of spline and two missing values added
            const auto sizeGridLine = gridLineIndex + 1 + 2 * (m_curvilinearParameters.m_refinement + 1) + 2;
            m_gridLine.resize(sizeGridLine);
            m_gridLineDimensionalCoordinates.resize(sizeGridLine);

            if (gridLineIndex > 0)
            {
                m_gridLine[gridLineIndex] = {constants::missing::doubleValue, constants::missing::doubleValue};
                m_gridLineDimensionalCoordinates[gridLineIndex] = constants::missing::doubleValue;
                gridLineIndex++;
            }

            m_leftGridLineIndex[s] = gridLineIndex;

            const UInt numM = MakeGridLine(s, gridLineIndex);

            gridLineIndex = gridLineIndex + numM + 1;
            m_gridLine[gridLineIndex] = Point();
            m_gridLineDimensionalCoordinates[gridLineIndex] = constants::missing::doubleValue;
            gridLineIndex++;

            // add other side of gridline
            m_rightGridLineIndex[s] = gridLineIndex;
            auto rightIndex = m_rightGridLineIndex[s] - 1;
            for (auto j = m_rightGridLineIndex[s] - 1; j >= m_leftGridLineIndex[s] && j != static_cast<UInt>(0) - 1; --j)
            {
                m_gridLine[rightIndex] = m_gridLine[j];
                m_gridLineDimensionalCoordinates[rightIndex] = m_gridLineDimensionalCoordinates[j];
                ++rightIndex;
            }

            gridLineIndex = rightIndex;

            m_numMSplines[s] = numM;
            m_numM = gridLineIndex;
        }

        std::cout << " end CurvilinearGridFromSplines::MakeAllGridLines " << std::endl;
    }

    UInt CurvilinearGridFromSplines::MakeGridLine(UInt splineIndex,
                                                  UInt startingIndex)
    {
        std::cout << " MakeGridLine " << splineIndex << "  " << startingIndex << std::endl;

        const auto endSplineAdimensionalCoordinate = static_cast<double>(m_splines->m_splineNodes[splineIndex].size()) - 1;
        const auto splineLength = m_splines->ComputeSplineLength(splineIndex,
                                                                 0.0,
                                                                 endSplineAdimensionalCoordinate,
                                                                 10,
                                                                 m_splinesToCurvilinearParameters.curvature_adapted_grid_spacing,
                                                                 m_maximumGridHeights[splineIndex]);

        // first estimation of nodes along m
        auto numM = 1 + static_cast<UInt>(std::floor(splineLength / m_splinesToCurvilinearParameters.average_width));
        // auto numM = 1 + static_cast<UInt>(std::floor(m_splines->m_splinesLength[splineIndex] / m_splinesToCurvilinearParameters.average_width));

        numM = std::min(numM, static_cast<UInt>(m_curvilinearParameters.m_refinement));

        m_gridLine[startingIndex] = m_splines->m_splineNodes[splineIndex][0];

        std::cout << "make grid lines: " << splineIndex << "  "
                  << splineLength << "  "
                  << m_splines->m_splinesLength[splineIndex] << "  "
                  << m_splines->m_splineNodes[splineIndex].size() << "  "
                  << m_splinesToCurvilinearParameters.average_width << "   "
                  << m_maximumGridHeights[splineIndex] << "  "
                  << startingIndex << "  "
                  << numM << "  "
                  << m_gridLine[startingIndex].x << ", "
                  << m_gridLine[startingIndex].y << "  "
                  << std::endl;

        auto currentMaxWidth = std::numeric_limits<double>::max();
        std::vector<double> distances(numM);
        while (currentMaxWidth > m_splinesToCurvilinearParameters.average_width)
        {
            currentMaxWidth = 0.0;

            for (UInt n = 0; n < numM; ++n)
            {
                distances[n] = splineLength * (n + 1.0) / static_cast<double>(numM);
            }

            auto [points, adimensionalDistances] = m_splines->ComputePointOnSplineFromAdimensionalDistance(splineIndex,
                                                                                                           m_maximumGridHeights[splineIndex],
                                                                                                           m_splinesToCurvilinearParameters.curvature_adapted_grid_spacing,
                                                                                                           distances);

            std::vector<double> computedDistances(distances.size());
            std::vector<Point> computedPoints(distances.size());

            Point start = m_splines->m_splineNodes[splineIndex][0];
            Point end;

            for (size_t i = 0; i < adimensionalDistances.size(); ++i)
            {
                end = m_splines->Evaluate(splineIndex, adimensionalDistances[i]);
                computedPoints[i] = end;
                computedDistances[i] = ComputeDistance(start, end, m_splines->m_projection);
                start = end;
            }

            // double h = 6.0 / adimensionalDistances.size();
            // points[0] = m_splines->Evaluate(splineIndex, adimensionalDistances[0]);

            // for (size_t i = 1; i < adimensionalDistances.size(); ++i)
            // {
            //     adimensionalDistances[i] = adimensionalDistances[i - 1] + h;
            //     points[i] = m_splines->Evaluate(splineIndex, adimensionalDistances[i]);
            // }

            std::cout << " distances[n] ";

            for (UInt n = 0; n < numM; ++n)
            {
                std::cout << computedDistances[n] << "   ";
            }

            std::cout << std::endl;
            std::cout << " distances[n] ";

            for (UInt n = 0; n < numM; ++n)
            {
                std::cout << distances[n] << "   ";
            }

            std::cout << std::endl;

            std::cout << " distances_step[n] ";

            for (UInt n = 1; n < numM; ++n)
            {
                std::cout << distances[n] - distances[n - 1] << "   ";
            }

            std::cout << std::endl;

            std::cout << " points ";

            for (size_t iii = 0; iii < points.size(); ++iii)
            {
                std::cout << points[iii].x << ", " << points[iii].y << " --- ";
            }

            std::cout << std::endl;

            std::cout << " points ";

            for (size_t iii = 0; iii < computedPoints.size(); ++iii)
            {
                std::cout << computedPoints[iii].x << ", " << computedPoints[iii].y << " --- ";
            }

            std::cout << std::endl;

            std::cout << " adimensionalDistances ";

            for (size_t iii = 0; iii < adimensionalDistances.size(); ++iii)
            {
                std::cout << adimensionalDistances[iii] << ",  ";
            }

            std::cout << std::endl;
            std::cout << std::endl;

            // If distances is initialised in size numM + 1 and with 0 in the zeroth position.
            // can this loop be moved outside the while loop, except the finding of th
            for (UInt n = 0; n < numM; ++n)
            {
                const auto index = startingIndex + n + 1;
                m_gridLineDimensionalCoordinates[index] = adimensionalDistances[n];
                m_gridLine[index] = points[n];
                currentMaxWidth = std::max(currentMaxWidth, ComputeDistance(m_gridLine[index - 1], m_gridLine[index], m_splines->m_projection));
            }

            // a gridline is computed
            if (currentMaxWidth < m_splinesToCurvilinearParameters.average_width ||
                numM == static_cast<UInt>(m_curvilinearParameters.m_refinement))
            {
                break;
            }

            // room for sub-division
            if (currentMaxWidth > m_splinesToCurvilinearParameters.average_width)
            {
                numM = std::min(std::max(static_cast<UInt>(m_curvilinearParameters.m_refinement / m_maximumGridHeights[splineIndex] * static_cast<double>(numM)),
                                         numM + static_cast<UInt>(1)),
                                static_cast<UInt>(m_curvilinearParameters.m_refinement));

                distances.resize(numM);
                adimensionalDistances.resize(numM);
                points.resize(numM);
            }
        }

        return numM;
    }

    void CurvilinearGridFromSplines::ComputeSplineProperties(const bool restoreOriginalProperties)
    {
        std::cout << " begin CurvilinearGridFromSplines::ComputeSplineProperties  " << std::endl;

        AllocateSplinesProperties();

        for (UInt s = 0; s < m_splines->GetNumSplines(); ++s)
        {
            GetSplineIntersections(s);
        }

        // select all non-cross splines only
        for (UInt s = 0; s < m_splines->GetNumSplines(); ++s)
        {
            m_type[s] = SplineTypes::crossing;
            // if more than 2 nodes, the spline must be a central spline
            if (m_splines->m_splineNodes[s].size() > 2)
            {
                m_type[s] = SplineTypes::central;
            }
        }
        // check the cross splines. The center spline is the middle spline that crosses the cross spline
        for (UInt s = 0; s < m_splines->GetNumSplines(); ++s)
        {
            // only crossing splines with one or more central spline
            if (m_splines->m_splineNodes[s].size() != 2 || m_numCrossingSplines[s] < 1)
            {
                continue;
            }

            auto middleCrossingSpline = std::min(m_numCrossingSplines[s] / 2, m_numCrossingSplines[s]);
            auto crossingSplineIndex = m_crossingSplinesIndices(s, middleCrossingSpline);

            // if m_numIntersectingSplines[s] is even, check if the middle spline has already been assigned as a bounding spline
            if (m_type[crossingSplineIndex] != SplineTypes::central && 2 * crossingSplineIndex == m_numCrossingSplines[s])
            {
                middleCrossingSpline = std::min(middleCrossingSpline + 1, m_numCrossingSplines[s] - 1);
                crossingSplineIndex = m_crossingSplinesIndices(s, middleCrossingSpline);
            }

            if (m_type[crossingSplineIndex] == SplineTypes::central)
            {
                // associate bounding splines with the middle spline
                for (UInt i = 0; i < middleCrossingSpline; ++i)
                {
                    const auto index = m_crossingSplinesIndices(s, i);
                    m_type[index] = SplineTypes::lateral; // lateral spline
                    m_centralSplineIndex[index] = -static_cast<int>(crossingSplineIndex);
                }
                for (auto i = middleCrossingSpline + 1; i < m_numCrossingSplines[s]; ++i)
                {
                    const auto index = m_crossingSplinesIndices(s, i);
                    m_type[index] = SplineTypes::lateral; // lateral spline
                    m_centralSplineIndex[index] = -static_cast<int>(crossingSplineIndex);
                }
            }
        }

        if (restoreOriginalProperties)
        {
            // restore original spline properties
            for (UInt s = 0; s < m_numOriginalSplines; ++s)
            {
                m_leftGridLineIndex[s] = m_leftGridLineIndexOriginal[s];
                m_rightGridLineIndex[s] = m_rightGridLineIndexOriginal[s];
                m_numMSplines[s] = m_mfacOriginal[s];
                m_maximumGridHeights[s] = m_maximumGridHeightsOriginal[s];
                m_type[s] = m_originalTypes[s];
            }

            // mark new splines as artificial cross splines
            for (UInt s = m_numOriginalSplines; s < m_splines->GetNumSplines(); ++s)
            {
                m_type[s] = SplineTypes::artificial;
            }
        }

        ComputeHeights();
        std::cout << " end CurvilinearGridFromSplines::ComputeSplineProperties  " << m_splines->GetNumSplines() << std::endl;
    }

    void CurvilinearGridFromSplines::ComputeHeights()
    {
        for (UInt i = 0; i < m_splines->GetNumSplines(); ++i)
        {
            // Heights should be computed only for center splines
            if (m_splines->m_splineNodes[i].size() <= 2)
            {
                continue;
            }
            for (UInt j = 0; j < m_numCrossingSplines[i]; ++j)
            {
                ComputeSubHeights(i, j);
            }
        }

        // compute m_maximumGridHeight
        for (UInt s = 0; s < m_splines->GetNumSplines(); ++s)
        {
            if (m_numCrossingSplines[s] == 0)
            {
                m_maximumGridHeights[s] = m_splinesToCurvilinearParameters.aspect_ratio * m_splines->m_splinesLength[s];
                std::cout << " m_maximumGridHeights[s] " << s << "  " << m_maximumGridHeights[s] << std::endl;
                continue;
            }

            std::cout << " m_numCrossingSplines[s] " << s << "   " << m_numCrossingSplines[s] << "  " << m_maximumGridHeights[s] << std::endl;

            double maximumHeight = 0.0;
            for (UInt c = 0; c < m_numCrossingSplines[s]; ++c)
            {
                double sumLeftHeights = 0.0;
                for (UInt ss = 0; ss < m_numCrossSplineLeftHeights(s, c); ++ss)
                {
                    sumLeftHeights += m_crossSplineLeftHeights(s, c)[ss];
                }
                double sumRightHeights = 0.0;
                for (UInt ss = 0; ss < m_numCrossSplineRightHeights(s, c); ++ss)
                {
                    sumRightHeights += m_crossSplineRightHeights(s, c)[ss];
                }
                maximumHeight = std::max(maximumHeight, std::max(sumLeftHeights, sumRightHeights));
            }

            m_maximumGridHeights[s] = maximumHeight;
        }
    }

    void CurvilinearGridFromSplines::ComputeSubHeights(UInt centerSplineIndex, UInt crossingSplineLocalIndex)
    {
        // find center spline index
        UInt centerSplineLocalIndex = 0;
        const auto crossingSplineIndex = m_crossingSplinesIndices(centerSplineIndex, crossingSplineLocalIndex); // js
        for (UInt s = 0; s < m_numCrossingSplines[crossingSplineIndex]; ++s)
        {
            if (m_crossingSplinesIndices(crossingSplineIndex, s) == centerSplineIndex)
            {
                centerSplineLocalIndex = s;
                break;
            }
        }

        // right part
        UInt numSubIntervalsRight = 0;
        UInt rightCenterSplineIndex = centerSplineLocalIndex;
        UInt leftCenterSplineIndex;
        m_crossSplineRightHeights(centerSplineIndex, crossingSplineLocalIndex).resize(m_maxNumCenterSplineHeights, 0);
        for (auto s = centerSplineLocalIndex; s < m_numCrossingSplines[crossingSplineIndex] - 1; ++s)
        {
            if (numSubIntervalsRight >= m_maxNumCenterSplineHeights)
            {
                break;
            }
            if (m_centralSplineIndex[m_crossingSplinesIndices(crossingSplineIndex, s + 1)] != -static_cast<int>(centerSplineIndex))
            {
                continue;
            }
            leftCenterSplineIndex = rightCenterSplineIndex;
            rightCenterSplineIndex = s + 1;
            m_crossSplineRightHeights(centerSplineIndex, crossingSplineLocalIndex)[numSubIntervalsRight] =
                m_splines->ComputeSplineLength(crossingSplineIndex,
                                               m_crossSplineCoordinates(crossingSplineIndex, leftCenterSplineIndex),
                                               m_crossSplineCoordinates(crossingSplineIndex, rightCenterSplineIndex));
            numSubIntervalsRight++;
        }

        const auto numSplineNodes = static_cast<UInt>(m_splines->m_splineNodes[crossingSplineIndex].size());
        m_crossSplineRightHeights(centerSplineIndex, crossingSplineLocalIndex)[numSubIntervalsRight] =
            m_splines->ComputeSplineLength(crossingSplineIndex,
                                           m_crossSplineCoordinates(crossingSplineIndex, rightCenterSplineIndex),
                                           static_cast<double>(numSplineNodes) - 1.0);
        numSubIntervalsRight++;
        std::fill(m_crossSplineRightHeights(centerSplineIndex, crossingSplineLocalIndex).begin() + numSubIntervalsRight,
                  m_crossSplineRightHeights(centerSplineIndex, crossingSplineLocalIndex).end(), 0.0);

        m_numCrossSplineRightHeights(centerSplineIndex, crossingSplineLocalIndex) = numSubIntervalsRight;

        // left part
        UInt numSubIntervalsLeft = 0;
        leftCenterSplineIndex = centerSplineLocalIndex;
        m_crossSplineLeftHeights(centerSplineIndex, crossingSplineLocalIndex).resize(m_maxNumCenterSplineHeights, 0);
        for (auto s = centerSplineLocalIndex; s >= 1; --s)
        {
            if (numSubIntervalsLeft >= m_maxNumCenterSplineHeights)
            {
                break;
            }
            if (m_centralSplineIndex[m_crossingSplinesIndices(crossingSplineIndex, s - 1)] != -static_cast<int>(centerSplineIndex))
            {
                continue;
            }
            rightCenterSplineIndex = leftCenterSplineIndex;
            leftCenterSplineIndex = s - 1;
            m_crossSplineLeftHeights(centerSplineIndex, crossingSplineLocalIndex)[numSubIntervalsLeft] =
                m_splines->ComputeSplineLength(crossingSplineIndex,
                                               m_crossSplineCoordinates(crossingSplineIndex, leftCenterSplineIndex),
                                               m_crossSplineCoordinates(crossingSplineIndex, rightCenterSplineIndex));
            numSubIntervalsLeft++;
        }

        m_crossSplineLeftHeights(centerSplineIndex, crossingSplineLocalIndex)[numSubIntervalsLeft] =
            m_splines->ComputeSplineLength(crossingSplineIndex,
                                           0.0,
                                           m_crossSplineCoordinates(crossingSplineIndex, leftCenterSplineIndex));

        numSubIntervalsLeft++;
        std::fill(m_crossSplineLeftHeights(centerSplineIndex, crossingSplineLocalIndex).begin() + numSubIntervalsLeft,
                  m_crossSplineLeftHeights(centerSplineIndex, crossingSplineLocalIndex).end(), 0.0);

        m_numCrossSplineLeftHeights(centerSplineIndex, crossingSplineLocalIndex) = numSubIntervalsLeft;

        // if not left oriented, swap
        if (!m_isLeftOriented(centerSplineIndex, crossingSplineLocalIndex))
        {
            m_numCrossSplineLeftHeights(centerSplineIndex, crossingSplineLocalIndex) = numSubIntervalsRight;
            m_numCrossSplineRightHeights(centerSplineIndex, crossingSplineLocalIndex) = numSubIntervalsLeft;

            const std::vector<double> leftSubIntervalsTemp(m_crossSplineLeftHeights(centerSplineIndex, crossingSplineLocalIndex));
            m_crossSplineLeftHeights(centerSplineIndex, crossingSplineLocalIndex) = m_crossSplineRightHeights(centerSplineIndex, crossingSplineLocalIndex);
            m_crossSplineRightHeights(centerSplineIndex, crossingSplineLocalIndex) = leftSubIntervalsTemp;
        }
    }

} // namespace meshkernel
