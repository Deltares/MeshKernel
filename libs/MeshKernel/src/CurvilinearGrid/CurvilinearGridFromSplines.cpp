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

namespace meshkernel
{
    void CurvilinearGridFromSplines::RootFinder(const std::array<double, MaxDegreeP1>& coefficients,
                                                int& degree,
                                                std::array<double, MaxDegree>& realRoots,
                                                std::array<double, MaxDegree>& imaginaryRoots) const
    {
        std::ranges::fill(realRoots, 0.0);
        std::ranges::fill(imaginaryRoots, 0.0);
        degree = 0;

        if (coefficients[0] != 0.0)
        {
            const double a = coefficients[0];
            const double b = coefficients[1];
            const double c = coefficients[2];

            const double discriminant = b * b - 4.0 * a * c;

            if (discriminant >= 0.0)
            {
                realRoots[1] = (-b + std::sqrt(discriminant)) / (2.0 * a);
                realRoots[0] = (-b - std::sqrt(discriminant)) / (2.0 * a);
            }
            else
            {
                const std::complex<double> r1 = (-b + std::sqrt(std::complex<double>(discriminant, 0.0))) / (2.0 * a);
                const std::complex<double> r0 = (-b - std::sqrt(std::complex<double>(discriminant, 0.0))) / (2.0 * a);

                realRoots[1] = r1.real();
                realRoots[0] = r0.real();

                imaginaryRoots[1] = r1.imag();
                imaginaryRoots[0] = r0.imag();
            }

            degree = 2;
        }
        else if (coefficients[0] != 0.0)
        {
            // linear or constant
            realRoots[0] = -coefficients[1] / coefficients[0];
            degree = 1;
        }
    }

    CurvilinearGridFromSplines::CurvilinearGridFromSplines(std::shared_ptr<Splines> splines,
                                                           const CurvilinearParameters& curvilinearParameters,
                                                           const SplinesToCurvilinearParameters& splinesToCurvilinearParameters)
        : m_splines(splines),
          m_curvilinearParameters(curvilinearParameters),
          m_splinesToCurvilinearParameters(splinesToCurvilinearParameters)
    {
        m_splinesToCurvilinearParameters.curvature_adapted_grid_spacing = 0;
        CheckCurvilinearParameters(curvilinearParameters);
        CheckSplinesToCurvilinearParameters(splinesToCurvilinearParameters);

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

        m_maximumGridHeights.resize(numSplines, constants::missing::doubleValue);
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
        for (auto layer = 1; layer <= m_curvilinearParameters.n_refinement; ++layer)
        {
            Iterate(layer);
        }

        if (m_splinesToCurvilinearParameters.remove_skinny_triangles == 1)
        {
            DeleteSkinnyTriangles();
        }

        return ComputeCurvilinearGridFromGridPoints();
    }

    UInt CurvilinearGridFromSplines::MoveGridNodes(const UInt i, const UInt j, const UInt firstLeftIndex, const UInt firstRightIndex)
    {
        constexpr double maxCosine = 0.93969;
        constexpr double squaredDistanceTolerance = 1e-4;
        constexpr double cosineTolerance = 1e-2;
        UInt numChanged = 0;

        const auto squaredCurrentDistance = ComputeSquaredDistance(m_gridPoints(j + 1, i),
                                                                   m_gridPoints(j + 1, firstRightIndex),
                                                                   m_splines->m_projection);
        const auto currentCosPhi = NormalizedInnerProductTwoSegments(m_gridPoints(j + 1, i),
                                                                     m_gridPoints(j, i),
                                                                     m_gridPoints(j + 1, i),
                                                                     m_gridPoints(j, firstRightIndex),
                                                                     m_splines->m_projection);

        if (squaredCurrentDistance < squaredDistanceTolerance && currentCosPhi > maxCosine)
        {

            // determine persistent node
            const auto leftCosPhi = NormalizedInnerProductTwoSegments(m_gridPoints(j - 1, i),
                                                                      m_gridPoints(j, i),
                                                                      m_gridPoints(j, i),
                                                                      m_gridPoints(j + 1, i),
                                                                      m_splines->m_projection);

            const auto rightCosPhi = NormalizedInnerProductTwoSegments(m_gridPoints(j - 1, firstRightIndex),
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

        return numChanged;
    }

    void CurvilinearGridFromSplines::DeleteSkinnyTriangles()
    {
        constexpr UInt numMaxIterations = 10;
        const UInt numN = static_cast<UInt>(m_gridPoints.rows()) - 2;
        constexpr double squaredDistanceTolerance = 1e-4;

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
                        numChanged += MoveGridNodes(i, j, firstLeftIndex, firstRightIndex);
                    }
                }

                if (numChanged == 0)
                {
                    break;
                }
            }
        }
    }

    std::tuple<Point, Point> CurvilinearGridFromSplines::GetCrossSplinePoints(const UInt s, const UInt i) const
    {
        const auto normal = NormalVectorOutside(m_gridLine[i], m_gridLine[i + 1], m_splines->m_projection);

        Point middle = 0.5 * (m_gridLine[i] + m_gridLine[i + 1]);
        Point x1{constants::missing::doubleValue, constants::missing::doubleValue};
        Point x2{constants::missing::doubleValue, constants::missing::doubleValue};

        if (m_splines->m_projection == Projection::cartesian)
        {
            x1 = middle - 2.0 * m_maximumGridHeights[s] * normal;
            x2 = middle + 2.0 * m_maximumGridHeights[s] * normal;
        }

        if (m_splines->m_projection == Projection::spherical)
        {
            const double factor = 1.0 / (constants::geometric::earth_radius * constants::conversion::degToRad);
            x1 = middle - (2.0 * m_maximumGridHeights[s] * factor) * normal;
            x2 = middle + (2.0 * m_maximumGridHeights[s] * factor) * normal;
        }

        return {x1, x2};
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

        for (UInt s = 0; s < m_numOriginalSplines; ++s)
        {
            // mirror only center splines
            if (m_type[s] != SplineTypes::central)
            {
                continue;
            }

            // construct the cross splines through the edges, along m discretization
            for (auto i = m_leftGridLineIndex[s]; i < m_leftGridLineIndex[s] + m_numMSplines[s]; ++i)
            {
                std::tie(newCrossSpline[0], newCrossSpline[1]) = GetCrossSplinePoints(s, i);

                std::cout << "intersection point: " << s << "  " << i << "  {" << newCrossSpline[0].x << ", " << newCrossSpline[0].y << "}, {"
                          << newCrossSpline[1].x << ", " << newCrossSpline[1].y << "}" << std::endl;

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

            std::cout << "grid line: (0, " << n << " ) = " << m_gridLine[n].x << ", " << m_gridLine[n].y << std::endl;

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
                m_validFrontNodes[i] = 0;
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

    bool CurvilinearGridFromSplines::SegmentCrossesCentreSpline(const Point& x1, const Point& x2) const
    {
        UInt centralSplineIndex = GetCentralSplineIndex();

        if (centralSplineIndex == constants::missing::uintValue)
        {
            throw AlgorithmError("No central spline found");
        }

        const std::vector<Point>& splinePoints = m_splines->m_splineNodes[centralSplineIndex];

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

    void CurvilinearGridFromSplines::SolveQuadratic(const std::array<double, MaxDegreeP1>& coefficients,
                                                    std::array<double, MaxDegree>& roots) const
    {
        int degree = MaxDegree;
        std::array<double, MaxDegreeP1> coeffs(coefficients);
        std::array<double, MaxDegree> realRoots{0.0, 0.0};
        std::array<double, MaxDegree> imagRoots{0.0, 0.0};
        roots.fill(constants::missing::doubleValue);

        for (UInt i = MaxDegree; i >= 1; --i)
        {

            if (std::abs(coeffs[MaxDegree - i]) > tolerance)
            {
                degree = i;
                break;
            }
        }

        // shuffle coefficients up if necessary
        if (degree < MaxDegree)
        {
            for (int j = 0; j <= degree; ++j)
            {
                coeffs[j] = coeffs[MaxDegree - degree + j];
            }

            for (int j = degree + 1; j < MaxDegreeP1; ++j)
            {
                coeffs[j] = 0.0;
            }
        }

        RootFinder(coeffs, degree, realRoots, imagRoots);

        if (degree <= 0)
        {
            return;
        }

        for (UInt i = 0; i < static_cast<UInt>(degree); ++i)
        {
            if (std::abs(imagRoots[i]) < 1.0e-4)
            {
                roots[i] = realRoots[i];
            }
        }
    }

    double CurvilinearGridFromSplines::ComputeNodeSegmentCrossingTime(const Point& x1, const Point& x3, const Point& x4,
                                                                      const Point& v1, const Point& v3, const Point& v4) const
    {
        // Extend limits by a relatively small amount
        const double lowerLimit = -8.0 * std::numeric_limits<double>::epsilon();
        const double upperLimit = 1.0 + 8.0 * std::numeric_limits<double>::epsilon();

        // Compute node differences
        const Point deltaPos13 = x3 - x1;
        const Point deltaPos34 = x4 - x3;
        const double cross13And34 = cross(deltaPos13, deltaPos34);

        if (cross13And34 < 0.0)
        {
            return 1.0e99;
        }

        // Compute velocity differences
        const Point deltaVel13 = v3 - v1;
        const Point deltaVel34 = v4 - v3;

        const double a = cross(deltaVel13, deltaVel34);
        const double b = cross(deltaPos13, deltaVel34) - cross(deltaPos34, deltaVel13);
        const double c = cross13And34;

        std::array<double, MaxDegreeP1> coeffs{a, b, c};
        std::array<double, MaxDegree> roots{constants::missing::doubleValue, constants::missing::doubleValue};
        std::array<double, MaxDegree> beta{constants::missing::doubleValue, constants::missing::doubleValue};

        SolveQuadratic(coeffs, roots);

        for (UInt i = 0; i < MaxDegree; ++i)
        {

            if (roots[i] == constants::missing::doubleValue || roots[i] < tolerance)
            {
                continue;
            }

            Point displacedPoint = deltaPos34 + deltaVel34 * roots[i];
            double det = lengthSquared(displacedPoint);

            if (std::abs(det) > tolerance)
            {
                beta[i] = -dot(deltaPos13 + deltaVel13 * roots[i], displacedPoint) / det;
            }
        }

        double time = 1.0e99;

        for (UInt i = 0; i < MaxDegree; ++i)
        {

            if (beta[i] >= lowerLimit && beta[i] <= upperLimit && roots[i] >= 0.0 && roots[i] != constants::missing::doubleValue)
            {
                time = std::min(time, roots[i]);
            }
        }

        return (time == constants::missing::doubleValue || time <= 0.0) ? 1.0e99 : time;
    }

    bool CurvilinearGridFromSplines::IncludeDirectNeighbours(const UInt layerIndex,
                                                             const UInt loopIndex,
                                                             const lin_alg::Matrix<UInt>& gridPointsIndices,
                                                             const UInt indexLeftOfLeft,
                                                             const UInt indexRightOfRight) const
    {

        UInt j1 = gridPointsIndices(loopIndex, 1);
        UInt j2 = gridPointsIndices(loopIndex + 1, 1);

        // Common check for layer index validity
        if (j1 < (layerIndex - 1) || j2 < (layerIndex - 1))
        {
            return true; // If layer indices are out of bounds, the neighbors should be included
        }

        UInt i1 = gridPointsIndices(loopIndex, 0);
        UInt i2 = gridPointsIndices(loopIndex + 1, 0);

        // Check index range conditions based on the relationship between indexRightOfRight and indexLeftOfLeft
        const bool isWithinRange = (indexRightOfRight >= indexLeftOfLeft)
                                       ? ((i1 > indexLeftOfLeft && i1 < indexRightOfRight) || (i2 > indexLeftOfLeft && i2 < indexRightOfRight))
                                       : (!(i1 >= indexRightOfRight && i1 <= indexLeftOfLeft) || !(i2 >= indexRightOfRight && i2 <= indexLeftOfLeft));

        // If the points are within the range, exclude the direct neighbors
        return !isWithinRange;
    }

    void CurvilinearGridFromSplines::ComputeMaximumTimeStep(const UInt layerIndex,
                                                            const lin_alg::RowVector<Point>& activeLayerPoints,
                                                            const std::vector<Point>& velocityVectorAtGridPoints,
                                                            const std::vector<Point>& frontGridPoints,
                                                            const std::vector<Point>& frontVelocities,
                                                            const lin_alg::Matrix<UInt>& gridPointsIndices,
                                                            const double timeStep,
                                                            double& otherTimeStep,
                                                            std::vector<double>& otherTimeStepMax) const
    {
        double maximumTimeStep = timeStep;

        std::vector<double> tmax(activeLayerPoints.size(), timeStep);

        for (UInt i = 0; i < activeLayerPoints.size() - 1; ++i)
        {
            if (!activeLayerPoints[i].IsValid() || !activeLayerPoints[i + 1].IsValid())
            {
                continue;
            }

            const Point x1 = activeLayerPoints[i];
            const Point x2 = activeLayerPoints[i + 1];

            const Point v1 = velocityVectorAtGridPoints[i];
            const Point v2 = velocityVectorAtGridPoints[i + 1];

            const double dL1 = ComputeDistance(x1, activeLayerPoints[i + 1], m_splines->m_projection);
            UInt indexLeft;
            UInt indexLeftOfLeft;
            UInt indexRight;
            UInt indexRightOfRight;
            UInt dummy;

            std::tie(indexLeft, dummy) = GetNeighbours(activeLayerPoints, i);
            std::tie(dummy, indexRight) = GetNeighbours(activeLayerPoints, i + 1);

            std::tie(indexLeftOfLeft, dummy) = GetNeighbours(activeLayerPoints, indexLeft);
            std::tie(dummy, indexRightOfRight) = GetNeighbours(activeLayerPoints, indexRight);

            for (UInt j = 0; j < frontGridPoints.size() - 1; ++j)
            {
                if (!frontGridPoints[j].IsValid() || !frontGridPoints[j + 1].IsValid())
                {
                    continue;
                }

                const Point x3 = frontGridPoints[j];
                const Point x4 = frontGridPoints[j + 1];

                const Point v3 = frontVelocities[j];
                const Point v4 = frontVelocities[j + 1];

                const double dL2 = ComputeDistance(x3, x4, m_splines->m_projection);

                const double d1 = ComputeDistance(x1, x3, m_splines->m_projection);
                const double d2 = ComputeDistance(x2, x3, m_splines->m_projection);
                const double d3 = ComputeDistance(x1, x4, m_splines->m_projection);
                const double d4 = ComputeDistance(x2, x4, m_splines->m_projection);

                if (d1 < tolerance || d2 < tolerance || d3 < tolerance || d4 < tolerance)
                {
                    continue;
                }

                if (!IncludeDirectNeighbours(layerIndex, j, gridPointsIndices, indexLeftOfLeft, indexRightOfRight))
                {
                    continue;
                }

                const double dmin = std::min({d1, d2, d3, d4});

                // get a lower bound for the cross time
                const double hlow2 = 0.25 * std::max(dmin * dmin - std::pow(0.5 * std::max(dL1, dL2), 2), 0.0);

                // check if the lower bounds is larger than the minimum found so far
                const double vv1 = std::sqrt(lengthSquared(v3 - v1));
                const double vv2 = std::sqrt(lengthSquared(v3 - v2));
                const double vv3 = std::sqrt(lengthSquared(v4 - v1));
                const double vv4 = std::sqrt(lengthSquared(v4 - v2));

                const double maxvv = std::max(std::max(vv1, vv2), std::max(vv3, vv4));

                if (std::sqrt(hlow2) > maxvv * std::min(tmax[i], tmax[i + 1]))
                {
                    continue;
                }

                const double t1 = ComputeNodeSegmentCrossingTime(x1, x3, x4, v1, v3, v4);
                const double t2 = ComputeNodeSegmentCrossingTime(x2, x3, x4, v2, v3, v4);
                const double t3 = ComputeNodeSegmentCrossingTime(x3, x1, x2, v3, v1, v2);
                const double t4 = ComputeNodeSegmentCrossingTime(x4, x1, x2, v4, v1, v2);

                const double tmin1234 = std::min(std::min(t1, t2), std::min(t3, t4));

                if (t1 == tmin1234)
                {
                    otherTimeStepMax[i] = std::min(otherTimeStepMax[i], tmin1234);
                    maximumTimeStep = std::min(maximumTimeStep, otherTimeStepMax[i]);
                }
                else if (t2 == tmin1234)
                {
                    otherTimeStepMax[i + 1] = std::min(otherTimeStepMax[i + 1], tmin1234);
                    maximumTimeStep = std::min(maximumTimeStep, otherTimeStepMax[i + 1]);
                }
                else if (t3 == tmin1234 || t4 == tmin1234)
                {
                    otherTimeStepMax[i] = std::min(otherTimeStepMax[i], tmin1234);
                    otherTimeStepMax[i + 1] = std::min(otherTimeStepMax[i + 1], tmin1234);
                    maximumTimeStep = std::min(maximumTimeStep, otherTimeStepMax[i]);
                    maximumTimeStep = std::min(maximumTimeStep, otherTimeStepMax[i + 1]);
                }

                if (tmin1234 == 0.0)
                {
                    break;
                }
            }
        }

        otherTimeStep = maximumTimeStep;
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

    void CurvilinearGridFromSplines::TranslateActiveLayerPoints(const std::vector<Point>& velocityVectorAtGridPoints,
                                                                const double localTimeStep,
                                                                lin_alg::RowVector<Point>& activeLayerPoints) const
    {
        std::cout << " CurvilinearGridFromSplines::TranslateActiveLayerPoints " << std::endl;

        for (UInt i = 0; i < velocityVectorAtGridPoints.size(); ++i)
        {
            if (m_validFrontNodes[i] == 1 && velocityVectorAtGridPoints[i].IsValid())
            {
                activeLayerPoints[i] += localTimeStep * velocityVectorAtGridPoints[i];
                std::cout << "activeLayerPoints[i] " << activeLayerPoints[i].x << ", " << activeLayerPoints[i].y << std::endl;
            }
            else
            {
                activeLayerPoints[i].x = constants::missing::doubleValue;
                activeLayerPoints[i].y = constants::missing::doubleValue;
            }
        }
    }

    void CurvilinearGridFromSplines::DisableValidFrontNodes(const std::vector<double>& otherTimeStepMax,
                                                            const double otherTimeStep,
                                                            const double localTimeStep,
                                                            std::vector<UInt>& newValidFrontNodes) const
    {
        if (otherTimeStep < localTimeStep)
        {
            for (UInt i = 0; i < newValidFrontNodes.size(); ++i)
            {
                if (otherTimeStepMax[i] - otherTimeStep < tolerance && (localTimeStep - otherTimeStepMax[i]) > tolerance)
                {
                    newValidFrontNodes[i] = 0;
                }
            }
        }

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
    }

    void CurvilinearGridFromSplines::InvalidateActiveLayerPoints(lin_alg::RowVector<Point>& activeLayerPoints) const
    {

        for (UInt i = 0; i < m_validFrontNodes.size(); ++i)
        {
            if (m_validFrontNodes[i] == 0)
            {
                activeLayerPoints[i] = {constants::missing::doubleValue, constants::missing::doubleValue};
            }
        }
    }

    void CurvilinearGridFromSplines::InvalidateGridNodes(const UInt layerIndex,
                                                         lin_alg::RowVector<Point>& activeLayerPoints,
                                                         std::vector<UInt>& newValidFrontNodes)
    {
        if (layerIndex < 2)
        {
            return;
        }

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
                const auto [currentLeftIndex, currentRightIndex] = GetNeighbours(activeLayerPoints, i);

                for (auto j = currentLeftIndex + 1; j < currentRightIndex; ++j)
                {
                    newValidFrontNodes[j] = 0;
                    m_gridPoints(layerIndex - 1, j) = {constants::missing::doubleValue, constants::missing::doubleValue};
                }
            }
        }
    }

    void CurvilinearGridFromSplines::GrowLayer(UInt layerIndex)
    {
        auto velocityVectorAtGridPoints = ComputeVelocitiesAtGridPoints(layerIndex - 1);

        lin_alg::RowVector<Point> activeLayerPoints(m_gridPoints.row(layerIndex - 1));

        for (UInt m = 0; m < velocityVectorAtGridPoints.size(); ++m)
        {
            if (!velocityVectorAtGridPoints[m].IsValid())
            {
                m_gridPoints(layerIndex - 1, m) = {constants::missing::doubleValue, constants::missing::doubleValue};
                activeLayerPoints[m] = {constants::missing::doubleValue, constants::missing::doubleValue};
            }
        }

        auto [frontVelocities, newFrontPoints, gridPointIndices] = CopyVelocitiesToFront(layerIndex - 1, velocityVectorAtGridPoints);

        double totalTimeStep = 0.0;
        double localTimeStep = 0.0;
        const auto numGridPoints = static_cast<UInt>(m_gridPoints.size());
        std::vector<UInt> newValidFrontNodes(numGridPoints);

        while (totalTimeStep < m_timeStep)
        {

            // Copy old front velocities
            newValidFrontNodes = m_validFrontNodes;
            InvalidateActiveLayerPoints(activeLayerPoints);

            const auto maximumGridLayerGrowTime = ComputeMaximumEdgeGrowTime(activeLayerPoints, velocityVectorAtGridPoints);
            localTimeStep = std::min(m_timeStep - totalTimeStep, *std::min_element(maximumGridLayerGrowTime.begin(), maximumGridLayerGrowTime.end()));
            double otherTimeStep = std::numeric_limits<double>::max();
            std::vector<double> otherTimeStepMax(newValidFrontNodes.size(), 1.0 + localTimeStep);

            if (m_splinesToCurvilinearParameters.check_front_collisions)
            {
                ComputeMaximumTimeStep(layerIndex,
                                       activeLayerPoints, velocityVectorAtGridPoints,
                                       newFrontPoints, frontVelocities, gridPointIndices,
                                       localTimeStep + 1.0, otherTimeStep, otherTimeStepMax);
            }

            DisableValidFrontNodes(otherTimeStepMax, otherTimeStep, localTimeStep, newValidFrontNodes);

            m_validFrontNodes = newValidFrontNodes;
            localTimeStep = std::min(localTimeStep, otherTimeStep);

            TranslateActiveLayerPoints(velocityVectorAtGridPoints, localTimeStep, activeLayerPoints);

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

                std::tie(frontVelocities, newFrontPoints, gridPointIndices) = CopyVelocitiesToFront(layerIndex, velocityVectorAtGridPoints);
            }
        }

        UInt numFrontPoints;
        std::tie(gridPointIndices, newFrontPoints, numFrontPoints) = FindFront();

        InvalidateGridNodes(layerIndex, activeLayerPoints, newValidFrontNodes);

        m_validFrontNodes = newValidFrontNodes;
    }

    std::vector<double> CurvilinearGridFromSplines::ComputeMaximumEdgeGrowTime(const lin_alg::RowVector<Point>& coordinates,
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

    bool CurvilinearGridFromSplines::CheckCornerNodes(const UInt p,
                                                      const lin_alg::RowVector<UInt>& indices,
                                                      const lin_alg::Matrix<UInt>& gridPointsIndices,
                                                      const bool increaseOffset) const
    {

        bool result = indices[0] == gridPointsIndices(p, 0) + (increaseOffset ? 1 : -1) &&
                      indices[1] == gridPointsIndices(p, 1) &&
                      m_validFrontNodes[indices[0]] == constants::missing::uintValue;

        result = result || (indices[0] == gridPointsIndices(p, 0) && indices[1] < gridPointsIndices(p, 1));

        return result;
    }

    std::tuple<std::vector<Point>, std::vector<Point>, lin_alg::Matrix<UInt>> CurvilinearGridFromSplines::CopyVelocitiesToFront(UInt layerIndex,
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
                const lin_alg::RowVector<UInt> previousIndices = gridPointsIndices.row(previous);
                const auto next = std::min(p + 1, numFrontPoints);
                const lin_alg::RowVector<UInt> nextIndices = gridPointsIndices.row(next);

                bool ll = CheckCornerNodes(p, previousIndices, gridPointsIndices, false /* sign */);
                bool lr = CheckCornerNodes(p, nextIndices, gridPointsIndices, true /* sign */);

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

        return {velocities, frontGridPoints, gridPointsIndices};
    }

    std::tuple<lin_alg::Matrix<UInt>,
               std::vector<Point>,
               UInt>
    CurvilinearGridFromSplines::FindFront() const
    {
        UInt const numGridPoints = static_cast<UInt>(m_gridPoints.size());
        lin_alg::Matrix<UInt> gridPointsIndices(numGridPoints, 2);
        gridPointsIndices.fill(constants::missing::uintValue);
        std::vector<Point> frontGridPoints(numGridPoints, Point());
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
                for (auto i = previousFrontPosition - 1; i >= currentFrontPosition; --i)
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

                std::cout << " 1 normalVectorLeft "
                          << m_gridPoints(layerIndex, m).x << ", " << m_gridPoints(layerIndex, m).y << " --- "
                          << m_gridPoints(layerIndex, currentLeftIndex).x << ", " << m_gridPoints(layerIndex, currentLeftIndex).y << " --- "
                          << m_gridPoints(layerIndex, currentRightIndex).x << ", " << m_gridPoints(layerIndex, currentRightIndex).y << " --- "
                          << normalVectorLeft.x << ", " << normalVectorLeft.y << " --- "
                          << std::endl;
            }
            else
            {
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

                std::cout << " 2 normalVectorLeft "
                          << m_gridPoints(layerIndex, m).x << ", " << m_gridPoints(layerIndex, m).y << " --- "
                          << m_gridPoints(layerIndex, currentLeftIndex).x << ", " << m_gridPoints(layerIndex, currentLeftIndex).y << " --- "
                          << m_gridPoints(layerIndex, currentRightIndex).x << ", " << m_gridPoints(layerIndex, currentRightIndex).y << " --- "
                          << normalVectorLeft.x << ", " << normalVectorLeft.y << " --- "
                          << normalVectorRight.x << ", " << normalVectorRight.y << " --- "
                          << std::endl;
            }

            if (currentLeftIndex == velocityVector.size() - 1)
            {
                continue;
            }

            const auto cosphi = dot(normalVectorLeft, normalVectorRight);

            if (cosphi < -1.0 + cosTolerance)
            {
                continue;
            }
            const auto leftVelocity = normalVectorLeft * m_edgeVelocities[currentLeftIndex];
            const auto rightVelocity = normalVectorRight * m_edgeVelocities[currentRightIndex - 1];
            const double rightLeftVelocityRatio = m_edgeVelocities[currentRightIndex - 1] / m_edgeVelocities[currentLeftIndex];

            if (rightLeftVelocityRatio - cosphi > eps && 1.0 / rightLeftVelocityRatio - cosphi > eps)
            {
                velocityVector[m] = (leftVelocity * (1.0 - rightLeftVelocityRatio * cosphi) +
                                     rightVelocity * (1.0 - 1.0 / rightLeftVelocityRatio * cosphi)) /
                                    (1.0 - cosphi * cosphi);
            }
            else if (cosphi - rightLeftVelocityRatio > eps)
            {
                velocityVector[m] = leftVelocity * rightLeftVelocityRatio / cosphi;
            }
            else
            {
                velocityVector[m] = rightVelocity * 1.0 / (rightLeftVelocityRatio * cosphi);
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
                if (localRightIndex + 1 == static_cast<int>(gridPoints.size()))
                {
                    break;
                }
            }
            else if (localRightIndex + 1 == static_cast<int>(gridPoints.size()))
            {
                localRightIndex = -1;
                circularConnection = false;
            }

            if (localRightIndex + 1 == static_cast<int>(gridPoints.size()) || !gridPoints[localRightIndex + 1].IsValid())
            {
                break;
            }
            localRightIndex++;
        }

        std::cout << " get neighbours: " << localLeftIndex << "  " << localRightIndex << std::endl;

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

            if (std::abs(crossProductIntersection) < m_splinesToCurvilinearParameters.min_cosine_crossing_angles)
            {
                crossing = false;
            }

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
                std::cout << " spline intersection: " << splineIndex << "  " << s << "   " << firstSplineRatio << "  " << secondSplineRatio << "  " << crossProductIntersection << std::endl;
            }
        }

        auto const sortedIndices = lin_alg::SortRow(m_crossSplineCoordinates.row(splineIndex));
        lin_alg::ReorderRow(m_crossSplineCoordinates.row(splineIndex), sortedIndices);
        lin_alg::ReorderRow(m_crossingSplinesIndices.row(splineIndex), sortedIndices);
        lin_alg::ReorderRow(m_isLeftOriented.row(splineIndex), sortedIndices);
    }

    void CurvilinearGridFromSplines::MakeAllGridLines()
    {
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

        for (size_t i = 0; i < m_gridLine.size(); ++i)
        {
            std::cout << "before grid line (" << i << ") = " << m_gridLine[i].x << ", " << m_gridLine[i].y << std::endl;
        }

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

        for (size_t i = 0; i < m_gridLine.size(); ++i)
        {
            std::cout << "initial grid line (" << i << ") = " << m_gridLine[i].x << ", " << m_gridLine[i].y << std::endl;
        }

        for (size_t i = 0; i < m_gridLineDimensionalCoordinates.size(); ++i)
        {
            std::cout << "initial grid line coords (" << i << ") = " << m_gridLineDimensionalCoordinates[i] << std::endl;
        }
    }

    UInt CurvilinearGridFromSplines::MakeGridLine(UInt splineIndex,
                                                  UInt startingIndex)
    {

        std::cout << " CurvilinearGridFromSplines::MakeGridLine " << std::endl;

        // first estimation of nodes along m
        auto numM = 1 + static_cast<UInt>(std::floor(m_splines->m_splinesLength[splineIndex] / m_splinesToCurvilinearParameters.average_width));
        numM = std::min(numM, static_cast<UInt>(m_curvilinearParameters.m_refinement));

        std::cout << " numM " << numM << std::endl;

        double parameterStart = 0.278436;
        double parameterEnd = 1.82739;

        [[maybe_unused]] const auto endSplineAdimensionalCoordinate = static_cast<double>(m_splines->m_splineNodes[splineIndex].size()) - 1.0;
        const auto splineStart = m_splines->ComputeSplineLength(splineIndex,
                                                                0.0,
                                                                parameterStart, // 0.0,
                                                                10,
                                                                m_splinesToCurvilinearParameters.curvature_adapted_grid_spacing,
                                                                m_maximumGridHeights[splineIndex]);

        const auto splineLength = m_splines->ComputeSplineLength(splineIndex,
                                                                 parameterStart,
                                                                 parameterEnd,
                                                                 // 0.16666666666666666666,
                                                                 // 1.833333333333333333333,
                                                                 // 0.0,
                                                                 // endSplineAdimensionalCoordinate,
                                                                 10,
                                                                 m_splinesToCurvilinearParameters.curvature_adapted_grid_spacing,
                                                                 m_maximumGridHeights[splineIndex]);

        std::cout << "spline length: " << splineStart << "  " << splineLength << "  "
                  << m_splines->Evaluate(splineIndex, 0.0).x << "  " << m_splines->Evaluate(splineIndex, 0.0).y << " -- "
                  << m_splines->Evaluate(splineIndex, 0.1666666666666666666666666).x << "  " << m_splines->Evaluate(splineIndex, 0.1666666666666666666666666).y << " -- "
                  << std::endl;

        m_gridLine[startingIndex] = m_splines->Evaluate(splineIndex, parameterStart);
        // m_gridLine[startingIndex] = m_splines->Evaluate(splineIndex, 0.1666666666666666666666666);
        // m_gridLine[startingIndex] = m_splines->m_splineNodes[splineIndex][0];

        auto currentMaxWidth = std::numeric_limits<double>::max();
        std::vector<double> distances(numM);
        while (currentMaxWidth > m_splinesToCurvilinearParameters.average_width)
        {
            currentMaxWidth = 0.0;
            std::cout << "distances ";

            for (UInt n = 0; n < numM; ++n)
            {
                distances[n] = splineStart + splineLength * (n + 1.0) / static_cast<double>(numM);
                std::cout << distances[n] << ",  ";
            }
            std::cout << std::endl;

            auto [points, adimensionalDistances] = m_splines->ComputePointOnSplineFromAdimensionalDistance(splineIndex,
                                                                                                           m_maximumGridHeights[splineIndex],
                                                                                                           m_splinesToCurvilinearParameters.curvature_adapted_grid_spacing,
                                                                                                           distances);
            std::cout << "adimensionalDistances ";

            for (size_t i = 0; i < adimensionalDistances.size(); ++i)
            {
                std::cout << adimensionalDistances[i] << ",  ";
            }

            std::cout << std::endl;
            std::cout << "points: ";

            // If distances is initialised in size numM + 1 and with 0 in the zeroth position.
            // can this loop be moved outside the while loop, except the finding of th
            for (UInt n = 0; n < numM; ++n)
            {
                std::cout << " {" << points[n].x << ", " << points[n].y << "} ";
            }

            std::cout << std::endl;
            std::cout << "other points " << currentMaxWidth << "  ";

            for (UInt n = 0; n < numM; ++n)
            {
                const auto index = startingIndex + n + 1;
                m_gridLineDimensionalCoordinates[index] = adimensionalDistances[n];

                m_gridLine[index] = points[n];
                currentMaxWidth = std::max(currentMaxWidth, ComputeDistance(m_gridLine[index - 1], m_gridLine[index], m_splines->m_projection));
                std::cout << " {" << m_gridLine[index - 1].x << ", " << m_gridLine[index - 1].y << "} "
                          << " {" << m_gridLine[index].x << ", " << m_gridLine[index].y << "} "
                          << ComputeDistance(m_gridLine[index - 1], m_gridLine[index], m_splines->m_projection)
                          << " -- ";
            }

            std::cout << std::endl;
            std::cout << " currentMaxWidth " << currentMaxWidth << "  " << m_splinesToCurvilinearParameters.average_width << "  "
                      << numM << "  " << m_curvilinearParameters.m_refinement
                      << std::endl;

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

        std::cout << "points for grid line: " << splineIndex << "  " << startingIndex << std::endl;

        for (UInt i = 0; i < numM; ++i)
        {
            std::cout << "grid line point " << i << " = {" << m_gridLine[i].x << ", " << m_gridLine[i].y << "}" << std::endl;
        }

        return numM;
    }

    void CurvilinearGridFromSplines::ComputeSplineProperties(const bool restoreOriginalProperties)
    {
        std::cout << " CurvilinearGridFromSplines::ComputeSplineProperties " << std::boolalpha << restoreOriginalProperties << std::endl;
        AllocateSplinesProperties();

        for (UInt s = 0; s < m_splines->GetNumSplines(); ++s)
        {
            std::cout << "GetSplineIntersections " << s << std::endl;
            GetSplineIntersections(s);
            std::cout << std::endl;
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
                    std::cout << " middle m_centralSplineIndex[index] " << s << "  " << index << "  " << m_centralSplineIndex[index] << std::endl;
                }
                for (auto i = middleCrossingSpline + 1; i < m_numCrossingSplines[s]; ++i)
                {
                    const auto index = m_crossingSplinesIndices(s, i);
                    m_type[index] = SplineTypes::lateral; // lateral spline
                    m_centralSplineIndex[index] = -static_cast<int>(crossingSplineIndex);
                    std::cout << " outer m_centralSplineIndex[index] " << s << "  " << index << "  " << m_centralSplineIndex[index] << std::endl;
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
                continue;
            }

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
