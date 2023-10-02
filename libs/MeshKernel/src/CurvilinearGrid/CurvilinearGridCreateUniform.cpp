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

#include "MeshKernel/Exceptions.hpp"

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridCreateUniform.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/RangeCheck.hpp>

#include <cmath>

namespace meshkernel
{

    CurvilinearGridCreateUniform::CurvilinearGridCreateUniform(Projection projection) : m_projection(projection)
    {
        if (m_projection != Projection::cartesian && m_projection != Projection::spherical)
        {
            throw meshkernel::NotImplementedError("Projection value: {} not supported", static_cast<int>(m_projection));
        }
    }

    CurvilinearGrid CurvilinearGridCreateUniform::Compute(const int numColumns,
                                                          const int numRows,
                                                          const double originX,
                                                          const double originY,
                                                          const double angle,
                                                          const double blockSizeX,
                                                          const double blockSizeY) const
    {
        range_check::CheckGreater(numColumns, 0, "Number of columns");
        range_check::CheckGreater(numRows, 0, "Number of rows");
        range_check::CheckInOpenInterval(angle, {-90.0, 90.0}, "Grid angle");
        range_check::CheckGreater(blockSizeX, 0.0, "X block size");
        range_check::CheckGreater(blockSizeY, 0.0, "Y block size");

        if (m_projection == Projection::spherical)
        {
            return CurvilinearGrid{ComputeSpherical(numColumns,
                                                    numRows,
                                                    originX,
                                                    originY,
                                                    angle,
                                                    blockSizeX,
                                                    blockSizeY),
                                   m_projection};
        }
        if (m_projection == Projection::cartesian)
        {
            return CurvilinearGrid{ComputeCartesian(numColumns,
                                                    numRows,
                                                    originX,
                                                    originY,
                                                    angle,
                                                    blockSizeX,
                                                    blockSizeY),
                                   m_projection};
        }
        throw NotImplementedError("Projection value {} not supported", static_cast<int>(m_projection));
    }

    lin_alg::Matrix<Point> CurvilinearGridCreateUniform::ComputeCartesian(const int numColumns,
                                                                          const int numRows,
                                                                          const double originX,
                                                                          const double originY,
                                                                          const double angle,
                                                                          const double blockSizeX,
                                                                          const double blockSizeY)

    {
        const auto angleInRad = angle * constants::conversion::degToRad;
        const auto cosineAngle = std::cos(angleInRad);
        const auto sinAngle = std::sin(angleInRad);

        const auto numM = numColumns + 1;
        const auto numN = numRows + 1;

        lin_alg::Matrix<Point> result(numN, numM);
        const auto blockSizeXByCos = blockSizeX * cosineAngle;
        const auto blockSizeYbySin = blockSizeY * sinAngle;
        const auto blockSizeXBySin = blockSizeX * sinAngle;
        const auto blockSizeYByCos = blockSizeY * cosineAngle;
        for (Eigen::Index n = 0; n < result.rows(); ++n)
        {
            for (Eigen::Index m = 0; m < result.cols(); ++m)
            {
                const double newPointXCoordinate = originX + m * blockSizeXByCos - n * blockSizeYbySin;
                const double newPointYCoordinate = originY + m * blockSizeXBySin + n * blockSizeYByCos;
                result(n, m) = {newPointXCoordinate, newPointYCoordinate};
            }
        }
        return result;
    }

    lin_alg::Matrix<Point> CurvilinearGridCreateUniform::ComputeSpherical(const int numColumns,
                                                                          const int numRows,
                                                                          const double originX,
                                                                          const double originY,
                                                                          const double angle,
                                                                          const double blockSizeX,
                                                                          const double blockSizeY)
    {
        lin_alg::Matrix<Point> result = ComputeCartesian(numColumns,
                                                         numRows,
                                                         originX,
                                                         originY,
                                                         angle,
                                                         blockSizeX,
                                                         blockSizeY);

        const auto numM = result.cols();
        const auto numN = result.rows();
        bool onPoles = false;
        constexpr double latitudePoles = 90.0;

        for (Eigen::Index n = 1; n < numN; ++n)
        {
            Eigen::Index lastRowOnPole = numM;
            for (Eigen::Index m = 0; m < numM; ++m)
            {
                const double adjustedLatitude = ComputeLatitudeIncrementWithAdjustment(blockSizeY, result(n - 1, m).y);
                result(n, m).y = adjustedLatitude;

                if (IsEqual(std::abs(adjustedLatitude), latitudePoles))
                {
                    onPoles = true;
                    lastRowOnPole = n;
                }
            }
            if (onPoles)
            {
                if (lastRowOnPole + 1 < result.rows())
                {
                    lin_alg::EraseRows(result, lastRowOnPole + 1, result.rows() - 1);
                }
                break;
            }
        }
        return result;
    }

    double CurvilinearGridCreateUniform::ComputeLatitudeIncrementWithAdjustment(double blockSize, double latitude)
    {
        constexpr double latitudeCloseToPoles = 88.0; // The latitude defining close to poles
        constexpr double minimumDistance = 2000;      // When the real distance along the latitude becomes smaller than minimumDistance and the location is close to the poles, snap the next point to the poles.

        const auto latitudeInRadiants = std::cos(constants::conversion::degToRad * latitude);
        const auto asp = latitudeInRadiants + (1.0 - latitudeInRadiants) * 0.3;
        const auto dy = blockSize * latitudeInRadiants * asp;

        double result = latitude + dy;
        // prevent too small dy increments in case we are on the poles
        const double difference = dy * constants::conversion::degToRad * constants::geometric::earth_radius;
        if (abs(result) > latitudeCloseToPoles && difference < minimumDistance)
        {
            const double sign = result < 0 ? -1.0 : 1.0;
            result = 90.0 * sign;
        }

        return result;
    }

    int CurvilinearGridCreateUniform::ComputeNumRows(double minY,
                                                     double maxY,
                                                     double blockSizeY,
                                                     Projection projection)
    {
        if (blockSizeY > std::abs(maxY - minY))
        {
            throw AlgorithmError("blockSizeY cannot be larger than mesh height");
        }

        if (projection == Projection::cartesian)
        {
            const int numM = static_cast<int>(std::ceil(std::abs(maxY - minY) / blockSizeY));
            return std::max(numM, 1);
        }

        double currentLatitude = minY;
        int result = 0;
        constexpr double latitudePoles = 90.0;

        while (currentLatitude < maxY)
        {
            currentLatitude = ComputeLatitudeIncrementWithAdjustment(blockSizeY, currentLatitude);
            result += 1;

            if (IsEqual(abs(currentLatitude), latitudePoles))
            {
                break;
            }
        }

        return result;
    }

    CurvilinearGrid CurvilinearGridCreateUniform::Compute(const double angle,
                                                          const double blockSizeX,
                                                          const double blockSizeY,
                                                          std::shared_ptr<Polygons> polygons,
                                                          UInt polygonIndex) const
    {

        range_check::CheckInOpenInterval(angle, {-90.0, 90.0}, "Grid angle");
        range_check::CheckGreater(blockSizeX, 0.0, "X block size");
        range_check::CheckGreater(blockSizeY, 0.0, "Y block size");

        if (polygons->IsEmpty())
        {
            throw AlgorithmError("Enclosures list is empty.");
        }

        if (polygons->GetProjection() != m_projection)
        {
            throw AlgorithmError("Polygon projection ({}) is not equal to curvilinear grid projection ({})",
                                 ToString(polygons->GetProjection()),
                                 ToString(m_projection));
        }

        // Compute the bounding box
        const auto boundingBox = polygons->GetBoundingBox(polygonIndex);
        const auto referencePoint = boundingBox.MassCentre();

        // Compute the max size
        const auto maxSize = std::max(boundingBox.Width(), boundingBox.Height());

        // Compute the lower left and upper right corners
        const Point lowerLeft(referencePoint.x - maxSize, referencePoint.y - maxSize);
        const Point upperRight(referencePoint.x + maxSize, referencePoint.y + maxSize);

        // Compute the number of rows and columns
        const int numColumns = std::max(static_cast<int>(std::ceil(std::abs(upperRight.x - lowerLeft.x) / blockSizeX)), 1);
        const int numRows = ComputeNumRows(lowerLeft.y, upperRight.y, blockSizeY, m_projection);

        // Rotated the lower left corner
        const auto lowerLeftMergedRotated = Rotate(lowerLeft, angle, referencePoint);

        // Set the origin
        const double originX = lowerLeftMergedRotated.x;
        const double originY = lowerLeftMergedRotated.y;

        CurvilinearGrid curvilinearGrid;
        switch (m_projection)
        {
        case Projection::spherical:
            curvilinearGrid = CurvilinearGrid{ComputeSpherical(numColumns,
                                                               numRows,
                                                               originX,
                                                               originY,
                                                               angle,
                                                               blockSizeX,
                                                               blockSizeY),
                                              m_projection};
            break;
        case Projection::cartesian:
            curvilinearGrid = CurvilinearGrid{ComputeCartesian(numColumns,
                                                               numRows,
                                                               originX,
                                                               originY,
                                                               angle,
                                                               blockSizeX,
                                                               blockSizeY),
                                              m_projection};
            break;
        default:
            throw NotImplementedError("Projection value {} not supported", static_cast<int>(m_projection));
        }

        // remove nodes outside the polygon
        curvilinearGrid.Delete(polygons, polygonIndex);

        return curvilinearGrid;
    }

    CurvilinearGrid CurvilinearGridCreateUniform::Compute(const double originX,
                                                          const double originY,
                                                          const double blockSizeX,
                                                          const double blockSizeY,
                                                          const double upperRightX,
                                                          const double upperRightY) const
    {
        range_check::CheckGreater(blockSizeX, 0.0, "X block size");
        range_check::CheckGreater(blockSizeY, 0.0, "Y block size");

        const int numColumns = static_cast<int>(std::ceil((upperRightX - originX) / blockSizeX));
        if (numColumns <= 0)
        {
            throw AlgorithmError("Number of columns cannot be <= 0");
        }

        const double angle = 0.0;
        const int numRows = ComputeNumRows(originY, upperRightY, blockSizeY, m_projection);

        CurvilinearGrid curvilinearGrid;
        switch (m_projection)
        {
        case Projection::spherical:

            curvilinearGrid = CurvilinearGrid{ComputeSpherical(numColumns,
                                                               numRows,
                                                               originX,
                                                               originY,
                                                               angle,
                                                               blockSizeX,
                                                               blockSizeY),
                                              m_projection};
            break;
        case Projection::cartesian:
            curvilinearGrid = CurvilinearGrid{ComputeCartesian(numColumns,
                                                               numRows,
                                                               originX,
                                                               originY,
                                                               angle,
                                                               blockSizeX,
                                                               blockSizeY),
                                              m_projection};
            break;
        case Projection::sphericalAccurate:
        default:
            throw meshkernel::NotImplementedError("Projection value {} not supported", static_cast<int>(m_projection));
        }

        return curvilinearGrid;
    }

} // namespace meshkernel