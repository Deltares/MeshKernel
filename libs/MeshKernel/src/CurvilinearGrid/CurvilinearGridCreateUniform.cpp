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
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Polygons.hpp>

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridCreateUniform;

CurvilinearGridCreateUniform::CurvilinearGridCreateUniform(Projection projection) : m_projection(projection)
{
}

CurvilinearGrid CurvilinearGridCreateUniform::Compute(const int numColumns,
                                                      const int numRows,
                                                      const double originX,
                                                      const double originY,
                                                      const double angle,
                                                      const double blockSizeX,
                                                      const double blockSizeY) const
{
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

    const std::string message = "Projection value: " + std::to_string(static_cast<int>(m_projection)) + " not supported";
    throw NotImplemented(message);
}

std::vector<std::vector<meshkernel::Point>> CurvilinearGridCreateUniform::ComputeCartesian(const int numColumns,
                                                                                           const int numRows,
                                                                                           const double originX,
                                                                                           const double originY,
                                                                                           const double angle,
                                                                                           const double blockSizeX,
                                                                                           const double blockSizeY)

{
    if (numColumns <= 0)
    {
        throw AlgorithmError("Number of columns cannot be <= 0");
    }
    if (numRows <= 0)
    {
        throw AlgorithmError("Number of rows cannot be <= 0");
    }
    if (blockSizeX <= 0.0)
    {
        throw AlgorithmError("BlockSizeX cannot be <= 0");
    }
    if (blockSizeY <= 0.0)
    {
        throw AlgorithmError("BlockSizeY cannot be <= 0");
    }

    const auto angleInRad = angle * constants::conversion::degToRad;
    const auto cosineAngle = std::cos(angleInRad);
    const auto sinAngle = std::sin(angleInRad);

    const auto numM = numColumns + 1;
    const auto numN = numRows + 1;

    std::vector<std::vector<Point>> result(numN, std::vector<Point>(numM));
    const auto blockSizeXByCos = blockSizeX * cosineAngle;
    const auto blockSizeYbySin = blockSizeY * sinAngle;
    const auto blockSizeXBySin = blockSizeX * sinAngle;
    const auto blockSizeYByCos = blockSizeY * cosineAngle;
    for (int n = 0; n < numN; ++n)
    {
        for (int m = 0; m < numM; ++m)
        {
            const double newPointXCoordinate = originX + m * blockSizeXByCos - n * blockSizeYbySin;
            const double newPointYCoordinate = originY + m * blockSizeXBySin + n * blockSizeYByCos;
            result[n][m] = {newPointXCoordinate, newPointYCoordinate};
        }
    }
    return result;
}

std::vector<std::vector<meshkernel::Point>> CurvilinearGridCreateUniform::ComputeSpherical(const int numColumns,
                                                                                           const int numRows,
                                                                                           const double originX,
                                                                                           const double originY,
                                                                                           const double angle,
                                                                                           const double blockSizeX,
                                                                                           const double blockSizeY)
{
    if (numColumns <= 0)
    {
        throw AlgorithmError("Number of columns cannot be <= 0");
    }
    if (numRows <= 0)
    {
        throw AlgorithmError("Number of rows cannot be <= 0");
    }
    if (blockSizeX <= 0.0)
    {
        throw AlgorithmError("BlockSizeX cannot be <= 0");
    }
    if (blockSizeY <= 0.0)
    {
        throw AlgorithmError("BlockSizeY cannot be <= 0");
    }

    std::vector result = ComputeCartesian(numColumns,
                                          numRows,
                                          originX,
                                          originY,
                                          angle,
                                          blockSizeX,
                                          blockSizeY);

    const auto numM = result[0].size();
    const auto numN = result.size();
    bool onPoles = false;
    constexpr double latitudePoles = 90.0;

    for (size_t n = 1; n < numN; ++n)
    {
        size_t lastRowOnPole = numM;
        for (size_t m = 0; m < numM; ++m)
        {
            const double adjustedLatitude = ComputeLatitudeIncrementWithAdjustment(blockSizeY, result[n - 1][m].y);
            result[n][m].y = adjustedLatitude;

            if (IsEqual(abs(adjustedLatitude), latitudePoles))
            {
                onPoles = true;
                lastRowOnPole = n;
            }
        }
        if (onPoles)
        {
            result.erase(result.begin() + lastRowOnPole + 1, result.end());
            break;
        }
    }
    return result;
}

double CurvilinearGridCreateUniform::ComputeLatitudeIncrementWithAdjustment(double blockSize, double latitude)
{
    constexpr double latitudeCloseToPoles = 88.0; // The latitude defining close to poles
    constexpr double minimumDistance = 2000;      // When the real distance along the latitude becomes smaller than minimumDistance and the location is close to the poles, snap the next point to the poles.

    const auto latitudeInRadiants = cos(constants::conversion::degToRad * latitude);
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

int CurvilinearGridCreateUniform::ComputeNumRowsSpherical(double minY,
                                                          double maxY,
                                                          double blockSizeY)
{
    if (blockSizeY <= 0.0)
    {
        throw AlgorithmError("blockSizeY cannot be <= 0");
    }
    if (blockSizeY > std::abs(maxY - minY))
    {
        throw AlgorithmError("blockSizeY cannot be larger than mesh height");
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
    if (blockSizeX <= 0.0)
    {
        throw AlgorithmError("BlockSizeX cannot be <= 0");
    }
    if (blockSizeY <= 0.0)
    {
        throw AlgorithmError("BlockSizeY cannot be <= 0");
    }

    if (polygons->GetProjection() != m_projection)
    {
        throw std::invalid_argument("CurvilinearGridCreateUniform::Compute polygon projection is not equal to CurvilinearGridCreateUniform projection ");
    }

    if (polygons->IsEmpty())
    {
        return {};
    }

    Point referencePoint;
    auto const& [startPolygonIndex, endPolygonIndex] = polygons->OuterIndices(polygonIndex);
    for (auto i = startPolygonIndex; i <= endPolygonIndex; ++i)
    {
        auto const& polygonNode = polygons->Node(i);
        if (polygonNode.IsValid())
        {
            referencePoint = polygonNode;
            break;
        }
    }

    // get polygon min/max in rotated (xi,eta) coordinates
    double xmin = std::numeric_limits<double>::max();
    double xmax = -xmin;
    double etamin = std::numeric_limits<double>::max();
    double etamax = -etamin;

    const auto angleInRad = angle * constants::conversion::degToRad;
    const auto cosineAngle = std::cos(angleInRad);
    const auto sinAngle = std::sin(angleInRad);
    Projection const polygonProjection = polygons->GetProjection();

    for (auto i = startPolygonIndex; i <= endPolygonIndex; ++i)
    {
        auto const& polygonNode = polygons->Node(i);
        if (polygonNode.IsValid())
        {
            const double dx = GetDx(referencePoint, polygonNode, polygonProjection);
            const double dy = GetDy(referencePoint, polygonNode, polygonProjection);
            double xi = dx * cosineAngle + dy * sinAngle;
            double eta = -dx * sinAngle + dy * cosineAngle;
            xmin = std::min(xmin, xi);
            xmax = std::max(xmax, xi);
            etamin = std::min(etamin, eta);
            etamax = std::max(etamax, eta);
        }
    }

    double xShift = xmin * cosineAngle - etamin * sinAngle;
    double yShift = xmin * sinAngle + etamin * cosineAngle;
    if (polygonProjection == Projection::spherical)
    {
        xShift = xShift / constants::geometric::earth_radius * constants::conversion::radToDeg;
        yShift = yShift / (constants::geometric::earth_radius * std::cos(referencePoint.y * constants::conversion::degToRad)) * constants::conversion::radToDeg;
    }

    const double originX = referencePoint.x + xShift;
    const double originY = referencePoint.y + yShift;
    const int numM = static_cast<int>(std::ceil((etamax - etamin) / blockSizeX) + 1);
    const int numColumns = numM > 0 ? numM - 1 : 0;
    const int numN = static_cast<int>(std::ceil((xmax - xmin) / blockSizeY) + 1);
    const int numRows = numN > 0 ? numN - 1 : 0;

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
        const std::string message = "Projection value: " + std::to_string(static_cast<int>(m_projection)) + " not supported";
        throw NotImplemented(message);
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
    if (blockSizeX <= 0.0)
    {
        throw AlgorithmError("BlockSizeX cannot be <= 0");
    }
    if (blockSizeY <= 0.0)
    {
        throw AlgorithmError("BlockSizeY cannot be <= 0");
    }

    const int numColumns = static_cast<int>(std::ceil((upperRightX - originX) / blockSizeX));
    if (numColumns <= 0)
    {
        throw AlgorithmError("Number of columns cannot be <= 0");
    }

    const double angle = 0.0;
    int numRows;
    CurvilinearGrid curvilinearGrid;
    switch (m_projection)
    {
    case Projection::spherical:

        numRows = ComputeNumRowsSpherical(originY, upperRightY, blockSizeY);
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
        numRows = static_cast<int>(std::ceil((upperRightY - originY) / blockSizeY));
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
        const std::string message = "Projection value: " + std::to_string(static_cast<int>(m_projection)) + " not supported";
        throw NotImplemented(message);
    }

    return curvilinearGrid;
}
