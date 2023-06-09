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
#include <MeshKernel/CurvilinearGrid/CurvilinearGridCreateUniform.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Polygons.hpp>

using meshkernel::CurvilinearGrid;
using meshkernel::CurvilinearGridCreateUniform;

CurvilinearGridCreateUniform::CurvilinearGridCreateUniform(const MakeGridParameters& makeGridParameters, Projection projection)
    : m_projection(projection),
      m_numM(makeGridParameters.num_columns + 1),
      m_numN(makeGridParameters.num_rows + 1),
      m_angle(makeGridParameters.angle),
      m_XGridBlockSize(makeGridParameters.block_size_x),
      m_YGridBlockSize(makeGridParameters.block_size_y),
      m_OriginXCoordinate(makeGridParameters.origin_x),
      m_OriginYCoordinate(makeGridParameters.origin_y),
      m_block_size_x(makeGridParameters.block_size_x),
      m_block_size_y(makeGridParameters.block_size_y)
{
}

CurvilinearGrid CurvilinearGridCreateUniform::Compute() const
{
    if (m_projection == Projection::spherical)
    {
        return CurvilinearGrid{computeSpherical(m_OriginXCoordinate,
                                                m_OriginYCoordinate,
                                                m_numM,
                                                m_numN,
                                                m_XGridBlockSize,
                                                m_YGridBlockSize,
                                                m_angle),
                               m_projection};
    }
    if (m_projection == Projection::cartesian)
    {
        return CurvilinearGrid{computeCartesian(m_OriginXCoordinate,
                                                m_OriginYCoordinate,
                                                m_numM,
                                                m_numN,
                                                m_XGridBlockSize,
                                                m_YGridBlockSize, m_angle),
                               m_projection};
    }
    return CurvilinearGrid();
}

std::vector<std::vector<meshkernel::Point>> CurvilinearGridCreateUniform::computeCartesian(double OriginXCoordinate,
                                                                                           double OriginYCoordinate,
                                                                                           size_t numM,
                                                                                           size_t numN,
                                                                                           double XGridBlockSize,
                                                                                           double YGridBlockSize,
                                                                                           double angle)

{

    const auto cosineAngle = std::cos(angle * constants::conversion::degToRad);
    const auto sinAngle = std::sin(angle * constants::conversion::degToRad);
    std::vector result(numN, std::vector<Point>(numM));
    for (size_t n = 0; n < numN; ++n)
    {
        for (size_t m = 0; m < numM; ++m)
        {
            const double newPointXCoordinate = OriginXCoordinate + m * XGridBlockSize * cosineAngle - n * YGridBlockSize * sinAngle;
            const double newPointYCoordinate = OriginYCoordinate + m * XGridBlockSize * sinAngle + n * YGridBlockSize * cosineAngle;
            result[n][m] = {newPointXCoordinate, newPointYCoordinate};
        }
    }
    return result;
}

std::vector<std::vector<meshkernel::Point>> CurvilinearGridCreateUniform::computeSpherical(double OriginXCoordinate,
                                                                                           double OriginYCoordinate,
                                                                                           size_t numM,
                                                                                           size_t numN,
                                                                                           double XGridBlockSize,
                                                                                           double YGridBlockSize,
                                                                                           double angle)
{
    std::vector result = computeCartesian(OriginXCoordinate,
                                          OriginYCoordinate,
                                          numM,
                                          numN,
                                          XGridBlockSize,
                                          YGridBlockSize,
                                          angle);

    constexpr double latitudeCloseToPoles = 88.0;
    constexpr double minimumDistance = 2000;
    bool onPoles = false;
    size_t lastRowOnPole;

    for (size_t n = 1; n < numN; ++n)
    {
        for (size_t m = 0; m < numM; ++m)
        {
            // adjust the latitude to preserve an aspect ratio of 1
            const auto c = cos(constants::conversion::degToRad * result[n - 1][m].y);
            const auto asp = c * 1.0 + (1.0 - c) * 0.3;
            const auto dy = XGridBlockSize * c * asp;

            const auto newPointYCoordinate = result[n - 1][m].y + dy;
            result[n][m].y = newPointYCoordinate;

            // prevent too small increments for dy
            const double distance = dy * constants::conversion::degToRad * constants::geometric::earth_radius;
            if (abs(result[n][m].y) > latitudeCloseToPoles && distance < minimumDistance)
            {
                double sign = result[n][m].y < 0 ? -1.0 : 1.0;
                result[n][m].y = 90.0 * sign;
                onPoles = true;
                lastRowOnPole = n;
            }
        }
        if (onPoles)
        {
            break;
        }
    }

    if (onPoles)
    {
        result.erase(result.begin() + lastRowOnPole + 1, result.end());
    }
    return result;
}

CurvilinearGrid CurvilinearGridCreateUniform::Compute(std::shared_ptr<Polygons> polygons, size_t polygonIndex) const
{
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

    const auto cosineAngle = std::cos(m_angle * constants::conversion::degToRad);
    const auto sinAngle = std::sin(m_angle * constants::conversion::degToRad);

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

    auto const OriginXCoordinate = referencePoint.x + xShift;
    auto const OriginYCoordinate = referencePoint.y + yShift;
    auto const numN = static_cast<size_t>(std::ceil((etamax - etamin) / m_block_size_x) + 1);
    auto const numM = static_cast<size_t>(std::ceil((xmax - xmin) / m_block_size_y) + 1);

    CurvilinearGrid curvilinearGrid;
    if (m_projection == Projection::spherical)
    {
        curvilinearGrid = CurvilinearGrid{computeSpherical(OriginXCoordinate,
                                                           OriginYCoordinate,
                                                           numM,
                                                           numN,
                                                           m_block_size_x,
                                                           m_block_size_y,
                                                           m_angle),
                                          m_projection};
    }
    else if (m_projection == Projection::cartesian)
    {
        curvilinearGrid = CurvilinearGrid{computeCartesian(OriginXCoordinate,
                                                           OriginYCoordinate,
                                                           numM,
                                                           numN,
                                                           m_block_size_x,
                                                           m_block_size_y,
                                                           m_angle),
                                          m_projection};
    }

    // remove nodes outside the polygon
    curvilinearGrid.Delete(polygons, polygonIndex);

    return curvilinearGrid;
}
