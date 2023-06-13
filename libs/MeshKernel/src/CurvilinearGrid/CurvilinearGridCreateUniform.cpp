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

CurvilinearGrid CurvilinearGridCreateUniform::Compute(const MakeGridParameters& makeGridParameters) const
{
    if (m_projection == Projection::spherical)
    {
        return CurvilinearGrid{ComputeSpherical(makeGridParameters), m_projection};
    }
    if (m_projection == Projection::cartesian)
    {
        return CurvilinearGrid{ComputeCartesian(makeGridParameters), m_projection};
    }

    const std::string message = "Projection value: " + std::to_string(static_cast<int>(m_projection)) + " not supported";
    throw NotImplemented(message);
}

std::vector<std::vector<meshkernel::Point>> CurvilinearGridCreateUniform::ComputeCartesian(const MakeGridParameters& makeGridParameters)

{
    const auto cosineAngle = std::cos(makeGridParameters.angle * constants::conversion::degToRad);
    const auto sinAngle = std::sin(makeGridParameters.angle * constants::conversion::degToRad);
    const auto numM = makeGridParameters.num_columns + 1;
    const auto numN = makeGridParameters.num_rows + 1;

    std::vector result(numN, std::vector<Point>(numM));
    for (int n = 0; n < numN; ++n)
    {
        for (int m = 0; m < numM; ++m)
        {
            const double newPointXCoordinate = makeGridParameters.origin_x + m * makeGridParameters.block_size_x * cosineAngle - n * makeGridParameters.block_size_y * sinAngle;
            const double newPointYCoordinate = makeGridParameters.origin_y + m * makeGridParameters.block_size_x * sinAngle + n * makeGridParameters.block_size_y * cosineAngle;
            result[n][m] = {newPointXCoordinate, newPointYCoordinate};
        }
    }
    return result;
}

std::vector<std::vector<meshkernel::Point>> CurvilinearGridCreateUniform::ComputeSpherical(const MakeGridParameters& makeGridParameters)
{
    std::vector result = ComputeCartesian(makeGridParameters);
    if (result.empty())
    {
        return result;
    }
    const auto numM = result[0].size();
    const auto numN = result.size();
    constexpr double latitudeCloseToPoles = 88.0;
    constexpr double minimumDistance = 2000;
    bool onPoles = false;

    for (size_t n = 1; n < numN; ++n)
    {
        size_t lastRowOnPole = numM;
        for (size_t m = 0; m < numM; ++m)
        {
            // adjust the latitude to preserve an aspect ratio of 1
            const auto c = cos(constants::conversion::degToRad * result[n - 1][m].y);
            const auto asp = c * 1.0 + (1.0 - c) * 0.3;
            const auto dy = makeGridParameters.block_size_x * c * asp;

            const auto newPointYCoordinate = result[n - 1][m].y + dy;
            result[n][m].y = newPointYCoordinate;

            // prevent too small dy increments in case we are on the poles
            const double distance = dy * constants::conversion::degToRad * constants::geometric::earth_radius;
            if (abs(result[n][m].y) > latitudeCloseToPoles && distance < minimumDistance)
            {
                const double sign = result[n][m].y < 0 ? -1.0 : 1.0;
                result[n][m].y = 90.0 * sign;
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

CurvilinearGrid CurvilinearGridCreateUniform::Compute(const MakeGridParameters& makeGridParameters,
                                                      std::shared_ptr<Polygons> polygons,
                                                      size_t polygonIndex) const
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

    const auto cosineAngle = std::cos(makeGridParameters.angle * constants::conversion::degToRad);
    const auto sinAngle = std::sin(makeGridParameters.angle * constants::conversion::degToRad);

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

    MakeGridParameters makeGridParametersInPolygon(makeGridParameters);

    makeGridParametersInPolygon.origin_x = referencePoint.x + xShift;
    makeGridParametersInPolygon.origin_y = referencePoint.y + yShift;
    const int numM = static_cast<int>(std::ceil((etamax - etamin) / makeGridParameters.block_size_x) + 1);
    makeGridParametersInPolygon.num_columns = numM > 0 ? numM - 1 : 0;
    const int numN = static_cast<int>(std::ceil((xmax - xmin) / makeGridParameters.block_size_y) + 1);
    makeGridParametersInPolygon.num_rows = numN > 0 ? numN - 1 : 0;

    CurvilinearGrid curvilinearGrid;
    if (m_projection == Projection::spherical)
    {
        curvilinearGrid = CurvilinearGrid{ComputeSpherical(makeGridParametersInPolygon),
                                          m_projection};
    }
    else if (m_projection == Projection::cartesian)
    {
        curvilinearGrid = CurvilinearGrid{ComputeCartesian(makeGridParametersInPolygon),
                                          m_projection};
    }

    // remove nodes outside the polygon
    curvilinearGrid.Delete(polygons, polygonIndex);

    return curvilinearGrid;
}
