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

CurvilinearGridCreateUniform::CurvilinearGridCreateUniform(const meshkernelapi::MakeGridParameters& makeGridParameters, Projection projection)
    : m_makeGridParameters(makeGridParameters), m_projection(projection)
{
}

CurvilinearGrid CurvilinearGridCreateUniform::Compute() const
{
    // regular grid
    auto numM = static_cast<size_t>(m_makeGridParameters.num_columns + 1);
    auto numN = static_cast<size_t>(m_makeGridParameters.num_rows + 1);
    const double XGridBlockSize = m_makeGridParameters.block_size_x;
    const double YGridBlockSize = m_makeGridParameters.block_size_y;
    const double cosineAngle = std::cos(m_makeGridParameters.angle * degrad_hp);
    const double sinAngle = std::sin(m_makeGridParameters.angle * degrad_hp);
    double OriginXCoordinate = m_makeGridParameters.origin_x;
    double OriginYCoordinate = m_makeGridParameters.origin_y;

    // No polygon, use MakeGridParameters as is

    std::vector<std::vector<Point>> gridNodes(numN, std::vector<Point>(numM));
    for (size_t n = 0; n < numN; ++n)
    {
        for (size_t m = 0; m < numM; ++m)
        {
            const double newPointXCoordinate = OriginXCoordinate + m * XGridBlockSize * cosineAngle - n * YGridBlockSize * sinAngle;
            double newPointYCoordinate = OriginYCoordinate + m * XGridBlockSize * sinAngle + n * YGridBlockSize * cosineAngle;
            if (m_projection == Projection::spherical && n > 0)
            {
                newPointYCoordinate = XGridBlockSize * cos(degrad_hp * gridNodes[n - 1][m].y);
            }
            gridNodes[n][m] = {newPointXCoordinate, newPointYCoordinate};
        }
    }
    return CurvilinearGrid{gridNodes, m_projection};
}

CurvilinearGrid CurvilinearGridCreateUniform::Compute(std::shared_ptr<Polygons> polygons, size_t polygonIndex) const
{
    if (polygons->m_projection != m_projection)
    {
        throw std::invalid_argument("CurvilinearGridCreateUniform::Compute polygon projection is not equal to CurvilinearGridCreateUniform projection ");
    }

    if (polygons->IsEmpty())
    {
        return {};
    }

    Point referencePoint{doubleMissingValue, doubleMissingValue};
    auto const& [startPolygonIndex, endPolygonIndex] = polygons->m_outer_polygons_indices[polygonIndex];
    for (auto i = startPolygonIndex; i <= endPolygonIndex; ++i)
    {
        if (polygons->m_nodes[i].IsValid())
        {
            referencePoint = polygons->m_nodes[i];
            break;
        }
    }

    // get polygon min/max in rotated (xi,eta) coordinates
    double xmin = std::numeric_limits<double>::max();
    double xmax = -xmin;
    double etamin = std::numeric_limits<double>::max();
    double etamax = -etamin;

    const double XGridBlockSize = m_makeGridParameters.block_size_x;
    const double YGridBlockSize = m_makeGridParameters.block_size_y;
    const double cosineAngle = std::cos(m_makeGridParameters.angle * degrad_hp);
    const double sinAngle = std::sin(m_makeGridParameters.angle * degrad_hp);

    for (auto i = startPolygonIndex; i <= endPolygonIndex; ++i)
    {
        if (polygons->m_nodes[i].IsValid())
        {
            const double dx = GetDx(referencePoint, polygons->m_nodes[i], polygons->m_projection);
            const double dy = GetDy(referencePoint, polygons->m_nodes[i], polygons->m_projection);
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
    if (polygons->m_projection == Projection::spherical)
    {
        xShift = xShift / earth_radius * raddeg_hp;
        yShift = yShift / (earth_radius * std::cos(referencePoint.y * degrad_hp)) * raddeg_hp;
    }

    auto const OriginXCoordinate = referencePoint.x + xShift;
    auto const OriginYCoordinate = referencePoint.y + yShift;
    auto const numN = static_cast<size_t>(std::ceil((etamax - etamin) / XGridBlockSize) + 1);
    auto const numM = static_cast<size_t>(std::ceil((xmax - xmin) / YGridBlockSize) + 1);

    std::vector<std::vector<Point>> gridNodes(numN, std::vector<Point>(numM));
    for (size_t n = 0; n < numN; ++n)
    {
        for (size_t m = 0; m < numM; ++m)
        {
            const double newPointXCoordinate = OriginXCoordinate + m * XGridBlockSize * cosineAngle - n * YGridBlockSize * sinAngle;
            double newPointYCoordinate = OriginYCoordinate + m * XGridBlockSize * sinAngle + n * YGridBlockSize * cosineAngle;
            if (polygons->m_projection == Projection::spherical && n > 0)
            {
                newPointYCoordinate = XGridBlockSize * cos(degrad_hp * gridNodes[n - 1][m].y);
            }
            gridNodes[n][m] = {newPointXCoordinate, newPointYCoordinate};
        }
    }

    CurvilinearGrid curvilinearGrid(gridNodes, m_projection);

    // remove nodes outside the polygon
    curvilinearGrid.Delete(polygons, polygonIndex);

    return curvilinearGrid;
}
