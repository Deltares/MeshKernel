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

#pragma once

#include "MeshKernel/CurvilinearGridCreateUniform.hpp"

#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGridRefinement.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Polygons.hpp>

meshkernel::CurvilinearGridCreateUniform::CurvilinearGridCreateUniform(const meshkernelapi::MakeMeshParameters& makeMeshParameters, std::shared_ptr<Polygons> polygons)
    : m_makeMeshParameters(makeMeshParameters),
      m_polygons(polygons)
{
}

meshkernel::CurvilinearGrid meshkernel::CurvilinearGridCreateUniform::Compute() const
{
    if (m_makeMeshParameters.GridType != 0)
    {
        throw std::invalid_argument("CurvilinearGridCreateUniform::Compute m_makeMeshParameters.GridType can only be 0");
    }

    // regular grid
    auto numM = static_cast<size_t>(m_makeMeshParameters.NumberOfColumns + 1);
    auto numN = static_cast<size_t>(m_makeMeshParameters.NumberOfRows + 1);
    const double XGridBlockSize = m_makeMeshParameters.XGridBlockSize;
    const double YGridBlockSize = m_makeMeshParameters.YGridBlockSize;
    const double cosineAngle = std::cos(m_makeMeshParameters.GridAngle * degrad_hp);
    const double sinAngle = std::sin(m_makeMeshParameters.GridAngle * degrad_hp);
    double OriginXCoordinate = m_makeMeshParameters.OriginXCoordinate;
    double OriginYCoordinate = m_makeMeshParameters.OriginYCoordinate;

    // in case a polygon is there, re-compute parameters
    if (!m_polygons->IsEmpty())
    {
        Point referencePoint{doubleMissingValue, doubleMissingValue};
        // rectangular grid in polygon

        for (const auto& node : m_polygons->m_nodes)
        {
            if (node.IsValid())
            {
                referencePoint = node;
                break;
            }
        }

        // get polygon min/max in rotated (xi,eta) coordinates
        double xmin = std::numeric_limits<double>::max();
        double xmax = -xmin;
        double etamin = std::numeric_limits<double>::max();
        double etamax = -etamin;
        for (const auto& node : m_polygons->m_nodes)
        {
            if (node.IsValid())
            {
                const double dx = GetDx(referencePoint, node, m_polygons->m_projection);
                const double dy = GetDy(referencePoint, node, m_polygons->m_projection);
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
        if (m_polygons->m_projection == Projection::spherical)
        {
            xShift = xShift / earth_radius * raddeg_hp;
            yShift = yShift / (earth_radius * std::cos(referencePoint.y * degrad_hp)) * raddeg_hp;
        }

        OriginXCoordinate = referencePoint.x + xShift;
        OriginYCoordinate = referencePoint.y + yShift;
        numN = static_cast<size_t>(std::ceil((etamax - etamin) / XGridBlockSize) + 1);
        numM = static_cast<size_t>(std::ceil((xmax - xmin) / YGridBlockSize) + 1);
    }

    std::vector<std::vector<Point>> gridNodes(numN, std::vector<Point>(numM));
    for (auto n = 0; n < numN; ++n)
    {
        for (auto m = 0; m < numM; ++m)
        {
            const double newPointXCoordinate = OriginXCoordinate + m * XGridBlockSize * cosineAngle - n * YGridBlockSize * sinAngle;
            double newPointYCoordinate = OriginYCoordinate + m * XGridBlockSize * sinAngle + n * YGridBlockSize * cosineAngle;
            if (m_polygons->m_projection == Projection::spherical && n > 0)
            {
                newPointYCoordinate = XGridBlockSize * cos(degrad_hp * gridNodes[n - 1][m].y);
            }
            gridNodes[n][m] = {newPointXCoordinate, newPointYCoordinate};
        }
    }

    // in case a polygon is there, remove nodes outside
    if (!m_polygons->IsEmpty())
    {
        std::vector<std::vector<bool>> nodeBasedMask(numN, std::vector<bool>(numM, false));
        std::vector<std::vector<bool>> faceBasedMask(numN - 1, std::vector<bool>(numM - 1, false));
        // mark points inside a polygon
        for (auto n = 0; n < numN; ++n)
        {
            for (auto m = 0; m < numM; ++m)
            {
                const bool isInPolygon = m_polygons->IsPointInPolygon(gridNodes[n][m], 0);
                if (isInPolygon)
                {
                    nodeBasedMask[n][m] = true;
                }
            }
        }

        // mark faces when at least one node is inside
        for (auto n = 0; n < numN - 1; ++n)
        {
            for (auto m = 0; m < numM - 1; ++m)
            {
                if (nodeBasedMask[n][m] || nodeBasedMask[n + 1][m] || nodeBasedMask[n][m + 1] || nodeBasedMask[n + 1][m + 1])
                {
                    faceBasedMask[n][m] = true;
                }
            }
        }

        //mark nodes that are member of a cell inside the polygon(s)
        for (auto n = 0; n < numN - 1; ++n)
        {
            for (auto m = 0; m < numM - 1; ++m)
            {
                if (faceBasedMask[n][m])
                {
                    nodeBasedMask[n][m] = true;
                    nodeBasedMask[n + 1][m] = true;
                    nodeBasedMask[n][m + 1] = true;
                    nodeBasedMask[n + 1][m + 1] = true;
                }
            }
        }

        // mark points inside a polygon
        for (auto n = 0; n < numN; ++n)
        {
            for (auto m = 0; m < numM; ++m)
            {
                if (!nodeBasedMask[n][m])
                {
                    gridNodes[n][m].x = doubleMissingValue;
                    gridNodes[n][m].y = doubleMissingValue;
                }
            }
        }
    }
    return CurvilinearGrid(std::move(gridNodes), m_polygons->m_projection);
}
