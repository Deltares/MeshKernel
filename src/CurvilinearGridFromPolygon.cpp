//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2020.
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

#include <vector>
#include <iostream>
#include <algorithm>
#include "Operations.cpp"
#include "Entities.hpp"
#include "Mesh.hpp"
#include "Polygons.hpp"
#include "LandBoundaries.hpp"
#include "CurvilinearGridFromPolygon.hpp"

GridGeom::CurvilinearGridFromPolygon::CurvilinearGridFromPolygon() :
    m_mesh(nullptr),
    m_polygon(nullptr)
{
}

GridGeom::CurvilinearGridFromPolygon::CurvilinearGridFromPolygon(std::shared_ptr<Mesh> mesh, 
                                                                 std::shared_ptr<Polygons> polygon) :
    m_mesh(mesh),
    m_polygon(polygon)
{

};

bool GridGeom::CurvilinearGridFromPolygon::Compute()
{
    bool successful = true;
    return true;
}