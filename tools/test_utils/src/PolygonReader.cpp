//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2025.
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

#include "TestUtils/PolygonReader.hpp"
#include <fstream>
#include <iostream>
#include <vector>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Point.hpp"

size_t readPolygonSize(const std::string& line)
{
    std::stringstream stream(line);
    size_t size;

    stream >> size;
    return size;
}

meshkernel::Point readPoint(const std::string& line)
{
    std::stringstream stream(line);

    meshkernel::Point point;

    stream >> point.x;
    stream >> point.y;
    return point;
}

std::unique_ptr<meshkernel::Polygons> ReadPolygons(const std::string& fileName, const meshkernel::Projection projection)
{
    std::unique_ptr<meshkernel::Polygons> polygons;

    std::string line;
    std::ifstream infile(fileName.c_str());
    std::vector<meshkernel::Point> points;

    size_t polygonCount = 0;

    while (std::getline(infile, line))
    {
        if (!line.empty())
        {
            if (line[0] == 'L')
            {
                ++polygonCount;

                if (polygonCount > 1)
                {
                    points.push_back(meshkernel::Point(meshkernel::constants::missing::doubleValue, meshkernel::constants::missing::doubleValue));
                }

                std::getline(infile, line);
                size_t size = readPolygonSize(line);

                size_t currentSize = points.size();

                for (size_t i = 0; i < size; ++i)
                {
                    std::getline(infile, line);
                    points.push_back(readPoint(line));
                }

                points.push_back(points[currentSize]);
            }
        }
    }

    if (points.size() > 0)
    {
        polygons = std::make_unique<meshkernel::Polygons>(points, projection);
    }

    return polygons;
}
