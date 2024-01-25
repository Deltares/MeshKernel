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

#include "MeshKernel/Entities.hpp"

#include <fstream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

std::vector<meshkernel::Sample> ReadSampleFile(std::string const& filePath);

template <typename values_type>
std::tuple<int, int, double, double, double, double, std::vector<values_type>> ReadAscFile(const std::string& filePath)
{
    // read sample file
    std::string line;
    std::ifstream infile(filePath.c_str());

    int numX = 0;
    int numY = 0;
    double xllCenter = 0.0;
    double yllCenter = 0.0;
    double cellSize = 0.0;
    double nodataValue;

    int numlines = 0;
    std::vector<std::vector<values_type>> rows;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        if (numlines == 0)
        {
            std::string s;
            iss >> s >> numX;
            numlines++;
            continue;
        }
        if (numlines == 1)
        {
            std::string s;
            iss >> s >> numY;
            numlines++;
            continue;
        }
        if (numlines == 2)
        {
            std::string s;
            iss >> s >> xllCenter;
            numlines++;
            continue;
        }
        if (numlines == 3)
        {
            std::string s;
            iss >> s >> yllCenter;
            numlines++;
            continue;
        }
        if (numlines == 4)
        {
            std::string s;
            iss >> s >> cellSize;
            numlines++;
            continue;
        }
        if (numlines == 5)
        {
            std::string s;
            iss >> s >> nodataValue;
            numlines++;
            continue;
        }

        rows.push_back(std::vector<values_type>());
        for (auto i = 0; i < numX; ++i)
        {
            values_type value;
            iss >> value;
            rows.back().push_back(value);
        }
        numlines++;
    }

    std::vector<values_type> values;
    for (size_t i = 0; i < rows.size(); ++i)
    {
        for (auto j = 0u; j < rows[i].size(); ++j)
        {
            values.push_back(rows[i][j]);
        }
    }

    return {numX, numY, xllCenter, yllCenter, cellSize, nodataValue, values};
}
