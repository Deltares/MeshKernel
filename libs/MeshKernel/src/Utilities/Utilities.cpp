//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2024.
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

#include "MeshKernel/Utilities/Utilities.hpp"
#include "MeshKernel/Exceptions.hpp"

#include <fstream>
#include <iomanip>

void meshkernel::Print(const std::vector<Point>& nodes, const std::vector<Edge>& edges, std::ostream& out)
{
    out << "nullId = " << constants::missing::uintValue << ";" << std::endl;
    out << "nullValue = " << constants::missing::doubleValue << ";" << std::endl;
    out << "nodex = zeros ( " << nodes.size() << ", 1);" << std::endl;
    out << "nodey = zeros ( " << nodes.size() << ", 1);" << std::endl;
    out << "edges = zeros ( " << edges.size() << ", 2);" << std::endl;

    for (UInt i = 0; i < nodes.size(); ++i)
    {
        out << "nodex (" << i + 1 << " ) = " << nodes[i].x << ";" << std::endl;
    }

    for (UInt i = 0; i < nodes.size(); ++i)
    {
        out << "nodey (" << i + 1 << " ) = " << nodes[i].y << ";" << std::endl;
    }

    out << "edges = zeros ( " << edges.size() << ", 2 );" << std::endl;

    for (UInt i = 0; i < edges.size(); ++i)
    {
        out << "edges ( " << i + 1 << ", 1 ) = " << edges[i].first + 1 << ";" << std::endl;
        out << "edges ( " << i + 1 << ", 2 ) = " << edges[i].second + 1 << ";" << std::endl;
    }
}

void meshkernel::Print(const std::vector<double>& xNodes,
                       const std::vector<double>& yNodes,
                       const std::vector<int>& edges,
                       std::ostream& out)
{
    // xNodes and yNodes should be the same size
    if (xNodes.size() != yNodes.size())
    {
        throw ConstraintError("x-node and y-nodes are not the same size, {} /= {}", xNodes.size(), yNodes.size());
    }

    out << "nullId = " << constants::missing::uintValue << ";" << std::endl;
    out << "nullValue = " << constants::missing::doubleValue << ";" << std::endl;
    out << "nodex = zeros ( " << xNodes.size() << ", 1);" << std::endl;
    out << "nodey = zeros ( " << xNodes.size() << ", 1);" << std::endl;
    out << "edges = zeros ( " << edges.size() << ", 2);" << std::endl;

    for (UInt i = 0; i < xNodes.size(); ++i)
    {
        out << "nodex (" << i + 1 << " ) = " << xNodes[i] << ";" << std::endl;
    }

    for (UInt i = 0; i < xNodes.size(); ++i)
    {
        out << "nodey (" << i + 1 << " ) = " << yNodes[i] << ";" << std::endl;
    }

    out << "edges = zeros ( " << edges.size() / 2 << ", 2 );" << std::endl;

    for (UInt i = 0; i < edges.size() / 2; ++i)
    {
        out << "edges ( " << i + 1 << ", 1 ) = " << edges[2 * i] + 1 << ";" << std::endl;
        out << "edges ( " << i + 1 << ", 2 ) = " << edges[2 * i + 1] + 1 << ";" << std::endl;
    }
}

void meshkernel::SaveVtk(const std::vector<Point>& nodes, const std::vector<std::vector<UInt>>& faces, const std::string& fileName)
{

    std::ofstream vtkFile(fileName.c_str());

    vtkFile.precision(18);

    std::string meshType = "UnstructuredGrid";
    std::string versionNumber = "1.0";
    std::array<UInt, 15> unsavedElements;

    unsavedElements.fill (0);

    UInt numberOfElements = 0;

    for (size_t i = 0; i < faces.size(); ++i)
    {
        if (faces[i].size() == 3 || faces[i].size() == 4)
        {
            ++numberOfElements;
        }

        if (faces[i].size() < 15) {
            ++unsavedElements[faces[i].size()];
        }

    }

    for (size_t i = 0; i < unsavedElements.size (); ++i){
        std::cout << "elements " << i <<  " = " << unsavedElements [i] << std::endl;
    }

    vtkFile << "<VTKFile type=\""
            << meshType
            << "\" version=\""
            << versionNumber
            << "\" byte_order=\"LittleEndian\" header_type=\"UInt64\">"
            << std::endl;
    vtkFile << "  <UnstructuredGrid>" << std::endl;

    vtkFile << "    <Piece"
            << "  NumberOfPoints=\""
            << nodes.size() << "\""
            << "  NumberOfCells=\""
            << numberOfElements
            << "\">"
            << std::endl;

    vtkFile << "      <Points>" << std::endl;
    vtkFile << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

    for (size_t i = 0; i < nodes.size(); ++i)
    {
        vtkFile << std::setw(8) << " " << nodes[i].x << "  " << std::setw(8) << nodes[i].y << "      " << 0.0 << std::endl;
    }

    vtkFile << "        </DataArray>" << std::endl;
    vtkFile << "      </Points>" << std::endl;

    //--------------------------------
    vtkFile << "      <Cells>" << std::endl;
    vtkFile << "        <DataArray type=\""
            << "Int64"
            << "\" Name=\""
            << "connectivity"
            << "\" format=\""
            << "ascii"
            << "\">"
            << std::endl;

    for (size_t i = 0; i < faces.size(); ++i)
    {

        if (faces[i].size() == 3 || faces[i].size() == 4)
        {
            for (size_t j = 0; j < faces[i].size(); ++j)
            {
                vtkFile << " " << std::setw(8) << faces[i][j];
            }

            vtkFile << std::endl;
        }
    }

    vtkFile << std::endl;
    vtkFile << "        </DataArray>" << std::endl;

    //--------------------------------
    vtkFile << "        <DataArray type=\""
            << "Int64"
            << "\" Name=\""
            << "offsets"
            << "\" format=\""
            << "ascii"
            << "\">"
            << std::endl;

    size_t sum = 0;

    for (size_t i = 0; i < faces.size(); ++i)
    {
        if (faces[i].size() == 3 || faces[i].size() == 4)
        {
            sum += faces[i].size();
            vtkFile << " " << std::setw(8) << sum;
        }
    }

    vtkFile << std::endl;
    vtkFile << "        </DataArray>" << std::endl;

    //--------------------------------
    vtkFile << "        <DataArray type=\""
            << "Int64"
            << "\" Name=\""
            << "types"
            << "\" format=\""
            << "ascii"
            << "\">"
            << std::endl;

    for (size_t i = 0; i < faces.size(); ++i)
    {
        UInt type = 0;

        if (faces[i].size() == 3)
        {
            type = 69;
        }
        else if (faces[i].size() == 4)
        {
            type = 70;
        }

        if (faces[i].size() == 3 || faces[i].size() == 4)
        {
            vtkFile << " " << std::setw(8) << type;
        }
    }

    vtkFile << std::endl;
    vtkFile << "        </DataArray>" << std::endl;
    vtkFile << "      </Cells>" << std::endl;

    vtkFile << "    </Piece>" << std::endl;
    vtkFile << "  </UnstructuredGrid>" << std::endl;
    vtkFile << "</VTKFile>" << std::endl;

    vtkFile.close();
}
