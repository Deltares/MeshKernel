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

#include "MeshKernel/BilinearInterpolationOnGriddedSamples.hpp"

#include <MeshKernel/Mesh2D.hpp>

using namespace meshkernel;

BilinearInterpolationOnGriddedSamples::BilinearInterpolationOnGriddedSamples(std::shared_ptr<Mesh2D> mesh,
                                                                             size_t numColumns,
                                                                             size_t numRows,
                                                                             double xOrigin,
                                                                             double yOrigin,
                                                                             double cellSize,
                                                                             std::vector<double>& values) : m_mesh(mesh),
                                                                                                            m_numColumns(numColumns),
                                                                                                            m_numRows(numRows),
                                                                                                            m_xOrigin(xOrigin),
                                                                                                            m_yOrigin(yOrigin),
                                                                                                            m_cellSize(cellSize),
                                                                                                            m_values(values) {}

void BilinearInterpolationOnGriddedSamples::Compute()
{
    m_nodeResults.resize(m_mesh->GetNumNodes());
    std::fill(m_nodeResults.begin(), m_nodeResults.end(), constants::missing::doubleValue);
    for (size_t n = 0; n < m_mesh->GetNumNodes(); ++n)
    {
        const auto node = m_mesh->m_nodes[n];
        m_nodeResults[n] = bilinearInterpolation(node);
    }

    m_edgeResults.resize(m_mesh->GetNumEdges());
    std::fill(m_edgeResults.begin(), m_edgeResults.end(), constants::missing::doubleValue);
    for (size_t e = 0; e < m_mesh->GetNumEdges(); ++e)
    {
        const auto& [first, second] = m_mesh->m_edges[e];
        m_edgeResults[e] = 0.5 * (m_nodeResults[first] + m_nodeResults[second]);
    }

    m_faceResults.resize(m_mesh->GetNumFaces(), constants::missing::doubleValue);
    std::fill(m_faceResults.begin(), m_faceResults.end(), constants::missing::doubleValue);
    for (size_t f = 0; f < m_mesh->GetNumFaces(); ++f)
    {
        m_faceResults[f] = bilinearInterpolation(m_mesh->m_facesMassCenters[f]);
    }
}

double BilinearInterpolationOnGriddedSamples::bilinearInterpolation(const Point& point) const
{
    double fractionalColumnIndex = GetFractionalNumberOfColumns(point);
    double fractionalRowIndex = GetFractionalNumberOfRows(point);

    const int columnIndex = static_cast<int>(fractionalColumnIndex);
    const int rowIndex = static_cast<int>(fractionalRowIndex);

    fractionalColumnIndex = fractionalColumnIndex - columnIndex;
    fractionalRowIndex = fractionalRowIndex - rowIndex;

    if (columnIndex < 0 || columnIndex >= m_numColumns || rowIndex < 0 || rowIndex >= m_numRows)
    {
        return constants::missing::doubleValue;
    }
    double result = fractionalColumnIndex * fractionalRowIndex * getGriddedValue(columnIndex + 1, rowIndex + 1) +
                    (1.0 - fractionalColumnIndex) * fractionalRowIndex * getGriddedValue(columnIndex, rowIndex + 1) +
                    (1.0 - fractionalColumnIndex) * (1.0 - fractionalRowIndex) * getGriddedValue(columnIndex, rowIndex) +
                    fractionalColumnIndex * (1.0 - fractionalRowIndex) * getGriddedValue(columnIndex + 1, rowIndex);

    return result;
}
