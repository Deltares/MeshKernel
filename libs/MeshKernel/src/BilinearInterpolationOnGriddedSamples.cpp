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

#include "MeshKernel/BilinearInterpolationOnGriddedSamples.hpp"

#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/MeshInterpolation.hpp>

using namespace meshkernel;

template <typename ValueType>
BilinearInterpolationOnGriddedSamples<ValueType>::BilinearInterpolationOnGriddedSamples(const Mesh2D& mesh,
                                                                                        UInt numXCoord,
                                                                                        UInt numYCoord,
                                                                                        const Point& origin,
                                                                                        double cellSize,
                                                                                        const void* values) : m_mesh(mesh),
                                                                                                              m_numXCoord(numXCoord),
                                                                                                              m_numYCoord(numYCoord),
                                                                                                              m_origin(origin),
                                                                                                              m_cellSize(cellSize),
                                                                                                              m_isCellSizeConstant(true)
{
    m_values = static_cast<const ValueType*>(values);
}

template <typename ValueType>
BilinearInterpolationOnGriddedSamples<ValueType>::BilinearInterpolationOnGriddedSamples(const Mesh2D& mesh,
                                                                                        const std::vector<double>& xCoordinates,
                                                                                        const std::vector<double>& yCoordinates,
                                                                                        const void* values) : m_mesh(mesh),
                                                                                                              m_numXCoord(static_cast<UInt>(xCoordinates.size())),
                                                                                                              m_numYCoord(static_cast<UInt>(yCoordinates.size())),
                                                                                                              m_xCoordinates(xCoordinates),
                                                                                                              m_yCoordinates(yCoordinates),
                                                                                                              m_isCellSizeConstant(false)
{
    m_values = static_cast<const ValueType*>(values);
}

template <typename ValueType>
void BilinearInterpolationOnGriddedSamples<ValueType>::Compute()
{
    const auto numNodes = m_mesh.GetNumNodes();
    const auto numEdges = m_mesh.GetNumEdges();
    const auto numFaces = m_mesh.GetNumFaces();

    m_nodeResults.resize(numNodes);
    std::ranges::fill(m_nodeResults, constants::missing::doubleValue);
    for (UInt n = 0; n < numNodes; ++n)
    {
        const auto node = m_mesh.m_nodes[n];
        m_nodeResults[n] = Interpolation(node);
    }

    m_edgeResults.resize(numEdges);
    std::ranges::fill(m_edgeResults, constants::missing::doubleValue);
    for (UInt e = 0; e < numEdges; ++e)
    {
        const auto& [first, second] = m_mesh.m_edges[e];
        m_edgeResults[e] = 0.5 * (m_nodeResults[first] + m_nodeResults[second]);
    }

    m_faceResults.resize(numFaces, constants::missing::doubleValue);
    std::ranges::fill(m_faceResults, constants::missing::doubleValue);
    for (UInt f = 0; f < numFaces; ++f)
    {
        m_faceResults[f] = Interpolation(m_mesh.m_facesMassCenters[f]);
    }
}

template <typename ValueType>
double BilinearInterpolationOnGriddedSamples<ValueType>::Interpolation(const Point& point) const
{

    double fractionalColumnIndex = GetFractionalNumberOfColumns(point);
    double fractionalRowIndex = GetFractionalNumberOfRows(point);

    double columnIndexTmp;
    fractionalColumnIndex = std::modf(fractionalColumnIndex, &columnIndexTmp);

    double rowIndexTmp;
    fractionalRowIndex = std::modf(fractionalRowIndex, &rowIndexTmp);

    if (columnIndexTmp < 0 || rowIndexTmp < 0)
    {
        return constants::missing::doubleValue;
    }

    UInt const columnIndex = static_cast<UInt>(columnIndexTmp);
    UInt const rowIndex = static_cast<UInt>(rowIndexTmp);

    if (columnIndex + 1 >= m_numXCoord || rowIndex + 1 >= m_numYCoord)
    {
        return constants::missing::doubleValue;
    }

    const auto result = fractionalColumnIndex * fractionalRowIndex * GetGriddedValue(columnIndex + 1, rowIndex + 1) +
                        (1.0 - fractionalColumnIndex) * fractionalRowIndex * GetGriddedValue(columnIndex, rowIndex + 1) +
                        (1.0 - fractionalColumnIndex) * (1.0 - fractionalRowIndex) * GetGriddedValue(columnIndex, rowIndex) +
                        fractionalColumnIndex * (1.0 - fractionalRowIndex) * GetGriddedValue(columnIndex + 1, rowIndex);
    return result;
}

template <typename ValueType>
[[nodiscard]] double BilinearInterpolationOnGriddedSamples<ValueType>::GetFractionalNumberOfColumns(const Point& point) const
{
    if (m_isCellSizeConstant)
    {
        return (point.x - m_origin.x) / m_cellSize;
    }
    double result = constants::missing::doubleValue;
    if (m_xCoordinates.size() < 2)
    {
        return result;
    }
    for (UInt i = 0; i < m_xCoordinates.size() - 1; ++i)
    {

        if (point.x >= m_xCoordinates[i] && point.x < m_xCoordinates[i + 1])
        {
            const double dx = m_xCoordinates[i + 1] - m_xCoordinates[i];
            result = static_cast<double>(i) + (point.x - m_xCoordinates[i]) / dx;
            break;
        }
    }
    return result;
}

template <typename ValueType>
double BilinearInterpolationOnGriddedSamples<ValueType>::GetFractionalNumberOfRows(const Point& point) const
{
    if (m_isCellSizeConstant)
    {
        return (point.y - m_origin.y) / m_cellSize;
    }

    double result = constants::missing::doubleValue;
    if (m_yCoordinates.size() < 2)
    {
        return result;
    }

    for (UInt i = 0; i < m_yCoordinates.size() - 1; ++i)
    {

        if (point.y >= m_yCoordinates[i] && point.y < m_yCoordinates[i + 1])
        {
            const double dy = m_yCoordinates[i + 1] - m_yCoordinates[i];
            result = static_cast<double>(i) + (point.y - m_yCoordinates[i]) / dy;
            break;
        }
    }
    return result;
}

template class BilinearInterpolationOnGriddedSamples<float>;
template class BilinearInterpolationOnGriddedSamples<double>;
template class BilinearInterpolationOnGriddedSamples<short>;
