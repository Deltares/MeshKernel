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

#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/MeshInterpolationInterface.hpp>

namespace meshkernel
{
    class BilinearInterpolationOnGriddedSamples : public MeshInterpolationInterface
    {
    public:
        /// @brief Interpolation based on averaging
        /// @param[in] mesh                            The input mesh
        BilinearInterpolationOnGriddedSamples(std::shared_ptr<Mesh2D> mesh,
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

        /// @brief Compute interpolation
        void Compute()
        {
            std::vector nodalInterpolation(m_mesh->GetNumNodes(), constants::missing::doubleValue);
            std::fill(nodalInterpolation.begin(), nodalInterpolation.end(), constants::missing::doubleValue);
            for (size_t n = 0; n < m_mesh->GetNumNodes(); ++n)
            {
                const auto node = m_mesh->m_nodes[n];
                nodalInterpolation[n] = bilinearInterpolation(node);
            }

            m_results.resize(m_mesh->GetNumFaces(), constants::missing::doubleValue);
            std::fill(m_results.begin(), m_results.end(), constants::missing::doubleValue);

            for (size_t f = 0; f < m_mesh->GetNumFaces(); ++f)
            {
                double maxVal = -std::numeric_limits<double>::max();
                double minVal = std::numeric_limits<double>::max();
                for (size_t n = 0; n < m_mesh->GetNumFaceEdges(f); ++n)
                {
                    const auto node = m_mesh->m_facesNodes[f][n];
                    const auto val = nodalInterpolation[node];
                    maxVal = std::max(maxVal, val);
                    minVal = std::min(minVal, val);
                }
                m_results[f] = MaskValue(f, minVal, maxVal);
            }
        }

        [[nodiscard]] double GetResults(size_t index) const override
        {
            return m_results[index];
        }

        [[nodiscard]] size_t GetResultsSize() const override { return m_results.size(); }


    private:
        double MaskValue(size_t f, double minVal, double maxVal)
        {
            if (minVal > 0.0)
            {
                return std::numeric_limits<double>::max();
            }
            if (maxVal >= 0.0 && minVal <= 0.0)
            {
                return 0.0;
            }
            return bilinearInterpolation(m_mesh->m_facesMassCenters[f]);
        }

        double bilinearInterpolation(const Point& node) const
        {
            double fractionalColumnIndex = GetFractionalColumnIndex(node);
            double fractionalRowIndex = GetFractionalRowIndex(node);

            const int columnIndex = static_cast<int>(fractionalColumnIndex);
            const int rowIndex = static_cast<int>(fractionalRowIndex);

            fractionalColumnIndex = fractionalRowIndex - columnIndex;
            fractionalRowIndex = fractionalRowIndex - rowIndex;

            if (columnIndex < 0 || columnIndex >= m_numColumns || rowIndex < 0 && rowIndex >= m_numRows)
            {
                return constants::missing::doubleValue;
            }
            double result = fractionalColumnIndex * fractionalRowIndex * getGriddedValue(columnIndex + 1, rowIndex + 1) +
                            (1.0 - fractionalColumnIndex) * fractionalRowIndex * getGriddedValue(columnIndex, rowIndex + 1) +
                            (1.0 - fractionalColumnIndex) * (1.0 - fractionalRowIndex) * getGriddedValue(columnIndex, rowIndex) +
                            fractionalColumnIndex * (1.0 - fractionalRowIndex) * getGriddedValue(columnIndex + 1, rowIndex);

            return result;
        }

        [[nodiscard]] double GetFractionalColumnIndex(const Point& p) const
        {
            return (p.x - m_xOrigin) / m_cellSize;
        }

        [[nodiscard]] double GetFractionalRowIndex(const Point& p) const
        {
            return (p.y - m_yOrigin) / m_cellSize;
        }

        [[nodiscard]] double getGriddedValue(size_t rowIndex, size_t columnIndex) const
        {
            const auto index = rowIndex * m_numColumns + columnIndex;
            return m_values[rowIndex];
        }

        const std::shared_ptr<Mesh2D> m_mesh; ///< Pointer to the mesh

        size_t m_numColumns;
        size_t m_numRows;
        double m_xOrigin;
        double m_yOrigin;
        double m_cellSize;
        std::vector<double> m_values;

        std::vector<double> m_results; ///< The interpolation results
    };
} // namespace meshkernel
