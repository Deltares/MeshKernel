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

#include <MeshKernel/Entities.hpp>

/// \namespace meshkernel
/// @brief Contains the logic of the C++ static library
namespace meshkernel
{
    /// @brief A class describing a network 1d.
    ///
    /// A network 1d contains a series of polylines and a series of fixed points at predefined chainages.
    class Network1D
    {
    public:
        /// @brief Default constructor
        Network1D() = default;

        /// @brief Construct a mesh1d by discretizing polyLines
        /// @param[in] polyLines The polylines to be discretize
        /// @param[in] fixedChainages The fixed chainages. These are locations where the first discretization points before and after should be at a distance equal to \ref offsetFixedChainages.
        /// @param[in] projection The projection to use
        explicit Network1D(std::vector<std::vector<Point>> const& polyLines,
                           std::vector<std::vector<double>> const& fixedChainages,
                           Projection projection);

        /// @brief Compute the chainages from fixed point locations
        /// @param[in] minChainageLength   The minimum face size. The distance between two discratization point must be no less than this length
        /// @param[in] offsetFromFixedChainages  The offset to use for fixed chainages
        void ComputeFixedChainages(double minChainageLength,
                                   double offsetFromFixedChainages);

        /// @brief Compute the chainages at a regular offset
        /// @param[in] polyLineIndex The index of the polyline we are operating on
        /// @param[in] offset The offset between points
        void ComputeOffsettedChainages(double offset);

        /// @brief Computes the discretization points from the chainages
        /// @return A discretization vector for each polyline
        [[nodiscard]] std::vector<std::vector<Point>> ComputeDiscretizationsFromChainages();

        Projection m_projection; ///< The projection used

    private:
        std::vector<std::vector<Point>> m_polyLines;       ///< The network poly lines
        std::vector<std::vector<double>> m_fixedChainages; ///< The fixed chainages at each poly line
        std::vector<std::vector<double>> m_chainages;      ///< The computed chainages
    };

} // namespace meshkernel
