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

#include "MeshKernel/Utilities/RTree.hpp"

#define BOOST_ALLOW_DEPRECATED_HEADERS
#include <boost/geometry.hpp>
#undef BOOST_ALLOW_DEPRECATED_HEADERS

namespace meshkernel
{
    /// @brief Factory class for creating instances of RTree with different projections.
    struct RTreeFactory
    {
        /// @brief Creates and returns a shared pointer to an instance of RTree with the specified projection.
        ///
        /// @param projection The projection type for the RTree.
        /// @return A shared pointer to the created RTree instance.
        ///
        /// This factory method creates an RTree instance based on the specified projection type.
        /// The supported projections are Cartesian, Spherical, and SphericalAccurate.
        ///
        /// @code
        /// auto cartesianRTree = RTreeFactory::create(Projection::cartesian);
        /// auto sphericalRTree = RTreeFactory::create(Projection::spherical);
        /// auto sphericalAccurateRTree = RTreeFactory::create(Projection::sphericalAccurate);
        /// @endcode
        ///
        /// @throws std::invalid_argument if an invalid projection value is provided.
        static std::shared_ptr<RTreeBase> create(Projection projection)
        {
            switch (projection)
            {
            case Projection::cartesian:
                return std::make_shared<RTree<bg::cs::cartesian>>();
            case Projection::spherical:
                return std::make_shared<RTree<bg::cs::geographic<bg::degree>>>();
            case Projection::sphericalAccurate:
                return std::make_shared<RTree<bg::cs::geographic<bg::degree>>>();
            default:
                throw std::invalid_argument("Invalid projection value: " + std::to_string(static_cast<int>(projection)));
            }
        }
    };

} // namespace meshkernel
