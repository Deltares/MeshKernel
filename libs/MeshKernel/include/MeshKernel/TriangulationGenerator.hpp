//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2026.
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

#include <memory>

#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Polygons.hpp"

namespace meshkernel
{

    /// \brief Generate a triangulation from a polygonal domain boundary.
    class TriangulationGenerator
    {
    public:
        /// \brief Destructor
        virtual ~TriangulationGenerator() = default;

        /// \brief Compute triangulation
        virtual std::unique_ptr<Mesh2D> Generate(const Polygons& polygon) const = 0;
    };

    /// \brief Generate a triangulation using the triangle.c function
    class SimpleTriangulationGenerator : public TriangulationGenerator
    {
    public:
        /// \brief Constructor
        SimpleTriangulationGenerator(const double factor) : scaleFactor_(factor) {}

        /// \brief Compute points within polygon using triangle
        std::vector<Point> GeneratePoints(const Polygons& polygon) const;

        /// \brief Compute triangulation using triangle
        std::unique_ptr<Mesh2D> Generate(const Polygons& polygon) const override;

    private:
        /// \brief The scale factor used when generating points in polygon
        const double m_scaleFactor;
    };

    /// \brief Generate a triangulation using the SEPRAN library
    class SepranTriangulationGenerator : public TriangulationGenerator
    {
    public:
        /// \brief Compute triangulation using SEPRAN library
        std::unique_ptr<Mesh2D> Generate(const Polygons& polygon) const override;

    private:
        /// \brief Generate a vector of references to polygons from the set of polygonal enclosures
        static std::vector<std::reference_wrapper<const Polygon>> GeneratePolygonReferences(const Polygons& polygon);

        /// \brief Construct the set of edges, elements and number of nodes per element to construct mesh2d
        ///
        /// From the set of elements (triples of node ids in a flat array) construct arra of edges and elements
        static std::tuple<std::vector<Edge>, std::vector<std::vector<UInt>>, std::vector<std::uint8_t>>
        GatherEdgesAndFaces(const std::vector<int>& triangulationElementNodes, const int numberOfElements);

        /// \brief Construct arra of points from flat array of double (x,y values store in adjacent pairs)
        static std::vector<Point> PointsFromFlatArray(const std::vector<double>& coordinates, const int numberOfPoints);
    };

} // namespace meshkernel
