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

extern "C"
{

    extern void mshoce_(
        const int* jnew,   // Input: Reset indicator (Fortran Logical pointer)
        double* coor,      // Output: Flattened array of node coordinates
        int* kmeshc,       // Output: Connectivity/topology grid matrix
        const int* inpelm, // Input: Element type identifier
        const int* nbound, // Input: Number of boundary elements
        double* bcord,     // Input: Coordinates of boundary control nodes

        int* kbndpt,            // Input: Type flags for boundary nodes
        int* boundary,          // Input: Edge-to-node connectivity map
        const int* numcurvboun, // Input: Total count of curved boundary segments
        int* npoint,            // Output: Count of generated points
        int* nelem,             // Output: Count of generated elements

        int* holeinfo,     // Input: Structural layout parameters for holes
        const int* nholes, // Input: Total count of internal holes
        const int* ncoar,  // Input: Quantity of sizing descriptors passed
        double* coar,      // Input: Target element sizing arrays
        int* userpoints,   // Input: Fixed target internal points

        int* isurnr,             // Output: Surface structural adjacency mapping register
        const int* numextcurves, // Input: Auxiliary alignment curve flags
        int* numnodextcurvs,     // Input: Node mappings for alignment paths

        int* curvenumbers, // Input: Curve curvature flags
        double* rinput,    // Input: Supplementary sizing/weighting matrices
        const int* nuspnt, // Input: Count of forced target control nodes
        const int* ndim    // Input: Domain spatial dimension identifier
    );
}

namespace meshkernel
{

    class TriangulationGenerator
    {
    public:
        virtual ~TriangulationGenerator() = default;

        virtual std::unique_ptr<Mesh2D> generate(const Polygons& polygon) const = 0;
    };

    /// \brief Generate a triangulation using the triangle.c function
    class SimpleTriangulationGenerator : public TriangulationGenerator
    {
    public:
        SimpleTriangulationGenerator(const double factor) : scaleFactor_(factor) {}

        std::unique_ptr<Mesh2D> generate(const Polygons& polygon) const override;

    private:
        const double scaleFactor_;
    };

    /// \brief Generate a triangulation using the SEPRAN library
    class SepranTriangulationGenerator : public TriangulationGenerator
    {
    public:
        std::unique_ptr<Mesh2D> generate(const Polygons& polygon) const override;

    private:
        static double minimumEdgeDelta(const std::vector<meshkernel::Point>& polygonNodes);

        static std::vector<std::reference_wrapper<const Polygon>> generatePolygonReferences(const Polygons& polygon);

        static std::tuple<std::vector<Edge>, std::vector<std::vector<UInt>>, std::vector<std::uint8_t>>
        gatherEdgesAndFaces(const std::vector<int>& kmeshc, const int numberOfElements);

        static std::vector<Point> pointsFromFlatArray(const std::vector<double>& coordinates, const int numberOfPoints);
    };

} // namespace meshkernel
