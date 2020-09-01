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
#include "Entities.hpp"
#include "CurvilinearParametersNative.hpp"
#include "SplinesToCurvilinearParametersNative.hpp"

namespace GridGeom
{
    class CurvilinearGrid;
    class Mesh;
    class LandBoundaries;

    class FlipEdges
    {
    public:

        /// <summary>
        /// Ctor
        /// </summary>
        /// <returns></returns>
        FlipEdges();

        /// <summary>
        /// 
        /// </summary>
        /// <param name="mesh">The mesh to use</param>
        /// <param name="landBoundary">The land boundaries</param>
        /// <returns></returns>
        FlipEdges(Mesh* mesh, LandBoundaries* landBoundary);

        /// <summary>
        /// Flip the edges
        /// </summary>
        /// <returns></returns>
        bool Compute();

    private:

        bool TriangulateFaces();

        bool TopologyFunctional( int edge,
                                 int& k1,
                                 int& k2,
                                 int& kl,
                                 int& kr,
                                 int& faceL,
                                 int& faceR,
                                 int& ntopo ) const;

        int OptimalNumberOfConnectedNodes(int nodeIndex) const;
       
        Mesh* m_mesh;                                      // A pointer to mesh
        LandBoundaries* m_landBoundaries;                   // A pointer to land boundaries
        bool m_triangulateFaces = false;
        bool m_projectToLandBoundary = false;

    };

}


