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

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif
#include <math.h> 

namespace GridGeom
{
    // missing value
    static constexpr double doubleMissingValue = -999.0;
    static constexpr int    intMissingValue    = -999;
    static constexpr double squareRootOfThree  = 1.73205080756887729352;
    static constexpr double oneThird           = 0.33333333333333333333;

    //geometric constants
    static constexpr double degrad_hp = M_PI / 180.0; // conversion factor from degrees to radians(pi / 180)
    static constexpr double raddeg_hp = 180.0 / M_PI; // conversion factor from radians to degrees(180 / pi)
    static constexpr double earth_radius = 6378137.0; // earth radius(m)
    static constexpr double one_over_earth_radius = 1.0 / earth_radius; //one over earth_radius(m-1);
    static constexpr double absLatitudeAtPoles = 0.0001;       // pole tolerance in degrees
    static constexpr double nearlyZero = 1e-16;                // used to determine if a length is zero

    //mesh constants
    static constexpr double minimumDeltaCoordinate = 1e-14;
    static constexpr int maximumNumberOfEdgesPerNode = 12;
    static constexpr int maximumNumberOfEdgesPerFace = 6;
    static constexpr int maximumNumberOfNodesPerFace = 8;
    static constexpr int maximumNumberOfConnectedNodes = maximumNumberOfEdgesPerNode * 4;
    static constexpr double minimumCellArea = 1e-12;
    static constexpr double weightCircumCenter = 1.0;
    static constexpr int numNodesQuads = 4;

    //orthogonalization 
    static constexpr double minimumEdgeLength = 1e-4;
    static constexpr double curvilinearToOrthogonalRatio= 0.5; //curvi - linear - like(0.0) or pure(1.0) orthogonalisation
    static constexpr double orthogonalizationToSmoothingFactor = 0.975; //Factor between grid smoothing and grid ortho resp (0.<=ATPF<=1.)

    // merging distance
    static constexpr double mergingDistance = 0.001;
    static constexpr double mergingDistanceSquared = mergingDistance * mergingDistance;

    // physical constants
    static constexpr double gravity = 9.81;

    // Operations averaging methods
    enum AveragingMethod 
    {
        SimpleAveraging = 1,
        ClosestPoint= 2,
        Max = 3,
        Min = 4,
        InverseWeightDistance = 5,
        MinAbs = 6,
        KdTree = 7
    };

}