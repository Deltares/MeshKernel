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

namespace meshkernel
{
    // missing value
    const double doubleMissingValue = -999.0;
    const int intMissingValue = -999;
    const double squareRootOfThree = 1.73205080756887729352;
    const double oneThird = 0.33333333333333333333;

    //geometric constants
    const double degrad_hp = M_PI / 180.0;                   // conversion factor from degrees to radians(pi / 180)
    const double raddeg_hp = 180.0 / M_PI;                   // conversion factor from radians to degrees(180 / pi)
    const double earth_radius = 6378137.0;                   // earth radius(m)
    const double one_over_earth_radius = 1.0 / earth_radius; //one over earth_radius(m-1);
    const double absLatitudeAtPoles = 0.0001;                // pole tolerance in degrees
    const double nearlyZero = 1e-16;                         // used to determine if a length is zero
	const int numNodesInTriangle = 3;  

    //mesh constants
    const double minimumDeltaCoordinate = 1e-14;
    const int maximumNumberOfEdgesPerNode = 12;
    const int maximumNumberOfEdgesPerFace = 6;
    const int maximumNumberOfNodesPerFace = 8;
    const int maximumNumberOfConnectedNodes = maximumNumberOfEdgesPerNode * 4;
    const double minimumCellArea = 1e-12;
    const double weightCircumCenter = 1.0;
    const int numNodesQuads = 4;

    //orthogonalization
    const double minimumEdgeLength = 1e-4;
    const double curvilinearToOrthogonalRatio = 0.5;         //curvi - linear - like(0.0) or pure(1.0) orthogonalisation
    const double orthogonalizationToSmoothingFactor = 0.975; //Factor between grid smoothing and grid ortho resp (0.<=ATPF<=1.)

    // merging distance
    const double mergingDistance = 0.001;
    const double mergingDistanceSquared = mergingDistance * mergingDistance;

    // physical constants
    const double gravity = 9.81;

    // Operations averaging methods
    enum class AveragingMethod
    {
        SimpleAveraging = 1,
        ClosestPoint = 2,
        Max = 3,
        Min = 4,
        InverseWeightDistance = 5,
        MinAbs = 6,
        KdTree = 7
    };

    // Interpolation locations
    enum class InterpolationLocation
    {
        Faces = 0,
        Nodes = 1,
        Edges = 2
    };

} // namespace meshkernel
