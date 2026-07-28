//----- AGPL --------------------------------------------------------------------
//
//  Copyright (C)  Stichting Deltares, 2017-2026.
//
//  This programme is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Affero General Public License as
//  published by the Free Software Foundation version 3.
//
//  This programme is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Affero General Public License for more details.
//
//  You should have received a copy of the GNU Affero General Public License
//  along with This programme. If not, see <http://www.gnu.org/licenses/>.
//
//  contact: delft3d.support@deltares.nl
//  Stichting Deltares
//  P.O. Box 177
//  2600 MH Delft, The Netherlands
//
//  All indications and logos of, and references to, "Delft3D",
//  "D-Flow Flexible Mesh" and "Deltares" are registered trademarks of Stichting
//  Deltares, and remain the property of Stichting Deltares. All rights reserved.
//
//-------------------------------------------------------------------------------
//
// C++20 translation of msho02, msho03, msho10, msho18, msho19, msho22
// from the SEPRAN library.

#pragma once

#include "SepranContext.hpp"

#include <span>

namespace sepran
{

// ---------------------------------------------------------------------------
// Array conventions used throughout sepran topology:
//
//  coor  – interleaved flat: coor[2*(k-1)] = x, coor[2*(k-1)+1] = y
//          k is a 1-based node index.
//
//  istart – CSR row-pointer array of length npoint.
//           istart[i-1] (1-based i) = total number of neighbours for nodes 1..i.
//           Neighbours of node i are at ibuur[istart[i-2] .. istart[i-1]-1]
//           (with istart[-1] := 0).
//
//  ibuur  – CSR adjacency list; ibuur[j] is the 1-based node index of a
//           neighbour.
//
//  kelem  – element (triangle) array, flat: kelem[3*e], kelem[3*e+1],
//           kelem[3*e+2] are the 1-based node indices of element e (0-based e).
//
//  kstapl – advancing-front edge list, flat pairs: kstapl[2*k]=n1,
//           kstapl[2*k+1]=n2 (1-based node indices), k is 0-based.
//
//  itri   – per-node front-membership counter; itri[k-1] for 1-based node k.
// ---------------------------------------------------------------------------

/// \brief Insert neighbour j into the CSR adjacency list of node i.
///
/// Translated from Fortran \c msho02 (SEPRAN, Niek Praagman, 1989-1997).
///
/// \param istart  CSR row-pointer array (length npoint, 1-based nodes).
/// \param ibuur   CSR adjacency list (pre-allocated).
/// \param i       Source node (1-based).
/// \param j       Neighbour node to insert (1-based).
void insertNeighbour(std::span<const int> istart,
                     std::span<int> ibuur,
                     int i,
                     int j);

/// \brief Compute the Euclidean distance between two nodes.
///
/// Translated from Fortran \c msho03 (SEPRAN, Niek Praagman, 1989-1994).
///
/// \param i     First node (1-based).
/// \param j     Second node (1-based).
/// \param coor  Node coordinate array.
/// \return      Euclidean distance ||node_i - node_j||.
double nodeDistance(int i, int j, std::span<const double> coor);

/// \brief Add a new triangle element and update the advancing-front arrays.
///
/// Translated from Fortran \c msho10 (SEPRAN, Niek Praagman, 1989-2003).
///
/// Appends element (i1, i2, jpn) to \p kelem, removes edges i2-jpn and
/// jpn-i1 from \p kstapl if they were already on the front, or adds them if
/// not, and updates \p itri accordingly.
///
/// \param kelem   Element array (flat, 3 entries per element).  Must be
///                pre-allocated to maximum capacity.
/// \param nelem   Current element count (in-out).
/// \param i1      First node of new triangle (1-based).
/// \param i2      Second node of new triangle (1-based).
/// \param jpn     Third node of new triangle (1-based).
/// \param kstapl  Advancing-front edge list (flat pairs, 1-based nodes).
/// \param kstap   Current length of kstapl (number of edges) (in-out).
/// \param itri    Per-node front-membership counter (1-based index − 1).
void addElement(std::span<int> kelem, int& nelem,
                int i1, int i2, int jpn,
                std::span<int> kstapl, int& kstap,
                std::span<int> itri);

/// \brief Compute the signed area of triangle (i1, i2, i3).
///
/// Translated from Fortran \c msho19 (SEPRAN, Niek Praagman, 1989-2005).
///
/// Returns a positive value for counter-clockwise orientation.
///
/// \param coor  Node coordinate array.
/// \param i1, i2, i3  Node indices (1-based).
/// \return Signed area = ((x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2)) / 2).
double triangleArea(std::span<const double> coor, int i1, int i2, int i3);

/// \brief Check triangle angles; return which angle (if any) exceeds 90°.
///
/// Translated from Fortran \c msho22 (SEPRAN, Niek Praagman, 1989-2010).
///
/// \param a  Length of side a.
/// \param b  Length of side b.
/// \param c  Length of side c.
/// \return   0 if all angles ≤ 90°; 1/2/3 indicating which vertex has the
///           obtuse angle.
int checkTriangleAngles(double a, double b, double c);

/// \brief Laplacian smoothing of interior nodes.
///
/// Translated from Fortran \c msho18 (SEPRAN, Niek Praagman, 1989-2010).
///
/// Iteratively moves interior nodes (indices nipnt+1 .. npoint, 1-based)
/// toward the barycentre of their neighbours, accepting moves that preserve
/// positive triangle areas.
///
/// \param kelem   Element array (flat, 3 entries per element).
/// \param nelem   Number of elements.
/// \param coor    Node coordinate array (in-out).
/// \param npoint  Total node count.
/// \param nipnt   Number of boundary (fixed) nodes (1-based indices 1..nipnt).
/// \param iwork   CSR row-pointer array for ibuurp (length npoint).
/// \param ibuurp  CSR adjacency list of all nodes (length iwork[npoint-1]).
/// \param coars   Reference coarseness / stopping distance (in-out; halved inside).
/// \param ctx     SEPRAN context.
void laplacianSmoothing(std::span<const int> kelem, int nelem,
                        std::span<double> coor, int npoint, int nipnt,
                        std::span<const int> iwork,
                        std::span<const int> ibuurp,
                        double& coars,
                        SepranContext& ctx);

} // namespace sepran
