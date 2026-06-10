// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 Stichting Deltares
//
// C++20 translation of msho01, msho09, msho11-msho14, msho17, msho20-msho21,
// msho24-msho30, msho33-msho35, msho36, msho38-msho42
// from the SEPRAN library ("Ingenieursbureau SEPRA", Niek Praagman, 1989-2010).

#pragma once

#include "SepranContext.hpp"

#include <span>
#include <string_view>

namespace sepran
{

// ---------------------------------------------------------------------------
// Array/index conventions (used throughout this module):
//   coor   – interleaved flat: coor[2*(k-1)] = x, coor[2*(k-1)+1] = y (1-based k)
//   kstapl – flat edge pairs: kstapl[2*s]=n1, kstapl[2*s+1]=n2 (0-based s)
//   kelem  – flat triples: kelem[3*e], [3*e+1], [3*e+2] (0-based e)
//   itri   – flat: itri[k-1] for 1-based node k
//   istart – CSR: istart[i-1] = cumulative neighbour count through node i
// ---------------------------------------------------------------------------

/// \brief Check boundary elements and fill per-node coarsenesses.
///
/// Translated from Fortran \c msho01 (SEPRAN, Niek Praagman, 1989-2008).
void msho01(std::span<const int>    kbound,
            int                     nbound,
            std::span<int>          istart,
            std::span<int>          ibuur,
            std::span<double>       coarse,
            std::span<const double> coor,
            int                     npoint,
            double&                 coarsemin,
            double&                 coarsemax,
            std::span<double>       coar,
            int                     ncoar,
            SepranContext&          ctx);

/// \brief Compute the unit perpendicular and midpoint of segment i→j.
///
/// Translated from Fortran \c msho09.
void msho09(std::span<const double> coor,
            int i, int j,
            double& e1, double& e2,
            double& xm, double& ym);

/// \brief Compute the cosine of the angle at vertex i2 between lines i2→i1 and i2→i3.
///
/// Translated from Fortran \c msho11.
void msho11(int i1, int i2, int i3,
            std::span<const double> coor,
            double& angle,
            const SepranContext& ctx);

/// \brief Find the best adjacent front edges at the endpoints of base edge i1–i2.
///
/// Translated from Fortran \c msho12.
void msho12(std::span<const double> coor,
            std::span<const int>    kstapl,
            int                     kstap,
            int i1, int i2,
            int& iex1, int& iex2,
            double& angle1, double& angle2,
            SepranContext& ctx);

/// \brief Check whether lines from i1 or i2 to new point (xn,yn) cross any front edge.
///
/// Translated from Fortran \c msho13.
/// kdrie: 0=none, -1=ambiguous, >0 = 1-based index of crossing edge.
void msho13(std::span<const double> coor,
            int i1, int i2,
            int& kdrie,
            std::span<const int> kstapl,
            int kstap,
            double xn, double yn,
            const SepranContext& ctx);

/// \brief Find the nearest front node to (xn,yn).
///
/// Translated from Fortran \c msho14.
void msho14(std::span<const double> coor,
            int& jpn,
            int npoint,
            std::span<const int> itri,
            int i1, int i2,
            double xn, double yn,
            double& dista);

/// \brief Commit a new node and triangle, update kstapl and itri.
///
/// Translated from Fortran \c msho15.
void msho15(std::span<double> coor,
            int& npoint,
            std::span<int> kelem,
            int& nelem,
            std::span<int> kstapl,
            int& kstap,
            std::span<int> itri,
            int i1, int i2,
            double xnx, double yny);

/// \brief Build the full CSR neighbour structure for all mesh nodes.
///
/// Translated from Fortran \c msho16.
void msho16(std::span<const int> kelem,
            int nelem, int npoint, int nipnt,
            std::span<int> iwork,
            std::span<int> ibuurp,
            int& leng,
            SepranContext& ctx);

/// \brief Insert neighbour ij into the CSR neighbour list for node ih.
///
/// Translated from Fortran \c msho17.
void msho17(std::span<int>       ibuurp,
            std::span<const int> iwork,
            int ih, int ij);

/// \brief Rotate kstapl left by one edge (first edge moved to end) and update itri.
///
/// Translated from Fortran \c msho20.
void msho20(std::span<int> kstapl,
            int kstap,
            std::span<int> itri);

/// \brief Compute reference triangle-area values for each coarseness-grid cell.
///
/// Translated from Fortran \c msho21.
void msho21(std::span<const double> cube,
            int ncube,
            std::span<double> refvol);

/// \brief Check whether segment i1–i2 crosses any front edge.
///
/// Translated from Fortran \c msho24.
/// icheck: -1 = Plaxis duplicate coincidence; 0 = clear; >0 = crossing edge index.
void msho24(std::span<const int>    kstapl,
            int kstap,
            std::span<const double> coor,
            int i1, int i2,
            int& icheck,
            const SepranContext& ctx);

/// \brief Move the shortest front edge to position 0 in kstapl.
///
/// Translated from Fortran \c msho25.
void msho25(std::span<int>          kstapl,
            int kstap,
            std::span<const double> coor,
            const SepranContext& ctx);

/// \brief Find the triangle (other than jelem) sharing edge i2–i3.
///
/// Translated from Fortran \c msho26.
void msho26(std::span<const int> kelem,
            int nelem,
            int i2, int i3,
            int& jelem,
            int& i4);

/// \brief Select the best diagonal in quadrilateral (i1,i2,i3,i4) using Delaunay criterion.
///
/// Translated from Fortran \c msho27.
void msho27(std::span<const double> coor,
            int& i1, int& i2, int& i3, int& i4,
            const SepranContext& ctx);

/// \brief Check whether any front node lies strictly inside triangle (i1,i2,i3).
///
/// Translated from Fortran \c msho28.
void msho28(std::span<const double> coor,
            int i1, int i2, int i3,
            double xn, double yn,
            int npoint,
            std::span<const int> itri,
            int& iperm,
            const SepranContext& ctx);

/// \brief Remove over-refined interior nodes with 3 or 4 neighbours.
///
/// Translated from Fortran \c msho29.
void msho29(std::span<int>          kelem,
            int&                    nelem,
            int                     npoint,
            std::span<const int>    iwork,
            std::span<const int>    ibuurp,
            int                     nipnt,
            std::span<const double> coor,
            int&                    icancel);

/// \brief Recover from triangulation error: remove ielem and rebuild front.
///
/// Translated from Fortran \c msho30.
void msho30(std::span<int>  kelem,
            int&            nelem,
            int             ielem,
            std::span<int>  kstapl,
            int&            kstap,
            int             npoint,
            std::span<int>  itri,
            int             nelmfix);

/// \brief Compute the triangle quality ratio 2·r_in/r_out.
///
/// Translated from Fortran \c msho33.
void msho33(std::span<const double> coor,
            double& ratio,
            int i1, int i2, int i3);

/// \brief Locate user-point curve endpoints for a given node.
///
/// Translated from Fortran \c msho34.
void msho34(int                     kpoint,
            std::span<const double> coor,
            int                     ncurvs,
            std::span<const int>    curves,
            std::span<const int>    kbndpt,
            int                     nbndpt,
            int                     numcurvboun,
            std::span<const int>    boundary,
            int                     nbound,
            int                     nholes,
            std::span<const double> userco,
            int                     nuspnt,
            std::span<const double> coaval,
            int&                    ius1,
            int&                    ius2);

/// \brief Check coarseness compatibility and clamp cmax if transition is infeasible.
///
/// Translated from Fortran \c msho39.
void msho39(double& cmax,
            double  cmin,
            double  dist,
            double  maxratio,
            int&    iallow);

/// \brief Compute the maximum allowed coarseness over distance dist from cmin.
///
/// Translated from Fortran \c msho41.
void msho41(double cmin, double& cmax, double dist, double maxratio);

/// \brief Check whether edge i–j is an internal boundary segment.
///
/// Translated from Fortran \c msho42.
void msho42(std::span<const int> kbndpt,
            int lenbnd,
            int i, int j,
            int& ja);

// ---------------------------------------------------------------------------
// Wrappers for msho04, msho05, msho06, msho07, msho08 – larger functions
// defined in SepranFront.cpp
// ---------------------------------------------------------------------------

/// \brief Compute bounding box of npoint nodes.
///
/// Translated from Fortran \c msho04.
void msho04(std::span<const double> coor,
            int npoint,
            double& xmin, double& xmax,
            double& ymin, double& ymax,
            const SepranContext& ctx);

/// \brief Compute extreme coarsenesses.
///
/// Translated from Fortran \c msho05.
void msho05(std::span<const double> dist,
            int npoint,
            double& dismin, double& dismax,
            const SepranContext& ctx);

/// \brief Determine per-quad coarsenesses and type flags.
///
/// Translated from Fortran \c msho06.
void msho06(int npoint,
            std::span<const double> coor,
            double dist, double xstart, double ystart,
            int nx, int ny,
            std::span<int> icube,
            std::span<const double> chelp,
            std::span<double> cube,
            std::span<int> jcube,
            std::span<const int> kbound,
            int nbound,
            std::span<const double> coar,
            int ncoar,
            int ncurvs,
            std::span<const int> curves,
            std::span<const double> cocurvs,
            const SepranContext& ctx);

/// \brief Point-in-polygon test (ray casting).
///
/// Translated from Fortran \c msho07.
void msho07(double xcub, double ycub,
            double xmini,
            std::span<const double> coor,
            std::span<const int> kbound,
            int nbound,
            int& ja,
            const SepranContext& ctx);

/// \brief Check and fix orientation of all front edges.
///
/// Translated from Fortran \c msho08.
void msho08(std::span<int> kstapl,
            int kstap,
            std::span<const double> coor,
            double xstart, double dismin,
            std::span<int> holeinfo,
            int nholes,
            bool check,
            SepranContext& ctx);

/// \brief Place fixed hexagonal clusters for user coarseness points.
///
/// Translated from Fortran \c msho35.
void msho35(int npoint,
            std::span<double> coor,
            double xstart, double ystart, double dismax,
            std::span<const double> coar, int ncoar,
            std::span<int> icube, int nx,
            std::span<int> kelem, int& nelem,
            std::span<int> kstapl, int& kstap,
            std::span<int> itri,
            int isurnr,
            std::span<const int> userpoints,
            std::span<int> kbndpt, int& nbndpt,
            SepranContext& ctx);

/// \brief Register internal curve nodes into coor, kstapl, kbndpt.
///
/// Translated from Fortran \c msho36.
void msho36(std::span<double> coor,
            int& npoint,
            int istep, int ncurvs,
            std::span<const int> curves,
            std::span<const double> cocurvs,
            std::span<int> kstapl, int& kstap,
            std::span<int> kbndpt, int& nbndpt,
            double coarsemin,
            std::span<int> extquanodes,
            std::span<double> coarse,
            SepranContext& ctx);

/// \brief Check/adjust coarseness smoothness.
///
/// Translated from Fortran \c msho38.
void msho38(int npoint,
            std::span<const double> coor,
            double dist, double xstart, double ystart,
            std::span<const int> kbndpt, int nbndpt,
            int numcurvboun,
            std::span<const int> boundary,
            int nbound, int nholes,
            int nx, int ny,
            std::span<const int> icube,
            std::span<const double> coarse,
            std::span<const double> cube,
            std::span<const int> jcube,
            std::span<const int> kstapl, int kstap,
            int ncurvs,
            std::span<const int> curves,
            std::span<const double> cocurv,
            int isurnr,
            std::span<const double> userco,
            int nuspnt,
            std::span<double> coaval,
            double tran,
            const SepranContext& ctx);

/// \brief Full mesh self-intersection check.
///
/// Translated from Fortran \c msho40.
/// iperm = 1 if intersections found, 0 if clean.
void msho40(std::span<const double> coor,
            std::span<const int> kelem,
            int npoint, int nelem,
            int& iperm);

} // namespace sepran
