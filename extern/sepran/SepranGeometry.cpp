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
// C++20 translation of mshcrossline.for, mshcrossline1.for, and msho75.for
// from the SEPRAN library.
//
// Original authors:
//   mshcrossline  – Guus Segal, 2000-2008
//   mshcrossline1 – Guus Segal, 2003
//   msho75        – Niek Praagman, 1996-2008

#include "SepranGeometry.hpp"

#include <algorithm>
#include <cmath>

namespace sepran
{

// ---------------------------------------------------------------------------
// crossLine  (mshcrossline.for)
// ---------------------------------------------------------------------------

void crossLine(const std::array<double, 2>& x1,
               const std::array<double, 2>& x2,
               const std::array<double, 2>& x3,
               const std::array<double, 2>& x4,
               double& fact1,
               double& fact2,
               double eps,
               const SepranContext& ctx)
{
    const double det = (x2[0] - x1[0]) * (x4[1] - x3[1])
                     - (x4[0] - x3[0]) * (x2[1] - x1[1]);

    const double amax = std::max({std::abs(x1[0]), std::abs(x1[1]),
                                  std::abs(x2[0]), std::abs(x2[1]),
                                  std::abs(x3[0]), std::abs(x3[1]),
                                  std::abs(x4[0]), std::abs(x4[1])});

    if (std::abs(det) <= ctx.epsmac * amax * 100.0)
    {
        fact1 = -ctx.rinfin;
        fact2 = -ctx.rinfin;
    }
    else
    {
        fact1 = ((x4[1] - x3[1]) * (x3[0] - x1[0])
               + (x3[0] - x4[0]) * (x3[1] - x1[1])) / det;

        if (fact1 >= -eps && fact1 <= 1.0 + eps)
        {
            fact2 = -((x1[1] - x2[1]) * (x3[0] - x1[0])
                    + (x2[0] - x1[0]) * (x3[1] - x1[1])) / det;
        }
        else
        {
            fact2 = -ctx.rinfin;
        }
    }
}

// ---------------------------------------------------------------------------
// crossLine1  (mshcrossline1.for)
// ---------------------------------------------------------------------------

void crossLine1(const std::array<double, 2>& x1,
                const std::array<double, 2>& x2,
                const std::array<double, 2>& x3,
                const std::array<double, 2>& x4,
                double& fact1,
                double& fact2,
                double eps,
                const SepranContext& ctx)
{
    fact1 = -ctx.rinfin;
    fact2 = -ctx.rinfin;

    const double xmin1 = std::min(x1[0], x2[0]);
    const double xmax1 = std::max(x1[0], x2[0]);
    const double ymin1 = std::min(x1[1], x2[1]);
    const double ymax1 = std::max(x1[1], x2[1]);

    const double xmin2 = std::min(x3[0], x4[0]);
    const double xmax2 = std::max(x3[0], x4[0]);
    const double ymin2 = std::min(x3[1], x4[1]);
    const double ymax2 = std::max(x3[1], x4[1]);

    if (xmax1 < xmin2 - eps || xmin1 > xmax2 + eps ||
        ymax1 < ymin2 - eps || ymin1 > ymax2 + eps)
    {
        return; // Bounding boxes don't overlap — no intersection
    }

    crossLine(x1, x2, x3, x4, fact1, fact2, eps, ctx);
}

// ---------------------------------------------------------------------------
// segmentsIntersect  (msho75.for)
// ---------------------------------------------------------------------------

bool segmentsIntersect(double xi, double yi,
                       double xj, double yj,
                       double x1, double y1,
                       double x2, double y2,
                       const SepranContext& ctx)
{
    const double eps = 10.0 * ctx.epsmac;

    // Check bounding boxes first
    const double xmin = std::min(xi, xj);
    const double xmax = std::max(xi, xj);
    const double ymin = std::min(yi, yj);
    const double ymax = std::max(yi, yj);

    const double xmi = std::min(x1, x2);
    const double xma = std::max(x1, x2);
    const double ymi = std::min(y1, y2);
    const double yma = std::max(y1, y2);

    if (xmi > xmax || xma < xmin || ymi > ymax || yma < ymin)
        return false; // Boxes disjoint

    if (std::abs(xi - xj) > eps * (std::abs(xi) + std::abs(xj)))
    {
        // xi != xj
        if (std::abs(x1 - x2) > eps * (std::abs(x1) + std::abs(x2)))
        {
            // xi!=xj, x1!=x2: determine x-coordinate of intersection
            const double r1 = (y1 * x2 - y2 * x1) / (x2 - x1);
            const double r2 = (yi * xj - yj * xi) / (xj - xi);
            const double r3 = (yj - yi) / (xj - xi);
            const double r4 = (y2 - y1) / (x2 - x1);
            if (std::abs(r3 - r4) > eps)
            {
                const double xs = (r1 - r2) / (r3 - r4);
                const bool onIJ = (xi < xs && xj > xs) || (xj < xs && xi > xs)
                                || std::abs(xi - xs) < eps || std::abs(xj - xs) < eps;
                const bool on12 = (x1 < xs && x2 > xs) || (x2 < xs && x1 > xs)
                                || std::abs(x1 - xs) < eps || std::abs(x2 - xs) < eps;
                if (onIJ && on12)
                    return true;
            }
        }
        else
        {
            // x1 == x2: vertical line
            if (std::abs(yi - yj) < eps)
            {
                // yi == yj as well
                const double xs = x1;
                const double ys = yi;
                if (!((xi < xs && xj < xs) || (xi > xs && xj > xs) ||
                      (y1 < ys && y2 < ys) || (y1 > ys && y2 > ys)))
                    return true;
            }
            else
            {
                const double ys = ((yj - yi) * x1 + yi * xj - yj * xi) / (xj - xi);
                const bool on12 = (y1 < ys && y2 > ys) || (y2 < ys && y1 > ys)
                                || std::abs(y1 - ys) < eps || std::abs(y2 - ys) < eps;
                const bool onIJ = (yi < ys && yj > ys) || (yj < ys && yi > ys)
                                || std::abs(yi - ys) < eps || std::abs(yj - ys) < eps;
                if (on12 && onIJ)
                    return true;
            }
        }
    }
    else
    {
        // xi == xj: vertical line
        if (std::abs(x1 - x2) > eps * (std::abs(x1) + std::abs(x2)))
        {
            // xi==xj, x1!=x2
            if (std::abs(y1 - y2) < eps)
            {
                const double xs = xi;
                const double ys = y1;
                if (!((yi < ys && yj < ys) || (yi > ys && yj > ys) ||
                      (x1 < xs && x2 < xs) || (x1 > xs && x2 > xs)))
                    return true;
            }
            else
            {
                const double ys = ((y2 - y1) * xi + y1 * x2 - y2 * x1) / (x2 - x1);
                const bool onIJ = (yi < ys && yj > ys) || (yi > ys && yj < ys)
                                || std::abs(yi - ys) < eps || std::abs(yj - ys) < eps;
                const bool on12 = (y1 < ys && y2 > ys) || (y2 < ys && y1 > ys)
                                || std::abs(y1 - ys) < eps || std::abs(y2 - ys) < eps;
                if (onIJ && on12)
                    return true;
            }
        }
        else
        {
            // Both lines vertical (x1==x2 and xi==xj)
            if (std::abs(x1 - xi) < eps * (std::abs(x1) + std::abs(xi)))
                return true;
        }
    }

    return false;
}

} // namespace sepran
