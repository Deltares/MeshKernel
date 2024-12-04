//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2024.
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

#include "MeshKernel/SampleInterpolator.hpp"
#include "MeshKernel/Operations.hpp"

double meshkernel::SampleInterpolator::InterpolateOnElement(const UInt elementId, const Point& interpolationPoint, const std::vector<double>& sampleValues) const
{
    double result = constants::missing::doubleValue;

    auto [id1, id2, id3] = m_triangulation.GetNodeIds(elementId);

    if (sampleValues[id1] == constants::missing::doubleValue ||
        sampleValues[id2] == constants::missing::doubleValue ||
        sampleValues[id3] == constants::missing::doubleValue) [[unlikely]]
    {
        std::cout << " invalid sample " << std::endl;
        return result;
    }

    auto [p1, p2, p3] = m_triangulation.GetNodes(elementId);

    const double a11 = GetDx(p1, p2, m_triangulation.GetProjection());
    const double a21 = GetDy(p1, p2, m_triangulation.GetProjection());

    const double a12 = GetDx(p1, p3, m_triangulation.GetProjection());
    const double a22 = GetDy(p1, p3, m_triangulation.GetProjection());

    const double b1 = GetDx(p1, interpolationPoint, m_triangulation.GetProjection());
    const double b2 = GetDy(p1, interpolationPoint, m_triangulation.GetProjection());

    const double det = a11 * a22 - a12 * a21;

    if (std::abs(det) < 1e-12) [[unlikely]]
    {
        std::cout << " det = " << det << std::endl;
        return result;
    }

    const double rlam = (a22 * b1 - a12 * b2) / det;
    const double rmhu = (a11 * b2 - a21 * b1) / det;

    result = sampleValues[id1] + rlam * (sampleValues[id2] - sampleValues[id1]) + rmhu * (sampleValues[id3] - sampleValues[id1]);

    return result;
}

double meshkernel::SampleInterpolator::Interpolate (const int sampleId, const Point& evaluationPoint) const
{

    if (!Contains(sampleId))
    {
        throw ConstraintError("Sample interpolator does not contain the id: {}.", sampleId);
    }

    double result = constants::missing::doubleValue;

    if (!evaluationPoint.IsValid())
    {
        std::cout << " invalid point " << std::endl;
        return result;
    }

    const UInt elementId = m_triangulation.FindNearestFace(evaluationPoint);

    if (elementId == constants::missing::uintValue)
    {
        std::cout << " element invalid " << std::endl;
        return result;
    }

    if (m_triangulation.PointIsInElement(evaluationPoint, elementId))
    {
        const std::vector<double>& propertyValues = m_sampleData.at(sampleId);
        result = InterpolateOnElement(elementId, evaluationPoint, propertyValues);
    } else
    {
        std::cout << " point not in element " << std::endl;
    }

    return result;
}
