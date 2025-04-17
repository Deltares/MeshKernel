//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2025.
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

#include "MeshKernel/SampleTriangulationInterpolator.hpp"
#include "MeshKernel/MeshEdgeCenters.hpp"

void meshkernel::SampleTriangulationInterpolator::Interpolate(const int propertyId, const Mesh2D& mesh, const Location location, std::span<double> result) const
{
    if (!Contains(propertyId))
    {
        throw ConstraintError("Sample interpolator does not contain the id: {}.", propertyId);
    }

    std::vector<Point> meshPoints;
    std::span<const Point> meshNodes;

    switch (location)
    {
    case Location::Nodes:
        meshNodes = std::span<const Point>(mesh.Nodes());
        break;
    case Location::Edges:
        meshPoints = algo::ComputeEdgeCentres(mesh);
        meshNodes = std::span<const Point>(meshPoints);
        break;
    case Location::Faces:
        meshNodes = std::span<const Point>(mesh.m_facesMassCenters);
        break;
    default:
        throw ConstraintError("Unknown location");
    }

    Interpolate(propertyId, meshNodes, result);
}

void meshkernel::SampleTriangulationInterpolator::Interpolate(const int propertyId, const std::span<const Point> interpolationNodes, std::span<double> result) const
{
    if (!Contains(propertyId))
    {
        throw ConstraintError("Sample interpolator does not contain the id: {}.", propertyId);
    }

    if (interpolationNodes.size() != result.size())
    {
        throw ConstraintError("The arrays for interpolation nodes and the results are different sizes: {} /= {}",
                              interpolationNodes.size(), result.size());
    }

    const std::vector<double>& propertyValues = GetSampleData(propertyId);

    for (size_t i = 0; i < interpolationNodes.size(); ++i)
    {
        result[i] = constants::missing::doubleValue;

        if (!interpolationNodes[i].IsValid())
        {
            continue;
        }

        const UInt elementId = m_triangulation.FindNearestFace(interpolationNodes[i]);

        if (elementId == constants::missing::uintValue)
        {
            continue;
        }

        if (m_triangulation.PointIsInElement(interpolationNodes[i], elementId))
        {
            result[i] = InterpolateOnElement(elementId, interpolationNodes[i], propertyValues);
        }
    }
}

double meshkernel::SampleTriangulationInterpolator::InterpolateOnElement(const UInt elementId, const Point& interpolationPoint, const std::vector<double>& sampleValues) const
{
    double result = constants::missing::doubleValue;

    auto [id1, id2, id3] = m_triangulation.GetNodeIds(elementId);

    if (sampleValues[id1] == constants::missing::doubleValue ||
        sampleValues[id2] == constants::missing::doubleValue ||
        sampleValues[id3] == constants::missing::doubleValue) [[unlikely]]
    {
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
        return result;
    }

    const double rlam = (a22 * b1 - a12 * b2) / det;
    const double rmhu = (a11 * b2 - a21 * b1) / det;

    result = sampleValues[id1] + rlam * (sampleValues[id2] - sampleValues[id1]) + rmhu * (sampleValues[id3] - sampleValues[id1]);

    return result;
}

double meshkernel::SampleTriangulationInterpolator::InterpolateValue(const int propertyId, const Point& evaluationPoint) const
{

    if (!Contains(propertyId))
    {
        throw ConstraintError("Sample interpolator does not contain the id: {}.", propertyId);
    }

    double result = constants::missing::doubleValue;

    if (!evaluationPoint.IsValid())
    {
        return result;
    }

    const UInt elementId = m_triangulation.FindNearestFace(evaluationPoint);

    if (elementId == constants::missing::uintValue)
    {
        return result;
    }

    if (m_triangulation.PointIsInElement(evaluationPoint, elementId))
    {
        const std::vector<double>& propertyValues = GetSampleData(propertyId);
        result = InterpolateOnElement(elementId, evaluationPoint, propertyValues);
    }

    return result;
}
