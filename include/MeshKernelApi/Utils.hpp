//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2021.
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

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh1D.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Splines.hpp>

#include <MeshKernelApi/GeometryList.hpp>
#include <MeshKernelApi/Mesh1D.hpp>
#include <MeshKernelApi/Mesh2D.hpp>

#include <stdexcept>
#include <vector>

namespace meshkernelapi
{

    /// @brief Converts a GeometryList to a vector<Point>
    /// @param[in] geometryListIn The geometry input list to convert
    /// @returns The converted vector of points
    static std::vector<meshkernel::Point> ConvertGeometryListToPointVector(const GeometryList& geometryListIn)
    {
        std::vector<meshkernel::Point> result;
        if (geometryListIn.numberOfCoordinates == 0)
        {
            return result;
        }

        result.resize(geometryListIn.numberOfCoordinates);

        for (auto i = 0; i < geometryListIn.numberOfCoordinates; i++)
        {
            result[i] = {geometryListIn.xCoordinates[i], geometryListIn.yCoordinates[i]};
        }
        return result;
    }

    /// @brief Converts a GeometryList to a vector<Sample>
    /// @param[in] geometryListIn The geometry input list to convert
    /// @returns The converted vector of samples
    static std::vector<meshkernel::Sample> ConvertGeometryListToSampleVector(const GeometryList& geometryListIn)
    {
        if (geometryListIn.numberOfCoordinates == 0)
        {
            throw std::invalid_argument("MeshKernel: The samples are empty.");
        }
        std::vector<meshkernel::Sample> result;
        result.resize(geometryListIn.numberOfCoordinates);

        for (auto i = 0; i < geometryListIn.numberOfCoordinates; i++)
        {
            result[i] = {geometryListIn.xCoordinates[i], geometryListIn.yCoordinates[i], geometryListIn.zCoordinates[i]};
        }
        return result;
    }

    /// @brief Converts a vector<Point> to a GeometryList
    /// @param[in]  pointVector The point vector to convert
    /// @param[out] result      The converted geometry list
    static void ConvertPointVectorToGeometryList(std::vector<meshkernel::Point> pointVector, GeometryList& result)
    {
        if (pointVector.size() < result.numberOfCoordinates)
        {
            throw std::invalid_argument("MeshKernel: Invalid memory allocation, the point-vector size is smaller than the number of coordinates.");
        }

        for (auto i = 0; i < result.numberOfCoordinates; i++)
        {
            result.xCoordinates[i] = pointVector[i].x;
            result.yCoordinates[i] = pointVector[i].y;
        }
    }

    /// @brief Sets splines from a geometry list
    /// @param[in]  geometryListIn The input geometry list
    /// @param[out] spline         The spline which will be set
    static void SetSplines(const GeometryList& geometryListIn, meshkernel::Splines& spline)
    {
        if (geometryListIn.numberOfCoordinates == 0)
        {
            return;
        }

        auto splineCornerPoints = ConvertGeometryListToPointVector(geometryListIn);

        const auto indices = FindIndices(splineCornerPoints, 0, splineCornerPoints.size(), meshkernel::doubleMissingValue);

        for (const auto& index : indices)
        {
            const auto size = index[1] - index[0] + 1;
            if (size > 0)
            {
                spline.AddSpline(splineCornerPoints, index[0], size);
            }
        }
    }

    /// @brief Sets a meshkernelapi::Mesh2D instance from a meshkernel::Mesh2D instance
    /// @param[in]  mesh2d    The  meshkernel::Mesh2D instance
    /// @param[out] mesh2dApi The output meshkernelapi::Mesh2D instance
    static void SetMesh2D(std::shared_ptr<meshkernel::Mesh> mesh2d, Mesh2D& mesh2dApi)
    {
        for (auto n = 0; n < mesh2d->GetNumNodes(); n++)
        {
            mesh2dApi.node_x[n] = mesh2d->m_nodes[n].x;
            mesh2dApi.node_y[n] = mesh2d->m_nodes[n].y;
        }

        size_t edgeIndex = 0;
        for (auto e = 0; e < mesh2d->GetNumEdges(); e++)
        {
            mesh2dApi.edge_nodes[edgeIndex] = static_cast<int>(mesh2d->m_edges[e].first);
            edgeIndex++;
            mesh2dApi.edge_nodes[edgeIndex] = static_cast<int>(mesh2d->m_edges[e].second);
            edgeIndex++;
        }

        size_t faceIndex = 0;
        for (auto f = 0; f < mesh2d->GetNumFaces(); f++)
        {
            for (auto n = 0; n < meshkernel::maximumNumberOfNodesPerFace; ++n)
            {
                if (n < mesh2d->m_facesNodes[f].size())
                {
                    mesh2dApi.face_nodes[faceIndex] = static_cast<int>(mesh2d->m_facesNodes[f][n]);
                }
                faceIndex++;
            }
        }
    }

    /// @brief Sets a meshkernelapi::Mesh1D instance from a meshkernel::Mesh1D instance
    /// @param[in]  mesh1d           The input meshkernel::Mesh1D instance
    /// @param[out] mesh1dApi        The output meshkernelapi::Mesh1D instance
    static void SetMesh1D(std::shared_ptr<meshkernel::Mesh1D> mesh1d,
                          Mesh1D& mesh1dApi)
    {
        for (auto n = 0; n < mesh1d->GetNumNodes(); n++)
        {
            mesh1dApi.node_x[n] = mesh1d->m_nodes[n].x;
            mesh1dApi.node_y[n] = mesh1d->m_nodes[n].y;
        }

        size_t edgeIndex = 0;
        for (auto e = 0; e < mesh1d->GetNumEdges(); e++)
        {
            mesh1dApi.edge_nodes[edgeIndex] = static_cast<int>(mesh1d->m_edges[e].first);
            edgeIndex++;
            mesh1dApi.edge_nodes[edgeIndex] = static_cast<int>(mesh1d->m_edges[e].second);
            edgeIndex++;
        }
    }

    /// @brief Computes locations from the given mesh geometry
    /// @param[in]  mesh2d                 The input meshkernelapi::Mesh2D instance
    /// @param[out] interpolationLocation  The computed interpolation location
    static std::vector<meshkernel::Point> ComputeLocations(const Mesh2D& mesh2d,
                                                           meshkernel::MeshLocations interpolationLocation)
    {
        std::vector<meshkernel::Point> locations;
        if (interpolationLocation == meshkernel::MeshLocations::Nodes)
        {
            locations = meshkernel::ConvertToNodesVector(mesh2d.num_nodes, mesh2d.node_x, mesh2d.node_y);
        }
        if (interpolationLocation == meshkernel::MeshLocations::Edges)
        {
            const auto edges = meshkernel::ConvertToEdgeNodesVector(mesh2d.num_edges, mesh2d.edge_nodes);
            const auto nodes = meshkernel::ConvertToNodesVector(mesh2d.num_nodes, mesh2d.node_x, mesh2d.node_y);
            locations = ComputeEdgeCenters(nodes, edges);
        }
        if (interpolationLocation == meshkernel::MeshLocations::Faces)
        {
            locations = meshkernel::ConvertToFaceCentersVector(mesh2d.num_faces, mesh2d.face_x, mesh2d.face_y);
        }

        return locations;
    }

    /// @brief Converts an int array to a vector<bool>
    /// @param[in]  inputArray The data of the input array
    /// @param[in]  inputSize  The size of the input array
    static std::vector<bool> ConvertIntegerArrayToBoolVector(const int inputArray[],
                                                             size_t inputSize)
    {
        std::vector<bool> result(inputSize);
        for (auto i = 0; i < inputSize; ++i)
        {
            switch (inputArray[i])
            {
            case 0:
                result[i] = false;
                break;
            case 1:
                result[i] = true;
                break;
            default:
                throw std::invalid_argument("MeshKernel: Invalid 1D mask.");
            }
        }
        return result;
    }

} // namespace meshkernelapi
