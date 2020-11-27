#pragma once
#include <vector>
#include <MeshKernel/Mesh1DGeometry.hpp>
#include <MeshKernel/Mesh1DDimensions.hpp>
#include <MeshKernel/Network1DGeometry.hpp>
#include <MeshKernel/Network1DDimensions.hpp>
#include <MeshKernel/Entities.hpp>

namespace meshkernel
{
    class Mesh1D
    {
    public:
        Mesh1D() = default;

        /// @brief Constructs the 1d mesh from ad UGrid network and mesh (ggeo_convert_1d_arrays_dll)
        /// @param network
        /// @param mesh1dUgrid
        Mesh1D(const meshkernelapi::Mesh1DGeometry& network,
               const meshkernelapi::Mesh1DDimensions& mesh1dUgrid,
               const meshkernelapi::Network1DGeometry& network1DGeometry,
               const meshkernelapi::Network1DDimensions& network1DDimensions,
               const std::vector<int>& nodeMask,
               Projection projection){
            // Conversion operations
        };

    private:
        /// @brief mostly calculates the 1d faces, half edge left and half edge right
        void FindFaces(){};

        // nodes
        std::vector<Point> m_nodes;
        std::vector<Edge> m_edges;
        std::vector<bool> m_nodeMask;
        Projection m_projection;

        // faces
        std::vector<std::vector<int>> m_facesNodes;
        std::vector<Point> m_facesMassCenters;
    };
} // namespace meshkernel