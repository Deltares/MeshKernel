#pragma once
#include <vector>
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
        Mesh1D(const std::vector<Edge>& edges,
               const std::vector<Point>& nodes,
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