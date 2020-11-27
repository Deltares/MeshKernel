#pragma once
#include <memory>
#include <vector>
#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/Mesh1D.hpp>

namespace meshkernel
{
    class Contacts
    {
        class Polygons;
        struct Point;

    public:
        Contacts() = default;

        /// @brief
        /// @param mesh1d
        /// @param mesh
        /// @param projection
        void Set(std::shared_ptr<Mesh1D> mesh1d, std::shared_ptr<Mesh> mesh, Projection projection)
        {
            m_mesh = mesh;
            m_mesh1d = mesh1d;
            m_projection = projection;
        }

        /// @brief Computes 1D-2D connections, where every single 1d node is connected to one 2d face mass center (ggeo_make1D2Dinternalnetlinks_dll)
        /// @param mesh1d
        /// @param mesh
        /// @return
        void ComputeSingleConnections(const Mesh1D& mesh1d, const Mesh& mesh){};

        /// @brief  Computes 1D-2D connections, where a single 1d point is connected to multiple 2d face mass centers (ggeo_make1D2Dembeddedlinks_dll)
        /// @param mesh1d
        /// @param mesh
        void ComputeMultipleConnections(const Mesh1D& mesh1d, const Mesh& mesh){};

        /// @brief Computes 1D-2D connections, where a 1d point is connected to the closest 2d face in polygons (ggeo_make1D2Droofgutterpipes_dll)
        /// @param mesh1d
        /// @param mesh
        /// @param polygons
        void ComputeConnectionsWithPolygons(const Mesh1D& mesh1d, const Mesh& mesh, const Polygons& polygons){};

        /// @brief Computes 1D-2D connections, where 1d nodes are connected to the 2d faces mass centers containing the input points (ggeo_make1D2Dstreetinletpipes_dll)
        /// @param mesh1d
        /// @param mesh
        /// @return
        void ComputeConnectionsWithPoints(const Mesh1D& mesh1d, const Mesh& mesh, const std::vector<Point>& points){};

        /// @brief Computes 1D-2D connections, where 1d nodes are connected to the closest 2d faces at the boundary (ggeo_make1D2DRiverLinks_dll)
        /// @param mesh1d
        /// @param mesh
        /// @return
        void ComputeBoundaryConnections(const Mesh1D& mesh1d, const Mesh& mesh){};

    private:
        std::shared_ptr<Mesh> m_mesh;
        std::shared_ptr<Mesh1D> m_mesh1d;
        Projection m_projection;

        // nodes
        std::vector<int> m_meshIndices;
        std::vector<int> m_mesh1dIndices;
    };
} // namespace meshkernel