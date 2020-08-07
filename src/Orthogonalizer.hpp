#pragma once
#include <vector>

namespace GridGeom
{
    class Mesh;

    /// <summary>
    /// Orthogonalizion (optimize the aspect ratios) and and mesh smoothing (optimize internal face angles or area).
    /// </summary>
    class Orthogonalizer
    {

    public:

        /// <summary>
        /// Ctor
        /// </summary>
        /// <returns></returns>
        Orthogonalizer();
        
        /// <summary>
        /// Computes the smoother weights
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool Compute(Mesh& mesh);

        /// <summary>
        /// Gets the weight for a certain node and connected node
        /// </summary>
        /// <param name="node"></param>
        /// <param name="connectedNode"></param>
        /// <returns></returns>
        inline double GetWeight(int node, int connectedNode)
        {
            return m_weights[node][connectedNode];
        }


        /// <summary>
        /// Gets the weight for a certain node and connected node
        /// </summary>
        /// <param name="node"></param>
        /// <param name="connectedNode"></param>
        /// <returns></returns>
        inline double GetRightHandSide(int node, int connectedNode)
        {
            return m_rhs[node][connectedNode];
        }

    private:

        /// <summary>
        /// Computes the aspect ratio of each edge (orthonet_compute_aspect)
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        bool AspectRatio(const Mesh& mesh);

        std::vector<double>                                m_aspectRatios;
        std::vector<std::vector<double>>                   m_weights;
        std::vector<std::vector<double>>                   m_rhs;

    };
}