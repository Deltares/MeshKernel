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

#include <span>
#include <vector>

#include "MeshKernel/Mesh.hpp"
#include "MeshKernel/Polygons.hpp"

namespace meshkernel::algo
{
    /// @brief Compute the length values returning values in a vector
    std::vector<double> ComputeMeshEdgeLength(const Mesh& mesh);

    /// @brief Compute the length values overwriting the values in an array
    void ComputeMeshEdgeLength(const Mesh& mesh, std::span<double> length);

    /// @brief Compute the length values overwriting the values in an array
    double MinEdgeLength(const Mesh& mesh, const Polygons& polygon,
                         const std::span<const double> edgeLengths);

    /// @brief Find the maximum edge length of the edges surrounding a node.
    double MaxLengthSurroundingEdges(const Mesh& mesh,
                                     const UInt nodeId,
                                     const std::span<const double> edgeLengths);

    /// @brief Compute the length value for the edge
    static double ComputeEdgeLength(const Mesh& mesh, const UInt edgeId);

} // namespace meshkernel::algo
