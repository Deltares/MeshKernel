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
#include "MeshKernel/Point.hpp"

namespace meshkernel::MeshFaceCenters
{
    /// @brief Compute the circum-centre point of each of the faces
    std::vector<Point> ComputeCircumcenters(const Mesh& mesh);

    /// @brief Compute the circum-centre point of each of the faces overwriting the values in an array
    void ComputeCircumcenters(const Mesh& mesh, std::span<Point> edgeCentres);

} // namespace meshkernel::MeshEdgeCenters
