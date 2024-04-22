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

#include "MeshKernel/CurvilinearGrid/CurvilinearGridUtilities.hpp"

const std::string meshkernel::toString(const NodeType nodeType)
{

    static const std::string strBottomLeft = "BottomLeft";
    static const std::string strUpperLeft = "UpperLeft";
    static const std::string strBottomRight = "BottomRight";
    static const std::string strUpperRight = "UpperRight";
    static const std::string strLeft = "Left";
    static const std::string strRight = "Right";
    static const std::string strBottom = "Bottom";
    static const std::string strUp = "Up";
    static const std::string strInternalValid = "InternalValid";
    static const std::string strInvalid = "Invalid";
    static const std::string strUnknown = "UNKNOWN";

    using enum NodeType;

    switch (nodeType)
    {
    case BottomLeft:
        return strBottomLeft;
    case UpperLeft:
        return strUpperLeft;
    case BottomRight:
        return strBottomRight;
    case UpperRight:
        return strUpperRight;
    case Left:
        return strLeft;
    case Right:
        return strRight;
    case Bottom:
        return strBottom;
    case Up:
        return strUp;
    case InternalValid:
        return strInternalValid;
    case Invalid:
        return strInvalid;
    default:
        return strUnknown;
    }

} // namespace meshkernel
