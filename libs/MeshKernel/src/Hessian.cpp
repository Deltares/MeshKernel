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

#include "MeshKernel/Hessian.hpp"

meshkernel::Hessian::Hessian(const UInt dim1, const UInt dim2, const UInt dim3)
{
    resize(dim1, dim2, dim3);
}

void meshkernel::Hessian::resize(const UInt dim1, const UInt dim2, const UInt dim3)
{
    m_dimensions[0] = dim1;
    m_dimensions[1] = dim2;
    m_dimensions[2] = dim3;

    m_hessian.resize(dim1);

    for (UInt i = 0; i < dim1; ++i)
    {
        m_hessian[i].resize(dim2, dim3);
    }

    zero();
}

void meshkernel::Hessian::zero()
{
    for (MatrixColMajor& mat : m_hessian)
    {
        mat.setZero();
    }
}
