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
