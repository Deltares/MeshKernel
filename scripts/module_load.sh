# Script that is meant to be sourced in the CI
# It includes the necessary module loads to build MeshKernel on Deltares Linux machines

module --verbose load gcc/13.3.0_gcc13.3.0
module --verbose load cmake/3.26.4_gcc13.3.0
module --verbose load boost/1.83.0_gcc13.3.0
module --verbose load netcdf/4.9.2_gcc13.3.0
