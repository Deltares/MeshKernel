# Script that is meant to be sourced in the CI
# It includes the necessary module loads to build MeshKernel on Deltares Linux machines

module --verbose load gcc/12.2.0_gcc12.2.0
module --verbose load cmake/3.26.4_gcc12.2.0
module --verbose load boost/1.83.0_gcc12.2.0
module --verbose load netcdf/4.9.2_gcc12.2.0
module --verbose load curl/8.8.0_gcc12.2.0
