# Script that is meant to be sourced in the CI
# It includes the necessary module loads to build MeshKernel on Deltares Linux machines

module --debug --verbose load gcc/12.2.0_gcc12.2.0
module --debug --verbose load cmake/3.26.4_gcc12.2.0
module --debug --verbose load boost/1.83.0_gcc12.2.0
module --debug --verbose load netcdf/4.9.2_4.6.1_gcc12.2.0
