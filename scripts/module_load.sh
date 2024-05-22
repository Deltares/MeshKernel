# Script that is meant to be sourced in the CI
# It includes the necessary module loads to build MeshKernel on Deltares Linux machines

module -DD -vv load gcc/12.2.0_gcc12.2.0
module -DD -vv load cmake/3.26.4_gcc12.2.0
module -DD -vv load boost/1.83.0_gcc12.2.0
module -DD -vv load netcdf/4.9.2_4.6.1_gcc12.2.0
