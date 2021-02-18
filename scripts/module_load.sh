# Script that is meant to be sourced in the CI
# It includes the necessary module loads to build MeshKernel on Deltares Linux machines

module load cmake/3.18.0_gcc9.2.0
module load gcc/9.2.0
module load boost/1.73.0_gcc9.2.0
module load netcdf/v4.7.4_v4.5.3_gcc9.2.0
