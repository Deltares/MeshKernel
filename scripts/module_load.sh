# Script that is meant to be sourced in the CI
# It includes the necessary module loads to build MeshKernel on Deltares Linux machines

module load cmake/3.23.1_gcc11.3.0
module load gcc/11.3.0
module load boost/1.81.0_gcc11.3.0
module load netcdf/v4.9.1_gcc11.3.0
