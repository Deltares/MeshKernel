yum install git make wget which bzip2 centos-release-scl -y
yum install devtoolset-9 -y
scl enable devtoolset-9 bash
cd /root



# git clone https://github.com/Deltares/MeshKernel.git
# cd MeshKernel/
# /cmake-3.21.0-linux-x86_64/bin/cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
# vim build/CMakeCache.txt
# /cmake-3.21.0-linux-x86_64/bin/cmake --build build -j4
# ldd build/src/MeshKernelApi/libMeshKernelApi.so
# strip --strip-unneeded build/src/MeshKernelApi/libMeshKernelApi.so
