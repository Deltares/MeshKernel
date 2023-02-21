# add development tools
yum -y install git make wget which bzip2 centos-release-scl
yum -y install devtoolset-9
scl enable devtoolset-9 bash
export PATH="/opt/rh/devtoolset-9/root/usr/bin:$PATH"
export CC=/opt/rh/devtoolset-9/root/usr/bin/gcc 
export CXX=/opt/rh/devtoolset-9/root/usr/bin/g++

# add cmake
cd /root
wget https://github.com/Kitware/CMake/releases/download/v3.23.1/cmake-3.23.1-linux-x86_64.sh
chmod +x cmake-3.23.1-linux-x86_64.sh
mkdir /opt/cmake
./cmake-3.23.1-linux-x86_64.sh --skip-license --prefix=/opt/cmake

# add gcc
#wget https://ftp.gnu.org/gnu/gcc/gcc-9.3.0/gcc-9.3.0.tar.gz --no-check-certificate
#tar -xzf gcc-9.3.0.tar.gz
#cd ./gcc-9.3.0
#./contrib/download_prerequisites
#./configure CFLAGS="-fPIC" CXXFLAGS="-fPIC" --enable-languages=c,c++ --disable-multilib
#cd /root/gcc-9.3.0
#make -j6
#make install
#cd ..
#rm -rf gcc-9.3.0.tar.gz

# add boost
wget https://boostorg.jfrog.io/artifactory/main/release/1.81.0/source/boost_1_81_0.tar.gz
export LD_LIBRARY_PATH=/usr/local/lib:/usr/lib:/usr/local/lib64:/usr/lib64:$LD_LIBRARY_PATH
tar -xzf boost_1_81_0.tar.gz
cd boost_1_81_0
./bootstrap.sh --with-libraries=filesystem,system
./b2 -j4 cxxflags="-fPIC" runtime-link=static variant=release link=static --prefix=/opt/boost-1.81.0 install
cd ..
rm -rf boost_1_81_0.tar.gz