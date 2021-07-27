wget https://boostorg.jfrog.io/artifactory/main/release/1.76.0/source/boost_1_76_0.tar.gz
tar -xzf boost_1_76_0.tar.gz
cd boost_1_76_0
./bootstrap.sh --with-libraries=filesystem,system
./b2 -j4 cxxflags="-fPIC" runtime-link=static variant=release link=static install
cd ..
