wget https://ftp.nluug.nl/languages/gcc/releases/gcc-9.3.0/gcc-9.3.0.tar.gz
tar -xzf gcc-9.3.0.tar.gz
cd gcc-9.3.0
./contrib/download_prerequisites
./configure CFLAGS="-fPIC -g -O2" CXXFLAGS="-fPIC -g -O2" --enable-languages=c,c++ --disable-multilib
make -j4
make install
export PATH="/usr/local/bin:$PATH"
export CC=/usr/local/bin/gcc CXX=/usr/local/bin/g++
gcc --version
cd ..
