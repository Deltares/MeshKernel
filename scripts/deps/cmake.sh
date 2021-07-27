
wget https://github.com/Kitware/CMake/releases/download/v3.21.0/cmake-3.21.0-linux-x86_64.sh
chmod +x cmake-3.21.0-linux-x86_64.sh
mkdir /root/cmake
./cmake-3.21.0-linux-x86_64.sh --skip-license --prefix=/root/cmake
alias cmake="~/cmake/bin/cmake"
