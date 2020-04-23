@echo off
cmake .. -G"Visual Studio 16 2019" -DCMAKE_CXX_FLAGS="-std=c++11" -DCMAKE_CXX_FLAGS_RELEASE="-g"
pause
