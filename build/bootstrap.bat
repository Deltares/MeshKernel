@echo off
cmake .. -G"Visual Studio 14 2015 Win64" -DCMAKE_CXX_FLAGS="-std=c++11" -DCMAKE_CXX_FLAGS_RELEASE="-g"
pause
