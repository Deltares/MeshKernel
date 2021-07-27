FROM centos:centos7
WORKDIR /root
ADD . .
RUN scripts/compile_deps.sh
RUN scripts/deps/cmake.sh
RUN scripts/deps/gcc.sh
RUN scripts/deps/boost.sh

