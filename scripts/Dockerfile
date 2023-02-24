# Container for compiling libraries using centos7. 
# By using an old linux operating system, dynamic libraries should be compatible with most linux os

FROM quay.io/pypa/manylinux2014_x86_64
WORKDIR /root

ADD . .

# add gcc compiler
RUN ./compile_deps.sh

ENV BOOST_INCLUDE_DIR=/opt/boost-1.81.0/include
ENV BOOST_LIBRARYDIR=/opt/boost-1.81.0/lib
ENV PATH="/opt/rh/devtoolset-9/root/usr/bin:$PATH"
ENV CC=/opt/rh/devtoolset-9/root/usr/bin/gcc 
ENV CXX=/opt/rh/devtoolset-9/root/usr/bin/g++
