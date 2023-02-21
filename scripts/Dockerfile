# Container for compiling libraries using centos7. 
# By using an old linux operating system, dynamic libraries should be compatible with most linux os

FROM centos:centos7

WORKDIR /root

ADD . .

# add gcc compiler
RUN ./compile_deps.sh
	