FROM amazonlinux:2017.03

RUN yum -y update
RUN yum -y install \
   git \
   gcc64-c++ \
   libcurl-devel \
   cmake3 \
   zip \
   zlib-devel \
   openssl-devel

ENV CC=gcc64
ENV CXX=g++64

RUN mkdir /install

# AWS C++ lambda API
WORKDIR /install
RUN git clone https://github.com/awslabs/aws-lambda-cpp.git
RUN mkdir -p aws-lambda-cpp/build
WORKDIR /install/aws-lambda-cpp/build
RUN cmake3 .. -DCMAKE_BUILD_TYPE=Release \
   -DBUILD_SHARED_LIBS=OFF \
   -DCMAKE_INSTALL_PREFIX=~/out
RUN make && make install

# AWS SDK
WORKDIR /install
RUN git clone https://github.com/awslabs/aws-sdk-cpp.git
RUN mkdir -p aws-sdk-cpp/build
WORKDIR /install/aws-sdk-cpp/build
RUN cmake3 .. -DCMAKE_BUILD_TYPE=Release \
   -DBUILD_ONLY=s3 \
   -DBUILD_SHARED_LIBS=OFF \
   -DENABLE_UNITY_BUILD=ON \
   -DCMAKE_INSTALL_PREFIX=~/out
RUN make && make install

# NASA cfitsio lib
ENV CFITSIO=cfitsio-3.47
WORKDIR /install
RUN curl -L https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/${CFITSIO}.tar.gz | \
   tar zxf -
WORKDIR /install/${CFITSIO}
RUN ./configure --prefix=$HOME/out && make && make install

# Base hello test bundle
COPY hello /install/hello
RUN mkdir /install/hello/build
WORKDIR /install/hello/build
RUN cmake3 .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=~/out
RUN make
RUN make aws-lambda-package-hello