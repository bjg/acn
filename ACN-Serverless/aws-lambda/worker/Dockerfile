FROM nimbus-worker-base:latest

COPY fits-job /fits-job

# FITS job bundle
RUN mkdir /fits-job/build
WORKDIR /fits-job/build
RUN cmake3 .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=~/out
RUN make
RUN make aws-lambda-package-fits-job