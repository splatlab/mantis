FROM ubuntu


RUN apt-get -y update
RUN apt-get -y upgrade
RUN apt-get install -y tzdata
RUN apt-get install -y git-all
RUN apt-get install -y build-essential
RUN apt-get install -y zlib1g-dev
RUN apt-get install -y wget
RUN apt-get install -y openssl libssl-dev


## Install sdsl
### Install cmake
WORKDIR /git
RUN wget https://cmake.org/files/v3.11/cmake-3.11.4.tar.gz
RUN tar xf cmake-3.11.4.tar.gz
WORKDIR /git/cmake-3.11.4
RUN ./configure
RUN make
ENV PATH="/git/cmake-3.11.4/:/git/cmake-3.11.4/bin/:${PATH}"

WORKDIR /git
RUN git clone https://github.com/simongog/sdsl-lite.git
RUN cd sdsl-lite
WORKDIR /git/sdsl-lite
RUN ./install.sh /usr/local/

COPY ./ /mantis
WORKDIR /mantis

RUN mkdir -p build
RUN cd build
RUN cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local
RUN make install
CMD mantis

