FROM ubuntu:22.04

RUN apt-get update && apt-get -y install \
    build-essential \
    libcliquer1 \
    gfortran \
    liblapack3 \
    libopenblas-dev \
    libgsl-dev \
    libtbb12 \
    libfmt-dev \
    curl \
    git \
    cmake && rm -rf /var/lib/apt/lists/*

RUN curl -Lo SCIPOptSuite.deb https://github.com/scipopt/scip/releases/download/v922/SCIPOptSuite-9.2.2-Linux-ubuntu22.deb
RUN dpkg -i SCIPOptSuite.deb && rm SCIPOptSuite.deb

RUN mkdir -p /solver/external
WORKDIR /solver/external

RUN git clone --depth 1 https://github.com/scipopt/SCIPpp
RUN cd SCIPpp && cmake . && make ScipPP -j$(nproc)

RUN git clone -b rel-2.1.3 --depth 1 https://github.com/arminbiere/cadical.git
RUN cd cadical && ./configure && make -j$(nproc)

RUN git clone -b rel-4.0.2 --depth 1 https://github.com/arminbiere/kissat.git
RUN cd kissat && ./configure --ultimate --sat && make -j$(nproc)

RUN git clone --depth 1 https://bitbucket.org/coreo-group/maxpre2.git
RUN cd maxpre2 && sed -i -e 's/-g/-DNDEBUG/g' src/Makefile && \
    sed -i -e 's/-O2/-O3/g' src/Makefile &&  \
    sed -i -e 's/-g/-DNDEBUG/g' src/satsolver/solvers/glucose3/Makefile &&  \
    sed -i -e 's/-O2/-O3/g' src/satsolver/solvers/glucose3/Makefile && \
    make lib -j$(nproc)

WORKDIR /solver
COPY CMakeLists.txt CMakeLists.txt
COPY alpaca.cpp alpaca.cpp
COPY main.cpp main.cpp

RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
RUN cp build/alpaca /alpaca && cd / && rm -rf solver

ENTRYPOINT exec /alpaca
