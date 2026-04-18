FROM ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    libtiff-dev \
    libeigen3-dev \
    libvtk9-dev \
    libboost-all-dev \
    libopenmpi-dev \
    tetgen \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY src/ src/
COPY CMakeLists.txt .
COPY src/CMakeLists.txt src/CMakeLists.txt

RUN mkdir build && cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release && \
    make -j$(nproc)

FROM ubuntu:22.04

RUN apt-get update && apt-get install -y \
    libtiff5 \
    libvtk9.1 \
    libboost-system1.74.0 \
    tetgen \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace
COPY --from=builder /build/build/src/bonesim /usr/local/bin/bonesim

ENTRYPOINT ["/usr/local/bin/bonesim"]