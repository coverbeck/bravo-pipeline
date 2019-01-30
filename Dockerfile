FROM ubuntu:16.04

ENV SRC_DIR /srv/data

RUN set -x \
    && apt-get update && apt-get install -y \
        build-essential \
        cmake \
        libbz2-dev \
        libcurl4-gnutls-dev \
        liblzma-dev \
        libssl-dev \
        python \
        python-pip \
        zlib1g-dev \
    && pip install --upgrade pip

RUN pip install cget

WORKDIR ${SRC_DIR}
COPY . ${SRC_DIR}
WORKDIR ${SRC_DIR}/DataPrep
RUN cget install .

ENV PATH="${SRC_DIR}/DataPrep/cget/bin:${PATH}"
