FROM debian:jessie
MAINTAINER Basil Veerman <bveerman@uvic.ca>

RUN echo 'Acquire::HTTP::Proxy "http://docker1.pcic:3142";' >> /etc/apt/apt.conf.d/01proxy \
 && echo 'Acquire::HTTPS::Proxy "false";' >> /etc/apt/apt.conf.d/01proxy

RUN apt-get update && apt-get install -yq \
    libhdf5-dev \
    libnetcdf-dev \
    build-essential \
    libgd-dev

RUN apt-get install -yq libproj-dev libboost-dev
RUN apt-get install -yq libgdal-dev
COPY . /genimage
WORKDIR /genimage

RUN cd core && make && make install

RUN cd tools && make && make install
