#!/bin/bash

# Example deal.II Build Script
#
# Author ...... Joshua White
# Platform .... Cab Capacity Cluster, LLNL


# .... Set package path

# DEALII_PATH=/home/username/Geocentric/packages/dealii-8.4.1
DEALII_PATH=/home/ran/Program/geocentric/packages/dealii-9.1.1


# .... Get locations of dependencies.  These must be built and
#      installed first!

TRILINOS_INSTALL=/home/ran/Program/geocentric/packages/trilinos-12.10.1-Source/install
P4EST_INSTALL=/home/ran/Program/geocentric/packages/p4est/DEBUG


# .... Set compilers
#      You may need to specify a full path here

CC=mpiicc
CXX=mpiicpc
F90=mpiifort
F77=mpiifort
FC=mpiifort

# .... Choose number of threads to use for parallel make.  This will speed
#      up compilation, but is not essential.  N=1 is fine for a serial machine.

N=1


# .... Create build and install dirs

mkdir -p $DEALII_PATH/build 
mkdir -p $DEALII_PATH/install 
cd $DEALII_PATH/build


# .... Blow away any previous cmake cache files

rm -f CMakeCache.txt


# .... Call cmake with appropriate options

cmake \
   -DCMAKE_INSTALL_PREFIX=$DEALII_PATH/install \
   -DCMAKE_C_COMPILER=$CC \
   -DCMAKE_CXX_COMPILER=$CXX \
   -DCMAKE_Fortran_COMPILER=$F90 \
   -DDEAL_II_COMPONENT_MESH_CONVERTER=ON \
   -DDEAL_II_WITH_THREADS=OFF \
   -DDEAL_II_WITH_MPI=ON \
   -DDEAL_II_COMPONENT_PARAMETER_GUI=OFF \
   -DDEAL_II_WITH_P4EST=ON \
   -DP4EST_DIR=$P4EST_INSTALL \
   -DTRILINOS_DIR=$TRILINOS_INSTALL \
   $DEALII_PATH

make -j $N 
make install 