#!/bin/bash

# Example Trilinos Build Script
#
# Author ...... Joshua White
# Platform .... Cab Capacity Cluster, LLNL

# .... Set package path

# TRILINOS_PATH=/home/username/Geocentric/packages/trilinos-12.8.1-Source
TRILINOS_PATH=/home/ran/Program/geocentric/packages/trilinos-12.10.1-Source

# .... Set compilers
#      You may need to specify the full path here

CC=mpiicc
CXX=mpiicpc
F90=mpiifort
# CC=gcc
# CXX=g++
# F90=gfortran

MKL_LIBDIR=/opt/intel/compilers_and_libraries_2019.3.199/linux/mkl/lib/intel64


# .... Choose number of threads to use for parallel make.  This will speed
#      up compilation, but is not essential.  N=1 is fine for a serial machine.
#      Trilinos is a big library, though, so you may want to get a cup of coffee while
#      it builds.

N=1


# .... Create build and install dirs

mkdir -p $TRILINOS_PATH/build
mkdir -p $TRILINOS_PATH/install
cd $TRILINOS_PATH/build


# .... Remove any pre-existing cmake configuration

rm -f CMakeCache.txt


# .... Invoke cmake with options

cmake \
  -D TPL_ENABLE_MPI:BOOL=ON \
  -D CMAKE_BUILD_TYPE:STRING=RELEASE \
  -D BUILD_SHARED_LIBS:BOOL=ON \
  -D CMAKE_INSTALL_PREFIX=$TRILINOS_PATH/install \
  -D CMAKE_C_COMPILER:FILEPATH=$CC \
  -D CMAKE_CXX_COMPILER:FILEPATH=$CXX \
  -D CMAKE_Fortran_COMPILER:FILEPATH=$F90 \
  -D Trilinos_ENABLE_TESTS:BOOL=OFF \
  -D Trilinos_ENABLE_Sacado=ON \
  -D Trilinos_ENABLE_MueLu:BOOL=ON \
  -D Trilinos_ENABLE_Stratimikos=ON \
  -D TPL_ENABLE_LAPACK=ON \
  -D TPL_BLAS_LIBRARIES="-L$MKL_LIBDIR -lmkl_rt -ldl -lpthread -lm" \
  -D TPL_LAPACK_LIBRARIES="-L$MKL_LIBDIR -lmkl_rt -ldl -lpthread -lm" \
  $TRILINOS_PATH

# .... (disabled) code to enable metis / parmetis

  #PARMETIS_PATH=/usr/gapps/caprock/packages/parmetis-4.0.3/install

  #-D METIS_INCLUDE_DIRS:PATH=$PARMETIS_PATH/include \
  #-D METIS_LIBRARY_DIRS:PATH=$PARMETIS_PATH/lib \
  #-D METIS_LIBRARY_NAMES:STRING="metis" \
  #-D TPL_ENABLE_ParMETIS:BOOL=ON \
  #-D ParMETIS_INCLUDE_DIRS:PATH=$PARMETIS_PATH/include \
  #-D ParMETIS_LIBRARY_DIRS:PATH=$PARMETIS_PATH/lib \
  #-D ParMETIS_LIBRARY_NAMES:STRING="parmetis" \
  #-D TPL_ENABLE_METIS:BOOL=ON \

# .... make and make install

make -j $N
make install

