#!/bin/bash

#############################################################
##  Set paths - Modify to fit to your installation
#############################################################
## Set main path to libraries installation
LIB_DIR=/home/ubuntu/local # <-------------------------------------- Modify

## Set paths
# TRILINOS_DIR=$LIB_DIR/trilinos
# VTK_DIR=$LIB_DIR/vtk
# GMSH_DIR=$LIB_DIR/gmsh
HIPERLIFE_DIR=$LIB_DIR/hiperlife/

## Set installation path
INSTALL_PATH=/home/ubuntu/mechanics_of_epithelial_domes/source_compiled # <--------------------- Modify



#############################################################
##  Configure
#############################################################

## Move to the build directory and delete
rm -rf $INSTALL_PATH/*
mkdir build
rm -rf build/* # <------------------------ Modify
cd build # <--------------------------------- Modify

## Configure
cmake \
    -D CMAKE_C_FLAGS="-march=native -O3"  \
    -D CMAKE_CXX_FLAGS="-march=native -O3"\
    -D CMAKE_INSTALL_PREFIX=$INSTALL_PATH \
    -D HL_BASE_PATH=$HIPERLIFE_DIR        \
    -D CMAKE_BUILD_TYPE=Release           \
    -D CMAKE_EXPORT_COMPILE_COMMANDS=YES  \
    ..

make -j 2 install
