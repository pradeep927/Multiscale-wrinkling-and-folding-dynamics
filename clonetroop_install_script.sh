module load CMake
module load hiperlife

## Set installation path
INSTALL_PATH=$HOME/local/NimeshPaperExecutables # <------------------ Modify


#############################################################
##  Configure
#############################################################

## Move to the build directory and delete
mkdir build
echo "You are about to remove all the contents of
 $(cd ./build/; pwd) \n Please, check that the path is correct."
rm -rf -I ./build/*
cd ./build

## Configure
cmake \
    -D CMAKE_INSTALL_PREFIX=$INSTALL_PATH       \
    -D CMAKE_CXX_COMPILER:STRING="mpicxx"       \
    -D CMAKE_C_COMPILER:STRING="mpicc"          \
    -D CMAKE_C_FLAGS:STRING="-O3 -march=native" \
    -D CMAKE_CXX_FLAGS:STRING="-O3 -march=native -std=c++14"   \
    -D MPI_BASE_PATH=${EBROOTOPENMPI}           \
    -D TRI_BASE_PATH=${EBROOTTRILINOS}          \
    -D VTK_BASE_PATH=${EBROOTVTK}               \
    -D GMSH_BASE_PATH=$GMSH_DIR                 \
    -D HPLFE_BASE_PATH=${EBROOTHIPERLIFE}       \
    -D CMAKE_BUILD_TYPE=RELEASE                 \
    -D CMAKE_EXPORT_COMPILE_COMMANDS=YES        \
    ..
