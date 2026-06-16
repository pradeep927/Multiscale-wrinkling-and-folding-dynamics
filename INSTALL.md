# Installation instructions
Installation guide for Activel gel vertex model of spherical domes
This document explains how to install and set up the vertex model application on your system.


## Downloading and installing hiperlife

Prior to following the instructions for installing the vertex model application, it is necessary to download the hiperlife library [[1](#References)] available [here](https://doi.org/10.5281/zenodo.14927572) (or [here](https://gitlab.com/hiperlife/hiperlife/-/tree/release-5.0.0-beta?ref_type=tags) for latest versions), and follow the INSTALL.md instructions therein. Note that installing the hiperlife library and it's dependencies can take over an hour.
For completeness, we summarize the necessary steps here.

#### Dependencies

By default, all dependencies are downloaded, compiled and installed directly by the hiperlife's build system.
Below is an explicit list of the dependencies.

General dependencies
```bash
sudo apt-get install git gcc g++ gfortran cmake libopenmpi-dev libboost-filesystem-dev liblapack-dev libscalapack-openmpi-dev
```

Required by Scotch
```bash
sudo apt-get install flex bison yacc
```

Required by VTK
```bash
sudo apt-get install qtbase5-dev libgl1-mesa-dev libopengl-dev libglx-dev
```
#### Installation

Once the library has been downloaded, create a build folder in the top directory
```bash
mkdir build
cd build
```
then execute cmake to generate the build system files (MakeFile, etc.)

```bash
cmake .. -DCMAKE_INSTALL_PREFIX=<path/to/install/hiperlife>
```
where <path/to/install/hiperlife> is where the library will be installed, e.g., in a new folder ~/.local/hiperlife.

Then compile with :

```bash
make -jX dependencies
```
where *X* is the number of cores (e.g. *8*).
Once all the dependencies are ready, hiperlife can be compiled and installed :
```bash
cmake ..
make -jX install
```


## Installing agvm_domes

We followed the instruction given [here](https://gitlab.com/hiperlife/hl-base-project) for creating a hiperlife-base project. 

This project uses **CMake** with multiple `CMakeLists.txt` files:
- `CMakeLists.txt` (root) → sets up the project and adds subdirectories.
- `src/CMakeLists.txt` → builds the core library (`aux`) and the main program.

Once the code has been downloaded, create a build folder in the top directory :
```bash
mkdir build
cd build
```
Then configure the project :

```bash
cmake -DHL_BASE_PATH=<path/to/install/hiperlife>  -DCMAKE_INSTALL_PREFIX=<path/to/install/app> -DCMAKE_BUILD_TYPE=RELEASE ..
```
where <path/to/install/hiperlife> is the path where the library has been installed (e.g. ~/.local/hiperlife) and <path/to/install/app> is the path where the application executable will be installed (e.g. in a new folder ~/.local/agvm_domes). Note that the name of the executable is *hlVXDomes*.

## Installation remarks

The installation was done with the following setup :
- Dell Intel® Xeon(R) Silver 4208 CPU @ 2.10GHz × 32
- 62.5 GiB memory 
- Ubuntu 20.04.6 LTS
- hiperlife-release-5.0.0-beta (also available [here](https://gitlab.com/hiperlife/hiperlife/-/tree/release-5.0.0-beta?ref_type=tags))

When installing hiperlife and the dependencies, here a few ad-hoc modifications that might be necessary :
1) Remember to tell your OS the location of the installed libraries. For UNIX systems you can add the paths of these libraries to the environment variable LD_LIBRARY_PATH. To avoid adding the paths to each terminal session, you can put these instructions in the ~/.bashrc file. As an example, the following code is an extract of the ~/.bashrc file in which we are appending the paths to the openmpi and hiperlife libraries and then exporting the environment variable LD_LIBRARY_PATH.


```bash
# Paths to libraries
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/openmpi/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/openmpi/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$HOME/.local/hiperlifelib/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$HOME/.local/hiperlife:$LD_LIBRARY_PATH
```

Same for the path to the agvm_domes executable
```bash
export PATH="$HOME/.local/agvm_domes/:$PATH"

```


Then the ~/.bashrc file needs to be sourced again
```bash
source ~/.bashrc
```


2) "sudo apt-get install flex bison" instead of "sudo apt-get install flex bison yacc"
3) Make sure to have an updated g++ compiler. Check your compiler version with 
```bash
g++ --version
```
4) If <path/to/install/hiperlife> was set to ~/.local/hiperlife, in ~/.local/hiperlife/lib/cmake/hiperlife/hiperlifeConfig.cmake:42 use
```bash
    set_and_check(hiperlife_MPI_DIR /usr/lib/x86_64-linux-gnu/openmpi)
```
instead of 
```bash
set_and_check(hiperlife_MPI_DIR /usr/lib/x86_64-linux-gnu/openmpi/openmpi)
```

## References

[1] Santos-Oliván, D., Vilanova, G., & Torres-Sánchez, A. (2025). hiperlife. Zenodo. https://doi.org/10.5281/zenodo.14927572