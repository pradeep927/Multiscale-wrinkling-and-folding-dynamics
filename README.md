
# Multiscale wrinkling and folding dynamics

This work presents a multiscale computational framework for simulating inflating and deflating epithelial shells [[1](#References)] using a continuum active shell model and an acitve gel vertex model for epithelial sheets [[2](#References)].  

## License 
Distributed under the GNU GENERAL PUBLIC LICENSE. See [LICENSE](LICENSE) for details.


## Hiperlife [[3](#References)]

We use a parallel finite element library called hiperlife (High Performance Library for Finite Elements) (https://zenodo.org/doi/10.5281/zenodo.14927572), which serves as the core numerical engine for the simulations presented in this study. This library depends on other open-source libraries. The library is openly distributed to the community and is available online at {https://gitlab.com/hiperlife/hiperlife}. The aim of this library is to provide a computational framework to address problems of cell and tissue mechanobiology for a wide range of cases and users, with special focus on curved surfaces. The hiperlife is written in C++, uses the Message Passing Interface (MPI) paradigm for parallelism, and is built on top of several packages of the Trilinos Project. The installation of the hiperlife libraries can be carried out by following the guidelines provided at: {https://gitlab.com/hiperlife/hiperlife/-/blob/dev/INSTALL.md}. Prior to following the instructions for installation, it is necessary to clone the libraries via  Linux command {{git clone git@gitlab.com:hiperlife/hiperlife.git}} or download from {https://gitlab.com/hiperlife/hiperlife}.
Installation of hiperlife takes time (up to several hours) due to some dependencies 


## STEP 1: Installation of hiperlife

The installation of the hiperlife libraries can be carried out by following the guidelines provided at: {https://gitlab.com/hiperlife/hiperlife/-/blob/dev/INSTALL.md?ref_type=heads}. Prior to following the instructions for installation, it is necessary to clone the libraries via  Linux command {{git clone git@gitlab.com:hiperlife/hiperlife.git}} or download the hiperlife repository from {https://gitlab.com/hiperlife/hiperlife} by clicking the dropdown icon with the name "code".

## STEP 2: Code Organization and Project Setup

For the present work, we implemented the problem-specific codes within the hiperlife ecosystem by creating a top-level project folder named mechanics_of_epithelial_domes. This project folder contains a global CMakeLists.txt file that controls the compilation and build process, together with two application subdirectories corresponding to the two continuum models investigated in this study.

The folder shell_model implements the phenomenological continuum shell model used to describe epithelial shell inflation, viscoelastic relaxation, and deflation-induced buckling. The folder active_gel_bilayer_shell implements the continuum active gel bilayer shell model derived from active gel theories of the actomyosin cortex and described in Supplementary Note 1. This model explicitly accounts for the active, viscoelastic, and turnover dynamics of apical, basal, and lateral cortical surfaces through a homogenized shell formulation.
The folder vertex_model implements a curved-surface epithelial vertex model in which the tissue is represented as a network of polygonal cells connected through shared cell-cell interfaces. The model incorporates the mechanical contributions of cell area elasticity, perimeter contractility, and junctional tension and is used to investigate epithelial tissue mechanics, morphogenesis, and topological rearrangements. The formulation is designed to operate on curved epithelial geometries and provides a complementary discrete description of tissue mechanics that can be directly compared with the continuum shell-based approaches developed in this work.


Each application resides in its own folder and contains a local CMakeLists.txt file that defines the corresponding executable and model-specific source files. The code used in this study is openly available at https://github.com/pradeep927/mechanics_of_epithelial_domes.

The build system relies on CMake to configure and manage the compilation. Two additional files are located in the project root directory: cmake.project.ubuntu.20.04.sh, a shell script that automates the configuration process on Ubuntu systems, and userConfig.cmake, a configuration file where the paths to hiperlife, Trilinos, and other required dependencies are specified.

The typical installation procedure is as follows. First, install hiperlife following the instructions provided above. For the simulations reported here, the hiperlife libraries were installed on a workstation running Ubuntu 22.04. Next, place cmake.project.ubuntu.20.04.sh and userConfig.cmake in the top-level project folder mechanics_of_epithelial_domes. From a terminal, make the configuration script executable and launch the build process using the commands below:


chmod 777 cmake.project.ubuntu.20.04.sh
./cmake.project.ubuntu.20.04.sh


Enter the build directory and compile wih the following commands:

cd build
make -j4 install


This process generates the executable binary files at:

/home/ubuntu/shell_models_continuum/source_compiled/bin/hlactive_gel_bilayer_shell
/home/ubuntu/shell_models_continuum/source_compiled/bin/hlshell_model
/home/ubuntu/shell_models_continuum/source_compiled/bin/hlvertex_model

Note that /home/ubuntu/ is the home directory in our case. 

## STEP 3: Running simulations

To execute a simulation, we create a dedicated folder {run_simulation} that contains the mesh files and the configuration files. Meshes are provided in VTK format (or in .txt format for the vertex model simulation) and typically correspond to epithelial footprints with prescribed geometries. 

The simulation parameters, including references to the mesh files, are specified in a configuration file named {config.cfg}. This file allows the user to set model parameters, time-stepping controls, solver tolerances, and material constants.

A simulation is launched in parallel using MPI as:

mpirun -n 4 /home/ubuntu/shell_models_continuum/source_compiled/bin/hlactive_gel_bilayer_shell config.cfg



Here, the option {-n 4} specifies the number of processors. This value can be adjusted according to the available computational resources and the problem size.


The solution convergence at each timestep has been printed to the file slurm-93325.out. It can be opened with any text editor or viewed directly in the terminal using the Linux command cat slurm-93325.out.

If ran on a Dell XPS 13 equipped with an Intel Core i7 8th Gen processor (4 cores, 8 threads), the simulation of 10,000 timesteps for a mesh with 3,000 element nodes requires approximately 24 hours using 4 MPI processes.


The solution at each time is printed in VTK format with the filename sol_dis.{timestep}.vtk. As an example, for the active-gel tissue bilayer model, the inflation process occurs from timestep 1 to 2000, followed by a hold period from timestep 2000 to 8000. Deflation then takes place over 4000 steps. At timestep 2000, the volume reaches 100% and the dome is fully inflated, and the volume remains constant till the end of 8000th step. Then the deflation is applied uniformly in 4000 total steps with a constant deltat such that: by timestep 10000, the volume is reduced by 50%, and by timestep 11000, it is reduced by 75%. For illustration, some representative VTK files have been stored in the folder named results. During the simulations, a data file named globalIntegrals.dat is also generated, which contains information about the area, volume, and elastic energy at each timestep.

  

## STEP 4: Postprocessing and Visualization

Now we need to postprocess the generated output files in VTK format (e.g., {sol_dis_*.vtk}), which store the displacement and other field variables at different time steps. These files can be directly visualized and post-processed using ParaView. To streamline postprocessing, we provide ParaView state files ({*.pvsm}), stored in the folder paraview_state_file, which loads the output VTK files, applies predefined visualization settings, and enables rapid analysis of simulation results, including deformation fields, wrinkling patterns, and stress distributions.

## References

[1] Nimesh Chahare, Adam Ouzeri, Thomas Wilson, Pradeep K. Bal, Tom Golde, Guillermo Vilanova, Pau Pujol-Vives, Pere Roca-Cusachs, Xavier Trepat, Marino Arroyo. Multiscale wrinkling dynamics in epithelial shells. bioRxiv 2025.06.30.662426; doi: https://doi.org/10.1101/2025.06.30.662426

[2] Ouzeri A., Kale S., Chahare N., Torres-Sanchez A., Santos-Olivan D., Trepat X., Arroyo M. (2026), PRX Life, https://doi.org/10.1103/4mx5-t5hx

[3] Santos-Oliván, D., Vilanova, G., & Torres-Sánchez, A. (2025). hiperlife. Zenodo. https://doi.org/10.5281/zenodo.14927572


