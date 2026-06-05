
<<<<<<< HEAD
We use a parallel finite element library called hiperlife (High Performance Library for Finite Elements) (https://zenodo.org/doi/10.5281/zenodo.14927572), which serves as the core numerical engine for the simulations presented in this study. This library depends on other open-source libraries. The library is openly distributed to the community and is available online at {https://gitlab.com/hiperlife/hiperlife}. The aim of this library is to provide a computational framework to address problems of cell and tissue mechanobiology for a wide range of cases and users, with special focus on curved surfaces. The hiperlife is written in C++, uses the Message Passing Interface (MPI) paradigm for parallelism, and is built on top of several packages of the Trilinos Project. 
=======
We use a parallel finite element library called hiperlife (High Performance Library for Finite Elements) (https://zenodo.org/doi/10.5281/zenodo.14927572), which serves as the core numerical engine for the simulations presented in this study. This library depends on other open-source libraries. The library is openly distributed to the community and is available online at {https://gitlab.com/hiperlife/hiperlife}. The aim of this library is to provide a computational framework to address problems of cell and tissue mechanobiology for a wide range of cases and users, with special focus on curved surfaces. The hiperlife is written in C++, uses the Message Passing Interface (MPI) paradigm for parallelism, and is built on top of several packages of the Trilinos Project. The installation of the hiperlife libraries can be carried out by following the guidelines provided at: {https://gitlab.com/hiperlife/hiperlife/-/blob/dev/INSTALL.md}. Prior to following the instructions for installation, it is necessary to clone the libraries via  Linux command {{git clone git@gitlab.com:hiperlife/hiperlife.git}} or download from {https://gitlab.com/hiperlife/hiperlife}.
Installation takes approximately 1 hour on a Dell XPS 13 equipped with an Intel Core i7 8th Gen processor running Ubuntu 22.04.
>>>>>>> 94a7c0c9ce9434b6d58336cb1afb9abb97deeb29

STEP 1: {Installation of hiperlife}

The installation of the hiperlife libraries can be carried out by following the guidelines provided at: {https://gitlab.com/hiperlife/hiperlife/-/blob/dev/INSTALL.md?ref_type=heads}. Prior to following the instructions for installation, it is necessary to clone the libraries via  Linux command {{git clone git@gitlab.com:hiperlife/hiperlife.git}} or download the hiperlife repository from {https://gitlab.com/hiperlife/hiperlife} by clicking the dropdown icon with the name "code".

STEP 2: {Code Organization and Project Setup}

For the present work, we implemented the problem-specific code within the hiperlife ecosystem by creating a top-level project folder named {wrinkling_epithelial_continuum}. This project folder contains a global {CMakeLists.txt} file, which controls the compilation, and one or more subprojects (applications), each corresponding to a specific simulation setup. For example, the application {inflate_hold_deflate} implements the model of interest here and resides in its own folder with a local \texttt{CMakeLists.txt}.
 The code used in this study is openly available at {https://github.com/pradeep927/mechanics_of_epithelial_domes}.

The build system relies on CMake to configure and manage the compilation. Two additional files are kept in the project root: {cmake.project.ubuntu.20.04.sh}, a shell script to automate the configuration process on Ubuntu systems; and \texttt{userConfig.cmake}, a configuration file where paths to hiperlife, Trilinos, and other dependencies are specified. The typical installation procedure is as follows: Install hiperlife following the instructions mentioned above. For our simulations, the hiperlife libraries were installed on a workstation running Ubuntu 22.04. Place these two files inside the top-level project folder {wrinkling_epithelial_continuum}. From the terminal, make the script executable and launch the build configuration through the following commands:

chmod 777 cmake.project.ubuntu.20.04.sh
./cmake.project.ubuntu.20.04.sh


Enter the build directory and compile wih the following commands:

cd build
make -j4 install


This process generates the executable binary at:

/home/ubuntu/wrinkling_epithelial_continuum/source_compiled/bin/hlinflate_hold_deflate

Note that /home/ubuntu/ is the home directory in our case. 

STEP 3: {Running Simulations}

To execute a simulation, we create a dedicated folder {run_simulation} that contains the mesh files and the configuration files. Meshes are provided in VTK format and typically correspond to epithelial footprints with prescribed geometries. For example, a mesh named {aspect_26} corresponds to a circular domain with aspect ratio 26. Boundary conditions are imposed by carefully fixing the exterior nodes of the domain using elastic springs, which mimic physical anchoring constraints.

The simulation parameters, including references to the mesh files, are specified in a configuration file named {config.cfg}. This file allows the user to set model parameters, time-stepping controls, solver tolerances, and material constants.

A simulation is launched in parallel using MPI as:

mpirun -n 4 /home/ubuntu/wrinkling_epithelial_continuum/source_compiled/bin/hlinflate_hold_deflate config.cfg



Here, the option {-n 4} specifies the number of processors. This value can be adjusted according to the available computational resources and the problem size.
<<<<<<< HEAD

The solution convergence at each timestep has been printed to the file slurm-55239.out. It can be opened with any text editor or viewed directly in the terminal using the Linux command cat slurm-55239.out.
=======
If run on a Dell XPS 13 equipped with an Intel Core i7 8th Gen processor (4 cores, 8 threads), the simulation of 5,000 timesteps for a mesh with 3,000 element nodes requires approximately 24 hours using 4 MPI processes.
{Postprocessing and Visualization}
>>>>>>> 94a7c0c9ce9434b6d58336cb1afb9abb97deeb29

The solution at each time is printed in VTK format with the filename sol_dis.{timestep}.vtk. The inflation process occurs from timestep 1 to 500, followed by a hold period from timestep 500 to 5000. Deflation then takes place over 4000 steps. At timestep 500, the volume reaches 100% and the dome is fully inflated, and the volume remains constant till the end of 5000th step. Then the deflation is applied uniformly in 4000 total steps with a constant deltat such that: by timestep 7000, the volume is reduced by 50%, and by timestep 8000, it is reduced by 75%. For illustration, some representative VTK files have been stored in the folder named Expected results. During the simulations, a data file named globalIntegrals.dat is also generated, which contains information about the area, volume, and elastic energy at each timestep.  

STEP 4: {Postprocessing and Visualization}

Now we need to postprocess the generated output files in VTK format (e.g., {sol_dis_*.vtk}), which store the displacement and other field variables at different time steps. These files can be directly visualized and post-processed using ParaView. To streamline postprocessing, we provide a ParaView state file ({paraview_state.pvsm}), which loads the output files, applies predefined visualization settings, and enables rapid analysis of simulation results, including deformation fields, wrinkling patterns, and stress distributions.
