

We use a parallel finite element library called hiperlife (High Performance Library for Finite Elements) (https://zenodo.org/doi/10.5281/zenodo.14927572), which serves as the core numerical engine for the simulations presented in this study. This library depends on other open-source libraries. The library is openly distributed to the community and is available online at {https://gitlab.com/hiperlife/hiperlife}. The aim of this library is to provide a computational framework to address problems of cell and tissue mechanobiology for a wide range of cases and users, with special focus on curved surfaces. The hiperlife is written in C++, uses the Message Passing Interface (MPI) paradigm for parallelism, and is built on top of several packages of the Trilinos Project. 

We use a parallel finite element library called hiperlife (High Performance Library for Finite Elements) (https://zenodo.org/doi/10.5281/zenodo.14927572), which serves as the core numerical engine for the simulations presented in this study. This library depends on other open-source libraries. The library is openly distributed to the community and is available online at {https://gitlab.com/hiperlife/hiperlife}. The aim of this library is to provide a computational framework to address problems of cell and tissue mechanobiology for a wide range of cases and users, with special focus on curved surfaces. The hiperlife is written in C++, uses the Message Passing Interface (MPI) paradigm for parallelism, and is built on top of several packages of the Trilinos Project. The installation of the hiperlife libraries can be carried out by following the guidelines provided at: {https://gitlab.com/hiperlife/hiperlife/-/blob/dev/INSTALL.md}. Prior to following the instructions for installation, it is necessary to clone the libraries via  Linux command {{git clone git@gitlab.com:hiperlife/hiperlife.git}} or download from {https://gitlab.com/hiperlife/hiperlife}.
Installation takes approximately 1 hour on a Dell XPS 13 equipped with an Intel Core i7 8th Gen processor running Ubuntu 22.04.


STEP 1: {Installation of hiperlife}

The installation of the hiperlife libraries can be carried out by following the guidelines provided at: {https://gitlab.com/hiperlife/hiperlife/-/blob/dev/INSTALL.md?ref_type=heads}. Prior to following the instructions for installation, it is necessary to clone the libraries via  Linux command {{git clone git@gitlab.com:hiperlife/hiperlife.git}} or download the hiperlife repository from {https://gitlab.com/hiperlife/hiperlife} by clicking the dropdown icon with the name "code".

STEP 2: {Code Organization and Project Setup}

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

STEP 3: {Running Simulations}

Each model under `run_simulation` uses four separate folders:

| Folder | Purpose |
| --- | --- |
| `input_short/` | Configuration and all mesh/connectivity files for a quick test |
| `input_full/` | Configuration and all mesh/connectivity files for the full simulation |
| `expected_output_short/` | Reference results from a tested short simulation |
| `expected_output_full/` | Reference results from the full simulation |

### Run manually with full executable paths

You can use your usual `mpirun` command. First copy the selected inputs into a
fresh working folder so the application cannot overwrite the reference files.
Load your HiPerLife/MPI library environment before running these commands.
The example below assumes the project is installed at
`/home/ubuntu/Multiscale-wrinkling-and-folding-dynamics-main`; replace this path
with your actual project path if different.

For a short active-gel bilayer simulation:

```bash
# Go to the project and create the parent directory for generated runs.
cd /home/ubuntu/Multiscale-wrinkling-and-folding-dynamics-main
mkdir -p runs

# Create a fresh folder, copy the configuration and mesh, and launch there.
# The commands stop if this run folder already exists.
mkdir runs/active_gel_short_001 &&
cp -R run_simulation/active_gel_bilayer_shell/input_short/. runs/active_gel_short_001/ &&
cd runs/active_gel_short_001 &&
mpirun -n 4 \
  /home/ubuntu/Multiscale-wrinkling-and-folding-dynamics-main/source_compiled/bin/hlactive_gel_bilayer_shell \
  config_shell.cfg
```

New results are written to `runs/active_gel_short_001/`. For another run, use a
new folder name, such as `active_gel_short_002`, in all three commands. For a full
simulation, copy from `input_full/` and use a new folder such as
`active_gel_full_001`. The continuum configuration filename is `config_shell.cfg`,
not `config.cfg`.

For a short shell-model simulation:

```bash
cd /home/ubuntu/Multiscale-wrinkling-and-folding-dynamics-main
mkdir -p runs
mkdir runs/shell_short_001 &&
cp -R run_simulation/shell_model/input_short/. runs/shell_short_001/ &&
cd runs/shell_short_001 &&
mpirun -n 4 \
  /home/ubuntu/Multiscale-wrinkling-and-folding-dynamics-main/source_compiled/bin/hlshell_model \
  config_shell.cfg
```

The vertex model follows the same procedure once its input folders are prepared:
copy the selected vertex inputs into a fresh run folder, enter that folder, and
launch `source_compiled/bin/hlvertex_model` using its full path and `config.cfg`.
Never launch a simulation inside an `input_*` or `expected_output_*` folder.

### Optional automatic launcher

The common launcher automates the same copy-and-run procedure and creates a
unique run folder automatically. Use it from the repository root:

```bash
./run_simulation.sh active_gel_bilayer_shell short 4
./run_simulation.sh shell_model short 4
./run_simulation.sh vertex_model short 4
```

Replace `short` with `full` for the full simulation. The final number is the number
of MPI processes; it defaults to 4 if omitted. Load your HiPerLife/MPI library
environment first. The launcher uses executables installed in
`source_compiled/bin`. For another installation, set the directory explicitly:

```bash
export HPLFEAPPS_BIN_DIR="$HOME/local/HPLFEApps/bin"
./run_simulation.sh shell_model short 4 --dry-run
./run_simulation.sh shell_model short 4
```

The required executable names are `hlactive_gel_bilayer_shell`, `hlshell_model`,
and `hlvertex_model`. Set `MPIEXEC` to another MPI launcher executable if needed;
the launcher must accept `-n NUMBER EXECUTABLE CONFIGURATION`.

Each launch copies the selected inputs into a unique folder under `runs/` and
runs the application there. New outputs and `run.log` are written in that folder;
reference inputs and expected outputs are left intact. `launch-command.txt`
records the command. `runs/` is excluded from version control. `--dry-run` checks
inputs, executable paths, and launcher availability without starting a simulation
or creating a run folder. It does not test binary compatibility or numerical results.

Continuum inputs use `config_shell.cfg` and the VTK mesh named by `prefixMesh`.
Vertex inputs use `config.cfg` and all required `vertexmesh_*.txt`, mesh, and
connectivity files, including `neighbourcells.txt` and `faceIDsinCell.txt`.
Prepare both vertex input folders before launching that model: the launcher does
not copy inputs from reference output folders or invent short-run parameters.

The current continuum short configurations request 50 active-gel steps and 35
shell steps. The launcher preserves the selected configuration exactly; it does
not alter physical parameters, time controls, or active-gel `v_target = 0.7`.

The solution convergence at each timestep has been printed to the file slurm-94546.out for the active gel tissue bilayer shell model, and slurm-90208.out for the shell model. It can be opened with any text editor or viewed directly in the terminal using the Linux command cat slurm-94546.out.

If run on a Dell XPS 13 equipped with an Intel Core i7 8th Gen processor (4 cores, 8 threads), the simulation of 10,000 timesteps for a mesh with 3,000 element nodes requires approximately 24 hours using 4 MPI processes.


The solution at each time is printed in VTK format with the filename sol_dis.{timestep}.vtk. As an example, for the active-gel tissue bilayer model, the inflation process occurs from timestep 1 to 2000, followed by a hold period from timestep 2000 to 6000. Deflation then takes place over 4000 steps. From step 9700 onward, the tissue relaxes while remaining adhered to the substrate toward the end of deflation. At timestep 2000, the volume reaches 100% and the dome is fully inflated, and the volume remains constant till the end of 6000th step. Then the deflation is applied uniformly in 4000 total steps with a constant deltat such that: by timestep 8000, the volume is reduced by 50%, and by timestep 9000, it is reduced by 75%. For illustration, some representative VTK files have been stored in the folder named results. During the simulations, a data file named globalIntegrals.dat is also generated, which contains information about the area, volume, and elastic energy at each timestep. For the shell model, the inflation process occurs from timestep 1 to 3000, followed by a hold period from timestep 3000 to 8000. Deflation then takes place over 4000 steps. 
  

STEP 4: {Postprocessing and Visualization}

Now we need to postprocess the generated output files in VTK format (e.g., {sol_dis_*.vtk}), which store the displacement and other field variables at different time steps. These files can be directly visualized and post-processed using ParaView. To streamline postprocessing, we provide ParaView state files ({*.pvsm}), stored in the folder `run_simulation/post_processing`, which loads the output VTK files, applies predefined visualization settings, and enables rapid analysis of simulation results, including deformation fields, wrinkling patterns, and stress distributions.
