#!/bin/bash
#SBATCH --job-name=35x35
##Select the partition:
#SBATCH --partition=mpi16c
##define the total number of threads to execute:
#SBATCH --ntasks=16
##define the number of nodes where the application will run
#SBATCH --nodes=1

#SBATCH --error=%x_%A.err
#SBATCH --output=%x_%A.out
#SBATCH --time=1-00
date

echo "Start:" `date`
echo "I ran on:"
#module load phdf5/1.10.1
#module load openmpi3

############################################
# Load ONLY compiler/MPI stack
############################################
module purge

module load GCC/11.2.0
module load OpenMPI/4.1.1-GCC-11.2.0
module load LAPACK/3.10.1-GCC-11.2.0
module load ScaLAPACK/2.1.0-gompi-2021b-fb

############################################
# Local HiPerLiFE installation
############################################
export HPLFE_ROOT=/clonetroop/kumar/local/new_hiperlife

export OMP_NUM_THREADS=1

############################################
# Run from the submission directory
############################################
cd $SLURM_SUBMIT_DIR

mpirun -np $SLURM_NTASKS /clonetroop/kumar/new_version_code/hiperlife_new/HPLFEApps/bin/hlvertex_model config.cfg



