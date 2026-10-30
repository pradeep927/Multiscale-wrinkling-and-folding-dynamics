#!/bin/bash
#SBATCH --job-name=35x35
##Select the desired partition:
#SBATCH --partition=standard
##Define the màximum number of cores you need to use:
##SBATCH -c 8
##Define the màximun memory you will use:
#SBATCH --mem=20G
##Define the màximum time your job will run.
##With format Days-Hours (1day-0hours)
#SBATCH --time=1-00

echo "Start:" `date`
echo "I ran on:"
cd $SLURM_SUBMIT_DIR
echo $SLURM_NODELIST

#module load hiperlife/devel
#module load vtk/8.2.0

## Remove openMP
export OMP_NUM_THREADS=1

mpiexec hlVertexModel VMconfig.cfg

echo "End:" `date`
