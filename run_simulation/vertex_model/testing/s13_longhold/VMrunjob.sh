#!/bin/bash

#SBATCH --job-name=20x20
#SBATCH --error=%x_%A.err
#SBATCH --output=%x_%A.out
#SBATCH --time=1-00
echo "Start:" `date`
echo "I ran on:"
echo $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR
echo $SLURM_NODELIST

#module load hiperlife/devel
#module load vtk/8.2.0


## Remove openMP
export OMP_NUM_THREADS=1

mpirun -n $SLURM_NTASKS hlvertex_model config.cfg 
