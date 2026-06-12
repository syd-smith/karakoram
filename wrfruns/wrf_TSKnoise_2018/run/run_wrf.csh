#!/bin/csh

#SBATCH --account=strong-kp
#SBATCH --partition=strong-kp 
#SBATCH --time=256:00:00         
#SBATCH --nodes=1             # number of cluster nodes, abbreviated by -N
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#SBATCH -J wrf_2018_noise
#SBATCH --ntasks=16           # number of MPI tasks, abbreviated by -A

#load modules
module purge
module use /uufs/chpc.utah.edu/common/home/u1301408/MyModules
module load wrf/2026.4.lua

#run program

unlimit stacksize
#ulimit -s unlimited
#limit stacksize unlimited
mpirun -np $SLURM_NTASKS ./wrf.exe
#unlimit stacksize

#mpirun -x OMP_NUM_THREADS=1 -np $SLURM_NTASKS ./real.exe
#mpirun -x OMP_NUM_THREADS=1 -np $SLURM_NTASKS ./wrf.exe

