#!/bin/bash
#PBS -N fftw_primer
#PBS -l nodes=1:ppn=24
#PBS -l walltime=00:10:00
#PBS -j oe             
  
module use $HOME/modulefiles 
module swap PrgEnv-cray PrgEnv-gnu
module swap gcc gcc/7.2.0
module load mpi/openmpi/3.0.0-gnu-7.2.0
# module load xpmem

# Change to the direcotry that the job was submitted from
# cd $PBS_O_WORKDIR

# Launch the parallel job to the allocated compute nodes
# mpirun -n 4 ../fftw/mpifftw  2048 4 6 1 0 > my_output_file 2>&1
mpirun -n 4  --map-by socket:span --bind-to core -report-bindings --display-map  -nooversubscribe ../fftw/mpifftw  $1 4 6 1 0 1
