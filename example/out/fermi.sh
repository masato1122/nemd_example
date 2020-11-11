#!/bin/sh
#PBS -q default
#PBS -l nodes=1:ppn=40
#PBS -j oe
#PBS -N test

export LANG=C
export OMP_NUM_THREADS=1
nprocs=40
cd $PBS_O_WORKDIR
rm test.o*

MPIRUN=/opt/intel/compilers_and_libraries_2018.3.222/linux/mpi/intel64/bin/mpirun
$MPIRUN -np ${nprocs} lmp_mpi < nemd0.in > log.txt

