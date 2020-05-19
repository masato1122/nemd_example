#!/bin/sh
#PBS -q default
#PBS -l nodes=1:ppn=24
#PBS -j oe
#PBS -N test_nemd

export LANG=C
export OMP_NUM_THREADS=1
nprocs=24
cd $PBS_O_WORKDIR

MPIRUN=/opt/intel/compilers_and_libraries_2018.3.222/linux/mpi/intel64/bin/mpirun

$MPIRUN -np ${nprocs} lmp_mpi < nemd0.in > nemd0.log
$MPIRUN -np ${nprocs} lmp_mpi < nemd1.in > nemd1.log
$MPIRUN -np ${nprocs} lmp_mpi < nemd2.in > nemd2.log

