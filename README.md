# nemd_example

Simple examples for non-equilibrium molecular dynamics (NEMD) simulation with LAMMPS package

# Requirement

python>=3.0, ase, pymatgen, pubchempy

# How to Use

To make LAMMPS input scripts:

``` discriptions
cd ~/example
PYTHONPATH=$PYTHONPATH:../nemd
python mk_graphite_nemd.py -n 4 -m 3
```

To make LAMMPS input scripts and run jobs:

``` simple way
cd .../example
sh mkneme.sh
cd ./{created_directory}
lmp_mpi < nemd0.in > log0.txt
lmp_mpi < nemd1.in > log1.txt
```

After the calculations, analyze results:

```
sh ../analyze_{direction_of_heat_flow}.sh
```

# Reference

* M. Ohnishi et al., Phys. Rev. B 95, 155405 (2017).

