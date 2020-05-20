# nemd_example

Simple example for NEMD simulation with LAMMPS package

# Installation

python>=3.0, ase, pymatgen

# How to Use

To make LAMMPS input scripts:

``` discriptions
cd ~/example
PYTHONPATH=$PYTHONPATH:../nemd
python mk_graphite_nemd.py -n 4 -m 3
```

To make LAMMPS input scripts and run jobs:

``` simple way
cd ~/example
sh do.sh
cd ./out
lmp_mpi < nemd0.in > nemd0.log
lmp_mpi < nemd1.in > nemd1.log
```

After the calculations are finished, analyze results:

```
sh analyze.sh
```

# Reference

* M. Ohnishi et al., Phys. Rev. B 95, 155405 (2017).

