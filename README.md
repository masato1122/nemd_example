# nemd_example

Simple example for NEMD simulation with LAMMPS package

# Installation

python>=3.0, ase, pymatgen

# How to Use

``` discriptions
cd ~/example
PYTHONPATH=$PYTHONPATH:../nemd
python mk_graphite_nemd.py -n 4 -m 3
```

``` simple way
cd ~/example
sh do.sh
cd ./out
lmp_mpi < nemd0.in > nemd0.log
lmp_mpi < nemd1.in > nemd1.log
```

after the calculation is finished,

```
sh analyze.sh
```

# Author

* Masato Ohnishi
* The University of Tokyo
* ohnishi@photon.t.u-tokyo.ac.jp

