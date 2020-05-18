PYTHONPATH=$PYTHONPATH:../

## for test
tnpt=100; tnemd=300
## for accurate simulation
#tnpt=500; tnemd=2000

python mk_graphite_nemd.py \
    -n 4 -m 3 \
    --nthermo 4 --ncenter 20 \
    --time_npt $tnpt --time_nemd $tnemd

cdir=./out
if [ ! -e $cdir ]; then
    mkdir $cdir
fi
mv *in *.xyz data.lammps $cdir
cp opt.tersoff $cdir


