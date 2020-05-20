nemd_dir=../../
PYTHONPATH=$PYTHONPATH:${nemd_dir}
tdir=${nemd_dir}/tools/

python $tdir/plot_profile.py
#exit

python $tdir/get_layer_temperature.py \
    --lmpdata  data.lammps \
    --lmpinput nemd2.in \
    --lmpdump  nemd2.dump \
    --outfile  temp_layer.txt

