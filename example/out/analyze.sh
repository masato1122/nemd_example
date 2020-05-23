nemd_dir=../../
PYTHONPATH=$PYTHONPATH:$nemd_dir
tdir=$nemd_dir/tools

for ical in 1 2; do
    
    python $tdir/plot_profile.py \
        --figname  fig_profile${ical}.png \
        --outfile temp_atom${ical}.txt \
        --lmpinput nemd1.in \
        --lmpdata  data.lammps \
        --lmpdump  nemd${ical}.dump

    python $tdir/get_layer_temperature.py \
        --lmpdata  data.lammps \
        --lmpinput nemd1.in \
        --lmpdump  nemd${ical}.dump \
        --outfile  temp_layer${ical}.txt
done

