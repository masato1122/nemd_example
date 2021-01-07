#!/bin/sh
WORKDIR=$(cd $(dirname $0); pwd)
tdir=$WORKDIR/../tools

python $tdir/plot_pressure.py \
    -f log0.txt \
    -i 2

python $tdir/plot_profile.py \
    --figname  fig_profile_npt.png \
    --outfile temp_atom_npt.txt \
    --lmpinput nemd0.in \
    --lmpdata  data.lammps \
    --log      log0.txt \
    --lmpdump  npt.dump

for ical in 0 1 2 3 4 5; do

    fdump=nemd${ical}.dump
    if [ ! -e $fdump ]; then
        continue
    fi

    python $tdir/plot_profile.py \
        --figname  fig_profile${ical}.png \
        --outfile temp_atom${ical}.txt \
        --lmpinput nemd${ical}.in \
        --lmpdata  data.lammps \
        --log      log${ical}.txt \
        --lmpdump  $fdump
    
    python $tdir/get_layer_temperature.py \
        --lmpdata  data.lammps \
        --lmpinput nemd${ical}.in \
        --lmpdump $fdump \
        --outfile  temp_layer${ical}.txt
    
    i1=${ical}
done

if [ ${i1} -lt 2 ]; then
    exit
fi
for mm in atom layer; do
    outfile=average_${mm}.txt
    python $tdir/get_average.py \
        --prefix temp_${mm} \
        --irange 2:${i1} \
        --outfile $outfile \
        --figname fig_ave_${mm}.png
done

