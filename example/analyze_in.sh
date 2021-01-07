#!/bin/sh
WORKDIR=$(cd $(dirname $0); pwd)
dir=$WORKDIR/../tools

python $dir/plot_pressure.py \
    -f log0.txt \
    -i 2

python $dir/plot_profile_bin.py \
    -d npt.dump \
    -i nemd0.in \
    --log log0.txt \
    --lbin 3.5 \
    --figname fig_profile0.png \
    --outdir out_npt

for i in `seq 1 10`; do
    fn=nemd${i}.dump
    if [ ! -e $fn ]; then
        continue
    fi
    python $dir/plot_profile_bin.py \
        -d nemd${i}.dump \
        -i nemd${i}.in \
        --log log${i}.txt \
        --lbin 3.5 \
        --outdir out_ave${i} \
        --figname fig_profile${i}.png

done

