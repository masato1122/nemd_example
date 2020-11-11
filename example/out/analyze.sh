PYTHONPATH=$PYTHONPATH:/home/ohnishi/work/nemd_2D/nemd_example
export PYTHONPATH

python ../../tools/plot_profile.py \
    --figname fig_profile.png \
    --lmpinput nemd0.in \
    --lmpdata data.lammps \
    --lmpdump nemd0.dump 


