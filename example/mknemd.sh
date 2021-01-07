#!/bin/sh
WORKDIR=$(cd $(dirname $0); pwd)

### Stractural Parameters
stage=0                          ## <0: FeCl3, =0: graphite, >0: FeCl3-GIC
#direction="out"; la=5; lb=5; lc=24;      # [nm]
#direction="in";  la=40; lb=5; lc=5
direction="in";  la=40; lb=2.4; lc=2.4

### temperatures [K]
Thot=305; Tcold=295

### damping parameters [ps]
tdamp=0.5
pdxy=10; pdz=1

## scaled length of thermostat
## length of thermostat should be longer than MFP*tdamp
if [ $direction == 'out' ]; then
    tnemd=2000
    lthermo=0.1
    lflow=$lc
else
    tnemd=1000
    if [ $stage -gt 0 ]; then
        mfp=20000   # m/s
    else
        mfp=8000   # m/s
    fi
    lflow=$la
    lpath=`echo "$mfp * $tdamp / 1000" | bc -l`  # nm
    lthermo=`echo "$lpath / $lflow" | bc -l`
    line_group="--lthermo ${lthermo}"
fi

## print length of middle regions
length_middle=`echo "$lflow * ( 0.5 - $lthermo )" | bc -l`
echo ""
echo " Length of middle region : ${length_middle} nm"
echo ""

## for both or out-of-plane and in-plane
tnpt=500
nloop=4

cdir=stage${stage}_${direction}_lz${lflow}nm
#cdir=test
fn=data.lammps
echo ""
echo " $cdir"
echo ""

python ../tools/mknemd_fecl3-gic.py \
    --datafile $fn \
    --stage $stage \
    --direction $direction \
    --la $la --lb $lb --lc $lc \
    --thot $Thot --tcold $Tcold \
    --time_npt $tnpt \
    --time_nemd $tnemd \
    --nloop $nloop \
    --tdamp $tdamp \
    --pdamp_xy $pdxy \
    --pdamp_z  $pdz \
    --thermostat langevin \
    --lthermo $lthermo

## save conditions
out2=PARAMETERS
cat >$out2<<EOF
stage     $stage
direction $direction

la        $la nm
lb        $lb nm
lc        $lc nm
Lmiddle   $length_middle nm

Thot      $Thot K
Tcold     $Tcold K

tnemd     $tnemd ps
tnpt      $tnpt ps
tdamp     $tdamp ps
pdxy      $pdxy ps
pdz       $pdz ps

lthermo   $lthermo (scaled length)
EOF

if [ ! -e $cdir ]; then
    mkdir $cdir
fi
mv *.in $fn POSCAR.nemd $out2 $cdir
cp opt.tersoff $cdir

