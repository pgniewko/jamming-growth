#!/bin/bash -x

# Set directories
run_dir=/Users/pawel/School/2DGrowth/output/growth
exe_file=/Users/pawel/bin/2dgrowth/bin/growth_2d


# geometric parameters
ar=1.01
Lx=10.0
Ly=10.0
AP=2.0

# cell parameters
celltype=2
divtype=4
P0=-1
att=0.0

# run parameters
rate0=0.002
desync=0.4
steps=10000
seed=1001
#seed=$RANDOM
distrem=4.5
skip=20



cd $run_dir
# output files
if [ ${celltype} -eq 1 ]
then
cell="ellipse"
elif [ ${celltype} -eq 2 ]
then
cell="budding"
elif [ ${celltype} -eq 3 ]
then
cell="disk"
fi

if [ ${P0%.*} -eq -1 ]
then 
suffix=${cell}_ar${ar}_div${divtype}_desync${desync}_seed$(($seed))_Lx${Lx}_Ly${Ly}_a${AP}_att${att}.dat
else
suffix=${cell}_ar${ar}_div${divtype}_desync${desync}_seed$(($seed))_Lx${Lx}_Ly${Ly}_a${AP}_att${att}_P${P0}.dat
fi
prodfile=prod_$suffix
pressfile=press_$suffix
phifile=phi_$suffix
flowfile=flow_$suffix
treefile=clones_$suffix
agefile=age_$suffix
divfile=div_$suffix

# run program
time $exe_file <<EOF
  $ar
  $Lx
  $Ly
  $AP
  $celltype
  $divtype
  $P0
  $att
  $rate0
  $desync
  $steps
  $seed
  $distrem
  $skip
  $prodfile
  $pressfile
  $phifile
  $flowfile
  $treefile
  $agefile
  $divfile
EOF
