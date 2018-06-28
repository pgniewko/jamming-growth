#! /bin/bash -x

MY_PATH=`pwd`
run_dir=$MY_PATH"/output"
exe_file=$MY_PATH"/../bin/jamming_by_growth"
cd $run_dir


att=0.0      # no attractions
cell="BUDD"  # type of the cell
celltype=2
rate0=0.002  # growth rate (in units of time-step)
skip=25      # save frame every $skip steps
AP=2.0       # NUMBER NOT USED
desync=0.4   # desynchronize growth rate
distrem=4.5  # NUMBER NOT USED
ar=1.01      # BUD OF THE SIZE 1% OF THE INITIAL SIZE
steps=10000  # NOT USED
divtype=4


P0=0.001
Lx=8
mm=1
nn=3
se=234


seed=$(($se+1000))
SEED2=$(( -1 * $seed)) # make seed negative

dphi=$mm"d-"$nn


Ly=$Lx

if [ ${P0%.*} -eq -1 ]
then 
suffix=${cell}_ar${ar}_div_${divtype}_desync${desync}_seed_$seed"_"Lx${Lx}_Ly${Ly}_a${AP}_att${att}.dat
else
suffix=${cell}_ar${ar}_div_${divtype}_desync${desync}_seed_$seed"_"Lx${Lx}_Ly${Ly}_a${AP}_att${att}_P${P0}.dat
fi

prodfile=prod_$suffix


time $exe_file <<EOF
  $ar
  $Lx
  $Ly
  $celltype
  $divtype
  $P0
  $att
  $rate0
  $desync
  $steps
  $SEED2
  $skip
  $dphi
  $prodfile
EOF


