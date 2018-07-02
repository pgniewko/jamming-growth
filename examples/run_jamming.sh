#! /bin/bash -x

MY_PATH=`pwd`
run_dir=$MY_PATH"/jamming_dir"
exe_file=$MY_PATH"/../bin/jamming_by_growth"
cd $run_dir


att=0.0      # no attractions
rate0=0.002  # growth rate (in units of time-step)
skip=25      # save frame every $skip steps
desync=0.4   # desynchronize growth rate
ar=1.01      # initial bud size is 1% of the mother cell
divtype=4    # random location of a new born bud


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
suffix=ar${ar}_div_${divtype}_desync${desync}_seed_${seed}_Lx${Lx}_Ly${Ly}_att${att}.dat
else
suffix=ar${ar}_div_${divtype}_desync${desync}_seed_${seed}_Lx${Lx}_Ly${Ly}_att${att}_P${P0}.dat
fi

version=1.0
prodfile=v${version}_${suffix}


time $exe_file <<EOF
  $ar
  $Lx
  $Ly
  $divtype
  $P0
  $att
  $rate0
  $desync
  $SEED2
  $skip
  $dphi
  $prodfile
EOF


