#! /bin/bash -x

MY_PATH=`pwd`
run_dir=$MY_PATH"/shear_dir"
exe_file=$MY_PATH"/../bin/shear_yeast_linearshear"
mkdir -p ${run_dir}
cd $run_dir

att=0.0      # no attractions
desync=0.4   # desynchronize growth rate
ar=1.01      # initial bud size is 1% of the mother cell
divtype=4    # random location of a new born bud

P0=0.001
Lx=10
mm=1
nn=3
se=234


seed=$(($se+1000))
SEED2=$(( -1 * $seed)) # make seed negative

Ly=$Lx
dphi=${mm}e-${nn}


if [ ${P0%.*} -eq -1 ]
then
suffix=ar${ar}_div_${divtype}_desync${desync}_seed_${seed}_Lx${Lx}_Ly${Ly}_att${att}.dat
else
suffix=ar${ar}_div_${divtype}_desync${desync}_seed_${seed}_Lx${Lx}_Ly${Ly}_att${att}_P${P0}.dat
fi

version=1.0
file_=LF_DPHI_v${version}_${suffix}


cp ${MY_PATH}/jamming_dir/${file_}  .

## PROCEED ONLY IF THE FILE IS LONGER THAN TWO LINES
if [[ $(wc -l <$file_) -ge 2 ]] ; then

time $exe_file <<EOF
   $Lx
   $Ly
   $att
   $file_
   1e-6
   1000
EOF

python $MY_PATH"/estimate_G.py" "G_data_"$file_ $dphi $P0

fi   

