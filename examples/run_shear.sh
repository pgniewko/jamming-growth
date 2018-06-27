#! /bin/bash -x

MY_PATH=`pwd`
run_dir=$MY_PATH"/shear_dir"
exe_file=$MY_PATH"/../bin/shear_yeast_linearshear"
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

INPUT_LIST="RUN_LIST_Lx15_51_to_100.txt"

INPUT_LINE=`head -n $SLURM_ARRAY_TASK_ID $INPUT_LIST | tail -n 1`

P0=0.001
Lx=10
mm=1
nn=3
se=234


seed=$(($se+1000))
SEED2=$(( -1 * $seed)) # make seed negative

Ly=$Lx
dphi=$mm"e-"$nn


if [ ${P0%.*} -eq -1 ]
then
suffix=${cell}_ar${ar}_div_${divtype}_desync${desync}_seed_$seed"_"Lx${Lx}_Ly${Ly}_a${AP}_att${att}.dat
else
suffix=${cell}_ar${ar}_div_${divtype}_desync${desync}_seed_$seed"_"Lx${Lx}_Ly${Ly}_a${AP}_att${att}_P${P0}.dat
fi


file_="LF_DPHI_prod_"$suffix


#### SPRAWDZ CZY PLIK JEST DLUZSZY NIZ 1 LINIA
cp $MY_PATH"/output"/$file_  .


if [[ $(wc -l <$file_) -ge 2 ]] ; then

time $exe_file <<EOF
   $Lx
   $Ly
   $divtype
   $P0
   $att
   $rate0
   $desync
   $steps
   $SEED2
   $distrem
   $skip
   $dphi
   $file_
   1e-6
   1000
EOF

python $MY_PATH"/estimate_G.py" "G_data_"$file_ $dphi $P0

fi   

