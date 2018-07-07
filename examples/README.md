GENERATE PACKINGS with `run_jamming.sh`
==================================================
Scirpt `run_jamming.sh` generates: i) jammed packings 
for growing budding yeast cells; ii) a packing at the 
volume fraction $\delta\phi$; iii) the growth trajectory 
of the colony starting from two cells and ending  at $\delta\phi$.

Arguments to be set in `./run_jamming.sh`:    
1. att - range of the attractive forces
2. rate0 - cellular growth rate: a cell division per unit of time 
3. skip - frequency of saving a population structure
4. desync - cell growth desyncronisation parameter: gr=(1+[U(-1,1)-1/2] x dsync) x rate0
5. ar - aspect ration of the cell at birth. For example if ar=1.01, it means that a bud is 1% of the size of the mother cell
6. divtype - Location of a new born bud at the cell division (look at the jamming code for details)
7. P0 - feedback strenght: k(i)~exp(-P(i)/P0). P0=-1 for no feedback
8. Lx - Linear size of the periodic box
9. $\delta\psi = mm \cdot 10^{-nn}$
    * mm - mantissa
    * nn - exponent
10. se - seed for pseudo-random number generator    

The output files are saved in `$run_dir`.
`LF_*` files contain coordinates of the cells at jamming
and $\delta\psi$. `NC_*` file contains contact statistics
at the jamming (the first line) and at the (the second line).
`STATS_LF_*` file contain the statistics for every cell
for jammed and compressed packings. `v1.0_*` file  contains
the coordinations time-series of the growing population.   
To visualize a simulation output use `vissim` tool from 
[vistools](https://github.com/pgniewko/vistools) library.

For example, to visualize the output from
the `run_jamming.sh` script type:   

```bash
${VISTOOLS_PATH}/simbox.py v1.0_ar1.01_div_4_desync0.4_seed_1234_Lx8_Ly8_att0.0_P0.001.dat 8.0 0.0
```


RUN `run_shear.sh`
==================================================
Scirpt `run_shear.sh` reads a packig from the file 
and shears it. The scripts reads the structure generated
from `run_jamming.sh` script automatically, so the paramters, 
even if not used in the calculations, should match the one given
in `run_jamming.sh`. The code performs 1000 steps with a shear
strain equal to 1e-4. Finally, shear modulus G is estimated with
`estimate_G.py` script.    

Thre are only two arguments that need to be set:
1. att - range of the attractive forces
2. Lx - Linear size of the periodic box

The output files are saved in `$run_dir`.
`SHEAR_TRAJ_LF_DPHI_*` file contains the shearing simulations,
`G_data_LF_DPHI_*` contains shear stress and strain data for
shear modulus G calculations, and `GFIT_G_data_LF_DPHI_` file
contains numberical estimation the the shear modulus.   

To visualize a simulation output use `vissim` tool from 
[vistools](https://github.com/pgniewko/vistools) library.

To visualize the output from the `run_shear.sh` script type:   
```bash
${VISTOOLS_PATH}/simbox.py SHEAR_TRAJ_LF_DPHI_v1.0_ar1.01_div_4_desync0.4_seed_1234_Lx8_Ly8_att0.0_P0.001.dat 8.0 0.0001
```

COPYRIGHT NOTICE
================
Copyright (C) 2018, Pawel Gniewek     
PG: pawel.gniewek@berkeley.edu
All rights reserved.   
License: BSD-3  

REFERENCES
===============
1. "Jamming by Growth", P.G., C.S., O.H. (2018)
