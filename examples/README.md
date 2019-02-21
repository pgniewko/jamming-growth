GENERATION OF JAMMED PACKINGS - `run_jamming.sh`
==================================================
Scirpt `run_jamming.sh` generates: i) a jammed packing 
of budding yeast cells; ii) a packing at the 
volume fraction $\delta\phi$; iii) a trajectory file that
contins the whole history of a growing budding yeast populations.

Arguments that can be set in `./run_jamming.sh`:    
1. att - range of the attractive forces; normally should be set to 0
2. rate0 - cellular growth rate: a cell division per unit of time 
3. skip - frequency of saving a population structure into a trajectory file
4. desync - cell growth desynchronisation parameter: gr=(1+[U(-1,1)-1/2] x dsync) x rate0
5. ar - aspect ratio of the cell at birth. For example if ar=1.01, it means that a new born bud is 1% of the size of the mother cell
6. divtype - Location of a new born bud at the cell division (look at the jamming source-code for details)
7. P0 - feedback strenght: k(i)=k0 x exp(-P(i)/P0). P0=-1 for no feedback
8. Lx - Linear size of the periodic box
9. $\delta\psi = mm \cdot 10^{-nn}$
    * mm - mantissa
    * nn - exponent
10. se - seed for pseudo-random number generator    

The output files are saved in `$run_dir`.
`LF_*` files contain coordinates of the cells at jamming
and at $\delta\psi$. `NC_*` file contains contact statistics
at the jamming (the first line) and at $\delta\psi$ (the second line).
`STATS_LF_*` file contains the statistics for each cell in the 
packing. 

`v1.0_*` file contains the trajectory of the growing population.   
This file can be rendered with `vissim` tool from 
[vistools](https://github.com/pgniewko/vistools) package.

For example, to visualize the output from
the `run_jamming.sh` script type:   

```bash
${VISTOOLS_PATH}/simbox.py v1.0_ar1.01_div_4_desync0.4_seed_1234_Lx8_Ly8_att0.0_P0.001.dat 8.0 0.0
```


SHEAR PACKINGS -  `run_shear.sh`
==================================================
Scirpt `run_shear.sh` reads cell coordinates from a file 
and shears the packing. The scripts requires the file generated
by `run_jamming.sh`, so the paramters (even though they are not used in the calculations), 
should match the one given in `run_jamming.sh`. 
The code performs 1000 steps with a shear strain equal to 1e-4. 
Finally, shear modulus G is estimated with `estimate_G.py` script.    

Thre are only two parameters that need to be set:
1. att - range of the attractive forces
2. Lx - Linear size of the periodic box

The output files are saved in `$run_dir`.
`SHEAR_TRAJ_LF_DPHI_*` file contains coordinates of sheard cells,
`G_data_LF_DPHI_*` contains shear stress and strain data used in
shear modulus G calculations, and `GFIT_G_data_LF_DPHI_` file
contains numerical estimation of the shear modulus.   

To visualize the simulation, use `vissim` tool from 
[vistools](https://github.com/pgniewko/vistools) package.

To visualize the output from the `run_shear.sh` script type:   
```bash
${VISTOOLS_PATH}/simbox.py SHEAR_TRAJ_LF_DPHI_v1.0_ar1.01_div_4_desync0.4_seed_1234_Lx8_Ly8_att0.0_P0.001.dat 8.0 0.0001
```

COPYRIGHT NOTICE
================
Copyright (C) 2018-2019, Pawel Gniewek     
PG: pawel.gniewek@berkeley.edu
All rights reserved.   
License: BSD-3  

REFERENCES
===============
1. [Jamming by growth](https://arxiv.org/abs/1810.01999), P.G., C.S., O.H. (2019)
