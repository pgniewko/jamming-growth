DESCRIPTION
==================================================

GENERATE PACKINGS with `run_jamming.sh`
==================================================
Scirpt `run_jamming.sh` generates jammed packings 
for growing budding yeast cells, a packing at the 
volume fraction $\delta\phi$, and the growth trajectory 
of the colony starting from two cells and ending 
with the packing at $\delta\psi$.

Arguments for `./run_jamming.sh`:    
1. att - 
2. rate0 - 
3. skip - 
4. desync - 
5. divtype - 
6. P0 - 
7. Lx - 
8. For ![dphi](https://latex.codecogs.com/gif.latex?%5Cdelta%5Cphi%3D%5C%24%5C%7Bmm%5C%7D%5Ccdot%2010%5E%7B-%5C%24%5C%7Bnn%5C%7D%7D)
    * mm -
    * nn -
9. se - 

The output files are saved in `$run_dir`.
`LF_*` files contain coordinates of the cells at jamming
and $\delta\psi$. `NC_*` file contains contact statistics
at the jamming (the first line) and at the (the second line).
`STATS_LF_*` file contain the statistics for every cell
for jammed and compressed packings. `v1.0_*` file  contains
the coordinations time-series of the growing population.:

To visualize a simulation output use `vissim` tool from 
[vistools](https://github.com/pgniewko/vistools) library.

For example, to visualize the output from 
the `run_jamming.sh` script type:

```bash
${VISTOOLS_PATH}/simbox.py v1.0_ar1.01_div_4_desync0.4_seed_1234_Lx8_Ly8_att0.0_P0.001.dat 8.0 0.0
```


RUN `run_shear.sh`
==================================================
To run ...

COPYRIGHT NOTICE
================
Copyright (C) 2018,  Pawel Gniewek     
PG: pawel.gniewek@berkeley.edu
All rights reserved.   
License: BSD-3  

REFERENCES
===============

