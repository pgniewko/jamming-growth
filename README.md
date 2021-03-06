[![DOI](https://zenodo.org/badge/138428274.svg)](https://zenodo.org/badge/latestdoi/138428274)

DESCRIPTION
==================================================
2D numerical simulations of budding yeast populations growing in confined environment.
This code accompanies the manuscript: "Biomechanical Feedback Strengthens Jammed Cellular Packings", by P. Gniewek, C. Schreck, and O. Hallatschek, PRL (2019).

GETTING THE CODE
==================================================
* To get the code:  
```
git clone git@github.com:pgniewko/jamming-growth.git
```

* To obtain the most recent version of the code:   
```
git pull origin master
```

COMPILING 
==================================================
* This repo contains two codes:   
    + `jamming_by_growth.f` - Generate jammed packings of growing budding yeast population    
    + `shear_yeast_linearshear.f` - Shear a packing of budding yeast cells   

* Building executable files:
    + Building only a code for jamming of budding yeast:    
        `make `       
    + Building all of the executables:    
        `make all`

* Compiling parameters can be changed in `config.mk` file:
    + FORTRAN compiler:
        `FC := gfortran`

* If the compilation is successful, binary files are located in `bin` directory.       

RUN 
==================================================
To learn how to run the compiled codes check out scripts in `examples` directory.


COPYRIGHT NOTICE
================
Copyright (C) 2018-2019,  Pawel Gniewek & Carl Schreck    
PG: pawel.gniewek@berkeley.edu    
CS: carl.schreck@berkeley.edu     
All rights reserved.    
License: BSD-3  

REFERENCES
===============
1. [Biomechanical feedback strengthens jammed cellular packings](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.122.208102), P.G., C.S., O.H. Phys. Rev. Lett. (2019)
2. [Self-driven jamming in growing microbial populations](https://www.nature.com/articles/nphys3741), M.D., J.H., C.S., P.G., L.H., S.H., O.H., Nature Phys (2016)
