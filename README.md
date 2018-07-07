DESCRIPTION
==================================================


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
* The repo contains two codes:   
    + `jamming_by_growth.f` - to generated jammed packings of growing population    
    + `shear_yeast_linearshear.f` - to shear the inanimate packing of dumb-bell cells   

* Building executable files:
    + To build just a code for budding yeast jamming type:    
        `make `       

    + To build all of the executables type in the command line:    
        `make all`

* If the compilation terminates with a success, the binary files are located in `bin` directory.       

RUN 
==================================================
To learn how to run the compiled codes check out scripts in `examples` directory.


COPYRIGHT NOTICE
================
Copyright (C) 2018,  Pawel Gniewek & Carl Schreck    
PG: pawel.gniewek@berkeley.edu    
CS: carl.schreck@berkeley.edu     
All rights reserved.    
License: BSD-3  

REFERENCES
===============
1. "Jamming by growth", P.G., C.S., O.H., 2018
2. [Self-driven jamming in growing microbial populations](https://www.nature.com/articles/nphys3741), M.D., J.H., C.S., P.G., L.H., S.H., O.H., Nature Phys (2016)
