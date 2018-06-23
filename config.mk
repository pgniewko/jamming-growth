# Author : Pawel Gniewek (UC Berkeley)
# Email  : pawel.gniewek@berkeley.edu
# License: BSD

# FORTRAN compiler
FC       := gfortran

# Relative include and library paths for compilation of the examples
INCLUDE  := -I/usr/local/include
LIB      := -L/usr/lib -L/usr/local/lib
DFLAGS   :=

#FTNFLAGS := -g -fcheck=all -Wall -O3
FTNFLAGS := -O3
#CXXFLAGS := -Wall -O3 $(INCLUDE) $(DBGFLAGS)
LDFLAGS  := $(LIB)
LDLIBS   :=

# Local dirs
SRC      := ./src
BIN      := ./bin
