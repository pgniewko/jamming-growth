# Author : Pawel Gniewek (UC Berkeley)
# Email  : gniewko.pablo@gmail.com 
# License: BSD 3

# FORTRAN compiler
FC       := gfortran

# Relative include and library paths for compilation of the examples
INCLUDE  := -I/usr/local/include
LIB      := -L/usr/lib -L/usr/local/lib
DFLAGS   :=

FTNFLAGS := -O3
LDFLAGS  := $(LIB)
LDLIBS   :=

# Local dirs
SRC      := ./src
BIN      := ./bin
