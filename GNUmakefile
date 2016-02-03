# -*- Mode: Makefile -*- 
OPT=TRUE
MPI=TRUE
DIM=3
USE_64=TRUE
USE_HDF=TRUE
# trace the chain of included makefiles
makefiles += UoB-benchmark

## Define the variables needed by Make.example

# the base name(s) of the application(s) in this directory
ebase = UoB-benchmark

# the location of the Chombo "lib" directory
CHOMBO_HOME = Chombo/lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = AMRElliptic AMRTools BoxTools BaseTools

# the locations of the source code directories
base_dir = .

# input file for 'run' target
INPUT = inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

# application-specific targets
