# block-structured-benchmark
A stress test for modest clusters based on the Chombo block-structured AMR toolkit


Prerequisites (on GNU/Linux): C++ and Fortran compilers, perl, csh, mpi and hdf5
Optional : openmp

Compiling:

1. Get Chombo 3.2 and place it in a Chombo subdirectory, e.g
> svn export https://anag-repo.lbl.gov/svn/Chombo/release/3.2 Chombo

2. Edit Make.defs.local. On some systems it may be better to start from on of the examples in Chombo/lib/mk/local/. 
For example,  Chombo/lib/mk/local/Make.defs.Archer is known to work on the UK Cray XC30 machine ARCHER

3: build 
> make all 

to produce an executable


Running the benchmark:

 

 
