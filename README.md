# block-structured-benchmark
A stress test for modest clusters based on the Chombo block-structured AMR toolkit


Prerequisites (on GNU/Linux): C++ and Fortran compilers, perl, csh, mpi and hdf5
Optional : openmp

Compiling:

1. Get Chombo 3.2 and place it in a Chombo subdirectory, e.g
> svn export https://anag-repo.lbl.gov/svn/Chombo/release/3.2 Chombo

2. Edit Chombo/lib/mk/local/Make.defs.local, based on Chombo/lib/mk/Make.defs.local.template. 
   CXX, FC, MPICXX, HDFLIBFLAGS, HDFMPILIBFLAGS,HDFINCFLAGS,HDFMPIINCFLAGS must be set.
   Optimization settings can be specified in cxxoptflags, foptflags

   On some systems it may be better to start from on of the examples in Chombo/lib/mk/local/. 
   For example,  Chombo/lib/mk/local/Make.defs.archer is known to work on the UK Cray XC30 machine ARCHER

3: build 
> make all 

to produce an executable named something like UoB-benchmark3d.Linux.64.$MPICXX.$FC.DEBUG.OPT.MPI.ex
(E.g, when mpicxx is mpic++ and FC is gfortran, you have UoB-benchmark3d.Linux.64.mpic++.gfortran.DEBUG.OPT.MPI.ex)


Running the benchmark:

Run the exectuable without command line arguments over a range of node counts, with one (and only one) mpi rank per core,
up to ~500 cores. Keep the files named *.pout.0*, *.timer* and *hdf5




 

 
