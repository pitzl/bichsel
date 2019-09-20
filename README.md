
https://www.slac.stanford.edu/~sudong/silicon/

Instructions for Hans Bichsel's CONV/FOLD programs
==================================================

                    Su Dong, Stanford
             Last Update: Apr/13/2016

This package contains the silicon ionization simulation using the Bichsel code.
The details of the Bichsel ionization calculations for silicon are described
in the paper: For the original paper (Rev. Mod. Phys. 60, p663, 1988):
  http://prola.aps.org/abstract/RMP/v60/i3/p663_1

The code is in FORTRAN! developed in the 1980's.

This distribution wrapped the package with user prompt for input particle.
A few plotting examples are added.
http://www.slac.stanford.edu/~sudong/silicon/bichsel_covfold.tar
contains two directories:
Fortran/  this is mostly the original code from Hans Bichsel for generating 
          the single collision cross sections and folding dE/dx spectra.
          The executables were built on SLAC Linux RHEL3.   

plot/    this contains a set of application programs plotting the collsion
         cross section and ionzation spectra. The plotting is implemented
         to use PAW kumac which runs the Hans Bichsel program within it. 
