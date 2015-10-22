# PATH INTEGRAL GROUND STATE MONTE CARLO

Path integral ground state Monte Carlo (PIGS) is a quantum Monte Carlo method intended to the study of quantum many body systems at zero temperature. This method solves the many body Schrödinger equation of the system by propagating an initial (or trial) wave function in imaginary time making that any non ground state contribution present in the initial guess of the wave function vanishes. The imaginary time propagation is implemented by successive application of the quantum propagator of the system. 

**VPI.F90**

This is the main program that contains the definitions of the variables needed to perform the simulation and the main 
Monte Carlo loop over blocks and steps. It also calls the routines intended to evaluate the desired observables.

**GLOBAL_MOD.F90**

Module that contains the definition of some global variables used along many functions of the code.

**VPI_MOD.F90**

This module contains all the routines concerning the Monte Carlo movements used to sample the path integral and all the
functions needed for the calculation. 

**List of functions**
  
    *GreenFunction
    *OBDMGuess

**List of routines**

  Initialization routines:

    * Jastrow_Table
    * init
    * CheckPoint
  
  Routines for the movement of the particles of the system

    * TranslateChain
    * Staging
    * MoveHead
    * MoveTail
    * TranslateHalfChain
    * StagingHalfChain
    * MoveHeadHalfChain
    * MoveTailHalfChain

  Routines for the evaluation of quantities for Metropolis algorithm

    * UpdateAction
    * UpdateWf
    * UpdatePot

**SAMPLE_MOD.F90**

This module contains all the functions required for the evaluation of the different observable quantities.

**List of functions**
    
    * Var

**List of routines**

  Routines for the evaluation of energies:

    * PotentialEnergy
    * LocalEnergy
    * Accumulate

  Routines for the evaluation of structural quantities

    * PairCorrelation
    * StructureFactor
    * Normalize
    * AccumGr
    * AccumSk
    * NormAvGr
    * NormAvSk

  Routines for the evaluation of OBDM
 
    * OBDM
    * Normalize
    * AccumNr
    * NormAvNr 

**PBC_MOD.F90**

Module that contains functions to implement the periodic boundary conditions in a cubic box (independent of 
dimensionality).

**List of routines**

    * BoundaryConditions
    * MinimumImage

**RANDOM_MOD.F90**

Module that contains the functions and routines used for the generation of random numbers: uniform and gaussian. The 
Twister Mersenne random number generator is used. 

**INTERPOLATE.F90**

Function that performs linear interpolation to find the value of a tabulated function in a given point.

**R8_GAMMA.F90**

Gamma function needed to evaluate the volume of a spherical shell in an arbitrary dimension.

Those are the generic parts of the code, in addition to this modules there are some system-dependent modules:

**SYSTEM_MOD.F90**

Module that contains the trial (or initial) model wave functions of the system (one and two-body terms) and the
potential (also one and two-body).

**BESSEL_MOD.F90**

This is a module containing the definition of Bessel functions and Modified Bessel functions of the first and second 
kind. This is a very specific module needed to build the two-body Jastrow factor of dipoles in two dimensions.

