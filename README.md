# Path integral ground state Monte Carlo

[TOC]

## Introduction

This is fortran 90 code implementing the path integral ground state (PIGS) Monte Carlo method. PIGS is a numerical method to study the properties of the ground state (T=0) of a quantum system. This method is exact (in the statistical sense) for bosons (particles with integer spin). 

PIGS method was formulated by Ceperley[^1] and called variational path integral (VPI), a method closely related to the finite temperature path integral Monte Carlo method (PIMC). Later, in 2000 Sarsa et al. [^2] performed a first realization of the method with applications to several test systems (Quantum harmonic oscillator and liquid Helium-4) showing that the PIGS method was a valid alternative to other ground state methods like for example diffusion Monte Carlo (DMC). After this work the method have gathered popularity in the toolbox of the quantum condensed matter physicists and have shown its value as an alternative to other methods due to the feature of being able to evaluate in a pure way (with no variational bias) observables that do not commute with the hamiltonian[^3] and to the fact that it can provide correct results independently of the trial wave function that is used[^4] (if the wave function has the correct bosonic simmetry). 

With the introduction of the Worm Algorithm[^5] [^6] (WA) PIGS method was able to evaluate in an exact way, and with the correct normalization, the one body density matrix (OBDM) of a many-body system which gives access to the study of the Bose-Einstein condensation phenomena in this kind of systems. 

## Theoretical background

PIGS is an imaginary time propagation (or projector) method, projector methods are a family of Monte Carlo methods that are able to solve the many-body Schrödinger equation by repeatedly applying a projection operator that wipes out the non-ground-state components of an initial (or trial) wave function. The result of this process is the ground-state wave function of the system. 

The starting point of all projector methods is the imaginary time dependent Schrödinger equation:

$-\frac{\partial \psi}{\partial \tau} = H \psi $

where H is the hamiltonian of the system. The formal solution of this equation is given by:

$ \psi(\tau) = e^{-H\tau}\psi(0)$

In the limit of $\tau \rightarrow \infty$, $\psi(\tau \rightarrow \infty)$ equals the ground state wave function $\phi_0$. The operator $e^{-H\tau}$ is the representation independent Green's function of the system and is, in general, unknown but some suitable short-time approximations can be found[^7] . 

 ## References

[^1]: D. Ceperley, Rev. Mod. Phys. 67, 2, (1995)
[^2]: A. Sarsa, K. E. Schmidt and W. R. Magro, J. Chem. Phys. **113**, 1366 (2000) 
[^3]: J. E. Cuervo, P-N. Roy and M. Boninsegni, J. Chem. Phys. **122**, 114504 (2005)
[^4]: M. Rossi, M. Nava, L. Reatto and D. E. Galli, J. Chem. Phys. **131**, 154108 (2009)
[^5]: M. Boninsegni, N. Prokof'ev and B. Svistunov, Phys. Rev. Lett. **96**, 070601 (2006)
[^6]: M. Boninsegni, N. Prokof'ev and B. Svistunov, Phys. Rev. E **74**, 036701 (2006)
[^7]: S. A. Chin and C. R. Chen, J. Chem. Phys. **114**, 1409, (2002) 

