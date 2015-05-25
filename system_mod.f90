module system_mod

use global_mod, only: rcut,dim,V0,Nmax,dr
use bessel_mod, only: Bessk

implicit none

contains

!-----------------------------------------------------------------------

  function LogPsi(opt,Rm,rij)

    real (kind=8)    :: Rm,rij
    real (kind=8)    :: LogPsi
    integer (kind=4) :: opt
    
    if (opt==0) then
       
       LogPsi = -0.5d0*(Rm/rij)**5
            
    else if (opt==1) then
       
       LogPsi = 2.5d0*(Rm/rij)**5/rij
       
    else if (opt==2) then

       LogPsi = -15.d0*(Rm/rij)**5/rij**2

    else 

       print *, 'The parameter opt in the function LogPsi crash!!!'
       stop
       
    end if

    return
  end function LogPsi

!-----------------------------------------------------------------------

  function Potential(xij,rij)

    implicit none

    real (kind=8) :: Potential,V0,rij,rij6

    real (kind=8), dimension(dim) :: xij

    V0   = 22.0228d0
    rij6 = rij**6

    Potential = V0*(1.d0/rij6-1.d0)/rij6

  end function Potential

!-----------------------------------------------------------------------  

end module system_mod
